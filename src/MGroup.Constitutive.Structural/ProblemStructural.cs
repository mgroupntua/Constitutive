using System.Collections.Generic;
using MGroup.Constitutive.Structural.Providers;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.Solution.AlgebraicModel;
using System.Linq;
using MGroup.Constitutive.Structural.BoundaryConditions;

//TODO: Usually the LinearSystem is passed in, but for GetRHSFromHistoryLoad() it is stored as a field. Decide on one method.
//TODO: Right now this class decides when to build or rebuild the matrices. The analyzer should decide that.
namespace MGroup.Constitutive.Structural
{
	public class ProblemStructural : IAlgebraicModelInterpreter, ITransientAnalysisProvider, INonTransientAnalysisProvider, INonLinearProvider
	{
		private IGlobalMatrix mass, damping, stiffness;
		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ISolver solver;
		private ElementStructuralStiffnessProvider stiffnessProvider = new ElementStructuralStiffnessProvider();
		private ElementStructuralMassProvider massProvider = new ElementStructuralMassProvider();
		private ElementStructuralDampingProvider dampingProvider = new ElementStructuralDampingProvider();
		private readonly IElementMatrixPredicate rebuildStiffnessPredicate = new MaterialModifiedElementMarixPredicate();

		public ProblemStructural(IModel model, IAlgebraicModel algebraicModel, ISolver solver)
		{
			this.model = model;
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			this.algebraicModel.BoundaryConditionsInterpreter = this;

			ActiveDofs.AddDof(StructuralDof.TranslationX);
			ActiveDofs.AddDof(StructuralDof.TranslationY);
			ActiveDofs.AddDof(StructuralDof.TranslationZ);
			ActiveDofs.AddDof(StructuralDof.RotationX);
			ActiveDofs.AddDof(StructuralDof.RotationY);
			ActiveDofs.AddDof(StructuralDof.RotationZ);
		}

		private IGlobalMatrix Mass
		{
			get
			{
				if (mass == null) BuildMass();
				return mass;
			}
		}

		private IGlobalMatrix Damping
		{
			get
			{
				if (damping == null) BuildDamping();
				return damping;
			}
		}

		private IGlobalMatrix Stiffness
		{
			get
			{
				if (stiffness == null)
					BuildStiffness();
				else
				{
					RebuildStiffness(); // This is the same but also resets the material modified properties. 
				}
				return stiffness;
			}
		}

		public ActiveDofs ActiveDofs { get; } = new ActiveDofs();

		private void BuildStiffness() => stiffness = algebraicModel.BuildGlobalMatrix(stiffnessProvider);

		private void RebuildStiffness()
		{
			algebraicModel.RebuildGlobalMatrixPartially(stiffness, model.EnumerateElements, stiffnessProvider, rebuildStiffnessPredicate);

			// Original code kept, in case we need to reproduce its behavior.
			//foreach (ISubdomain subdomain in model.Subdomains)
			//{
			//    if (subdomain.MaterialsModified)
			//    {
			//        stiffness[subdomain.ID] = solver.BuildGlobalMatrix(subdomain, stiffnessProvider);
			//        subdomain.ResetMaterialsModifiedProperty();
			//    }
			//}
		}

		private void BuildMass() => mass = algebraicModel.BuildGlobalMatrix(massProvider);

		//TODO: With Rayleigh damping, C is more efficiently built using linear combinations of global K, M, 
		//      instead of building and assembling element k, m matrices.
		private void BuildDamping() => damping = algebraicModel.BuildGlobalMatrix(dampingProvider);

		#region IAnalyzerProvider Members
		public void ClearMatrices()
		{
			damping = null;
			stiffness = null;
			mass = null;
		}

		public void Reset()
		{
			// TODO: Check if we should clear material state - (goat) removed that, seemed erroneous
			//foreach (ISubdomain subdomain in model.Subdomains)
			//	foreach (IElement element in subdomain.Elements)
			//		element.ElementType.ClearMaterialState();

			damping = null;
			stiffness = null;
			mass = null;
		}

		public void GetProblemDofTypes()
		{
			//model.AllDofs.AddDof(StructuralDof.TranslationX);
			//model.AllDofs.AddDof(StructuralDof.TranslationY);
			//model.AllDofs.AddDof(StructuralDof.TranslationZ);
			//model.AllDofs.AddDof(StructuralDof.RotationX);
			//model.AllDofs.AddDof(StructuralDof.RotationY);
			//model.AllDofs.AddDof(StructuralDof.RotationZ);
		}
		#endregion

		#region IImplicitIntegrationProvider Members

		public void LinearCombinationOfMatricesIntoEffectiveMatrix(TransientAnalysisCoefficients coefficients)
		{
			//TODO: when the matrix is mutated, the solver must be informed via observers (or just flags).
			IGlobalMatrix matrix = Stiffness;
			matrix.LinearCombinationIntoThis(coefficients.ZeroOrderDerivativeCoefficient, Mass, coefficients.SecondOrderDerivativeCoefficient);
			matrix.AxpyIntoThis(Damping, coefficients.FirstOrderDerivativeCoefficient);
			solver.LinearSystem.Matrix = matrix;
		}

		public void LinearCombinationOfMatricesIntoEffectiveMatrixNoOverwrite(TransientAnalysisCoefficients coefficients)
		{
			//TODO: when the matrix is mutated, the solver must be informed via observers (or just flags).
			IGlobalMatrix matrix = Stiffness.Copy();
			matrix.LinearCombinationIntoThis(coefficients.ZeroOrderDerivativeCoefficient, Mass, coefficients.SecondOrderDerivativeCoefficient);
			matrix.AxpyIntoThis(Damping, coefficients.FirstOrderDerivativeCoefficient);
			solver.LinearSystem.Matrix = matrix;
		}

		public void ProcessRhs(TransientAnalysisCoefficients coefficients, IGlobalVector rhs)
		{
			// Method intentionally left empty.
		}
		//Transient
		public IGlobalVector GetSecondOrderDerivativeVectorFromBoundaryConditions(double time)
		{
			var boundaryConditions = model.EnumerateBoundaryConditions(model.EnumerateSubdomains().First().ID).ToArray();
			foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IStructuralDofType>>())
			{
				boundaryCondition.CurrentTime = time;
			}

			IGlobalVector accelerations = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(id =>
			{
				var boundaryConditions = model.EnumerateBoundaryConditions(id).ToArray();
				foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IStructuralDofType>>())
				{
					boundaryCondition.CurrentTime = time;
				}

				return boundaryConditions
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalAccelerationBoundaryCondition>();
			},
				accelerations);
			algebraicModel.AddToGlobalVector(boundaryConditions
					.SelectMany(x => x.EnumerateDomainBoundaryConditions())
					.OfType<IDomainAccelerationBoundaryCondition>(),
				accelerations);

			return accelerations;
		}

		public IGlobalVector GetFirstOrderDerivativeVectorFromBoundaryConditions(double time)
		{
			var boundaryConditions = model.EnumerateBoundaryConditions(model.EnumerateSubdomains().First().ID).ToArray();
			foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IStructuralDofType>>())
			{
				boundaryCondition.CurrentTime = time;
			}

			IGlobalVector velocities = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(id =>
			{
				var boundaryConditions = model.EnumerateBoundaryConditions(id);
				foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IStructuralDofType>>().ToArray())
				{
					boundaryCondition.CurrentTime = time;
				}

				return boundaryConditions
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalVelocityBoundaryCondition>();
			},
				velocities);
			algebraicModel.AddToGlobalVector(boundaryConditions
					.SelectMany(x => x.EnumerateDomainBoundaryConditions())
					.OfType<IDomainVelocityBoundaryCondition>(),
				velocities);

			return velocities;
		}

		public IGlobalVector GetRhs(double time)
		{
			solver.LinearSystem.RhsVector.Clear(); //TODO: this is also done by model.AssignLoads()
			AssignRhs();
			algebraicModel.AddToGlobalVector(EnumerateEquivalentNeumannBoundaryConditions, solver.LinearSystem.RhsVector);

			return solver.LinearSystem.RhsVector.Copy();
		}

		public IGlobalVector SecondOrderDerivativeMatrixVectorProduct(IGlobalVector vector)
		{
			IGlobalVector result = algebraicModel.CreateZeroVector();
			Mass.MultiplyVector(vector, result);
			return result;
		}
		//
		public IGlobalVector FirstOrderDerivativeMatrixVectorProduct(IGlobalVector vector)
		{
			IGlobalVector result = algebraicModel.CreateZeroVector();
			Damping.MultiplyVector(vector, result);
			return result;
		}


		public IEnumerable<INodalNeumannBoundaryCondition<IDofType>> EnumerateEquivalentNeumannBoundaryConditions(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateEquivalentNodalNeumannBoundaryConditions(model.EnumerateElements(subdomainID)))
				.OfType<INodalLoadBoundaryCondition>()
				.Where(x => model.EnumerateBoundaryConditions(subdomainID)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalDisplacementBoundaryCondition>()
					.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false);

		#endregion

		#region IStaticProvider Members
		public void CalculateMatrix()
		{
			if (stiffness == null) BuildStiffness();
			solver.LinearSystem.Matrix = stiffness;
		}
		#endregion

		#region INonLinearProvider Members
		public double CalculateRhsNorm(IGlobalVector rhs) => rhs.Norm2();

		public void ProcessInternalRhs(IGlobalVector solution, IGlobalVector rhs) { }
		#endregion

		public void AssignRhs()
		{
			solver.LinearSystem.RhsVector.Clear();
			algebraicModel.AddToGlobalVector(id =>
				model.EnumerateBoundaryConditions(id)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalLoadBoundaryCondition>()
					.Where(x => model.EnumerateBoundaryConditions(id)
						.SelectMany(x => x.EnumerateNodalBoundaryConditions())
						.OfType<INodalDisplacementBoundaryCondition>()
						.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false),
				solver.LinearSystem.RhsVector);
		}

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering() =>
			model.EnumerateSubdomains()
				.SelectMany(x => model.EnumerateBoundaryConditions(x.ID)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions()).OfType<INodalDisplacementBoundaryCondition>()
					.OrderBy(x => x.Node.ID)
					.GroupBy(x => (x.Node.ID, x.DOF))
					.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount))))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateNodalBoundaryConditions()).OfType<INodalDisplacementBoundaryCondition>()
				.OrderBy(x => x.Node.ID)
				.GroupBy(x => (x.Node.ID, x.DOF))
				.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount)))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));
	}
}
