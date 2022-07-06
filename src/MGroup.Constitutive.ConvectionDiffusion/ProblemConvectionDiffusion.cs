using System.Collections.Generic;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Solution;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.Constitutive.ConvectionDiffusion.Providers;
using System.Linq;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;


namespace MGroup.Constitutive.ConvectionDiffusion
{
	public class ProblemConvectionDiffusion : IAlgebraicModelInterpreter, ITransientAnalysisProvider, INonTransientAnalysisProvider, INonLinearProvider
	{
		private IGlobalMatrix convection, diffusion, production, firstTimeDerivativeMatrix;
		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ISolver solver;
		private ElementProductionProvider productionProvider = new ElementProductionProvider();
		private ElementDiffusionProvider diffusionProvider = new ElementDiffusionProvider();
		private ElementConvectionProvider convectionProvider = new ElementConvectionProvider();
		private ElementFirstTimeDerivativeMatrixProvider fistTimeDerivativeMatrixProvider = new ElementFirstTimeDerivativeMatrixProvider();

		// TODO: Is this right?
		private readonly IElementMatrixPredicate rebuildDiffusionPredicate = new MaterialModifiedElementMarixPredicate();

		public ProblemConvectionDiffusion(IModel model, IAlgebraicModel algebraicModel, ISolver solver)
		{
			this.model = model;
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			algebraicModel.BoundaryConditionsInterpreter = this;

			ActiveDofs.AddDof(ConvectionDiffusionDof.UnknownVariable);
		}

		private IGlobalMatrix Convection
		{
			get
			{
				if (convection == null) BuildConvection();
				return convection;
			}
		}

		private IGlobalMatrix Diffusion
		{
			get
			{
				if (diffusion == null) BuildDiffusion();
				return diffusion;
			}
		}

		private IGlobalMatrix Production
		{
			get
			{
				if (production == null) BuildProduction();
				return production;
			}
		}

		private IGlobalMatrix FirstTimeDerivativeMatrix
		{
			get
			{
				if (firstTimeDerivativeMatrix == null) BuildFirstOrderDerivativeMatrix();
				return firstTimeDerivativeMatrix;
			}
		}

		public ActiveDofs ActiveDofs { get; } = new ActiveDofs();

		private void BuildConvection() => convection = algebraicModel.BuildGlobalMatrix(convectionProvider);

		private void BuildDiffusion() => convection = algebraicModel.BuildGlobalMatrix(diffusionProvider);

		private void BuildProduction() => convection = algebraicModel.BuildGlobalMatrix(productionProvider);

		private void BuildFirstOrderDerivativeMatrix() => convection = algebraicModel.BuildGlobalMatrix(fistTimeDerivativeMatrixProvider);

		// TODO:Is this right?
		private void RebuildDiffusion() => algebraicModel.RebuildGlobalMatrixPartially(diffusion, 
			model.EnumerateElements, diffusionProvider, rebuildDiffusionPredicate);

		#region IAnalyzerProvider Members
		public void ClearMatrices()
		{
			convection = null;
			diffusion = null;
			production = null;
			firstTimeDerivativeMatrix = null;
		}

		public void Reset()
		{
			ClearMatrices();
		}

		public void GetProblemDofTypes()
		{
			// Contents were commented out, function no longer used
		}
		#endregion

		#region IImplicitIntegrationProvider Members

		public void LinearCombinationOfMatricesIntoEffectiveMatrix(TransientAnalysisCoefficients coefficients)
		{
			IGlobalMatrix matrix = Diffusion;
			matrix.AddIntoThis(Convection);
			matrix.AddIntoThis(Production);
			matrix.AxpyIntoThis(FirstTimeDerivativeMatrix, coefficients.FirstOrderDerivativeCoefficient);
			solver.LinearSystem.Matrix = matrix;
		}

		public void LinearCombinationOfMatricesIntoEffectiveMatrixNoOverwrite(TransientAnalysisCoefficients coefficients)
		{
			IGlobalMatrix matrix = Diffusion.Copy();
			matrix.AddIntoThis(Convection);
			matrix.AddIntoThis(Production);
			matrix.AxpyIntoThis(FirstTimeDerivativeMatrix, coefficients.FirstOrderDerivativeCoefficient);
			solver.LinearSystem.Matrix = matrix;
		}

		public void ProcessRhs(TransientAnalysisCoefficients coefficients, IGlobalVector rhs)
		{
			// Method intentionally left empty.
		}

		public IGlobalVector GetSecondOrderDerivativeVectorFromBoundaryConditions(double time) => algebraicModel.CreateZeroVector();

		public IGlobalVector GetFirstOrderDerivativeVectorFromBoundaryConditions(double time)
		{
			var boundaryConditions = model.EnumerateBoundaryConditions(model.EnumerateSubdomains().First().ID).ToArray();
			foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IConvectionDiffusionDofType>>())
			{
				boundaryCondition.CurrentTime = time;
			}

			IGlobalVector velocities = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(
			id =>
			{
				var boundaryConditions = model.EnumerateBoundaryConditions(id);
				foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IConvectionDiffusionDofType>>().ToArray())
				{
					boundaryCondition.CurrentTime = time;
				}

				// TODO: Is this right?
				return boundaryConditions
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalFirstTimeDerivativeBoundaryCondition>();
			},
			velocities);

			// TODO: Is this right?
			algebraicModel.AddToGlobalVector(boundaryConditions
				.SelectMany(x => x.EnumerateDomainBoundaryConditions())
				.OfType<IDomainFirstTimeDerivativeBoundaryCondition>(),
			velocities);

			return velocities;
		}

		public IGlobalVector GetRhs(double time)
		{
			solver.LinearSystem.RhsVector.Clear();
			AssignRhs();
			IGlobalVector result = solver.LinearSystem.RhsVector.Copy();

			return result;
		}

		public IGlobalVector SecondOrderDerivativeMatrixVectorProduct(IGlobalVector vector) => algebraicModel.CreateZeroVector();

		public IGlobalVector FirstOrderDerivativeMatrixVectorProduct(IGlobalVector vector)
		{// TODO: is the matrix right?
			IGlobalVector result = algebraicModel.CreateZeroVector();
			FirstTimeDerivativeMatrix.MultiplyVector(vector, result);
			return result;
		}

		#endregion

		#region IStaticProvider Members

		public void CalculateMatrix()
		{// TODO: is the matrix right?
			if (firstTimeDerivativeMatrix == null) BuildFirstOrderDerivativeMatrix();
			solver.LinearSystem.Matrix = firstTimeDerivativeMatrix;
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
					.OfType<INodalUnknownVariableFluxBoundaryCondition>()
					.Where(x => model.EnumerateBoundaryConditions(id)
						.SelectMany(x => x.EnumerateNodalBoundaryConditions())
						.OfType<INodalUnknownVariableBoundaryCondition>()
						.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false), 
					solver.LinearSystem.RhsVector);
		}

		public IEnumerable<INodalNeumannBoundaryCondition<IDofType>> EnumerateEquivalentNeumannBoundaryConditions(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateEquivalentNodalNeumannBoundaryConditions(model.EnumerateElements(subdomainID)))
				.OfType<INodalUnknownVariableFluxBoundaryCondition>()
				.Where(x => model.EnumerateBoundaryConditions(subdomainID)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalUnknownVariableBoundaryCondition>()
					.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false);

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering() =>
			model.EnumerateSubdomains()
				.SelectMany(x => model.EnumerateBoundaryConditions(x.ID)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions()).OfType<INodalUnknownVariableBoundaryCondition>()
					.OrderBy(x => x.Node.ID)
					.GroupBy(x => (x.Node.ID, x.DOF))
					.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount))))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateNodalBoundaryConditions()).OfType<INodalUnknownVariableBoundaryCondition>()
				.OrderBy(x => x.Node.ID)
				.GroupBy(x => (x.Node.ID, x.DOF))
				.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount)))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));
	}
}
