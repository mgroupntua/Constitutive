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
using MGroup.Constitutive.Thermal.Providers;
using System.Linq;
using MGroup.Constitutive.Thermal.BoundaryConditions;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.Constitutive.Thermal.InitialConditions;

//TODO: Usually the LinearSystem is passed in, but for GetRHSFromHistoryLoad() it is stored as a field. Decide on one method.
//TODO: I am not too fond of the provider storing global sized matrices.
namespace MGroup.Constitutive.Thermal
{
	public class ProblemThermal : IAlgebraicModelInterpreter, ITransientAnalysisProvider, INonTransientAnalysisProvider, INonLinearProvider
	{
		private bool shouldRebuildConductivityMatrixForZeroOrderDerivativeMatrixVectorProduct = false;
		private IGlobalMatrix capacity, conductivity;
		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ISolver solver;
		private ElementConductivityProvider conductivityProvider = new ElementConductivityProvider();
		private ElementCapacityProvider capacityProvider = new ElementCapacityProvider();
		private readonly IElementMatrixPredicate rebuildConductivityPredicate = new MaterialModifiedElementMarixPredicate();

		public ProblemThermal(IModel model, IAlgebraicModel algebraicModel, ISolver solver)
		{
			this.model = model;
			this.algebraicModel = algebraicModel;
			this.solver = solver;
			algebraicModel.BoundaryConditionsInterpreter = this;

			ActiveDofs.AddDof(ThermalDof.Temperature);
		}

		private IGlobalMatrix Capacity
		{
			get
			{
				if (capacity == null) BuildCapacity();
				return capacity;
			}
		}

		private IGlobalMatrix Conductivity
		{
			get
			{
				if (conductivity == null) BuildConductivity();
				//else RebuildConductivityMatrices();
				return conductivity;
			}
		}

		public ActiveDofs ActiveDofs { get; } = new ActiveDofs();

		private void BuildConductivity() => conductivity = algebraicModel.BuildGlobalMatrix(conductivityProvider);
		
		private void RebuildConductivity() => algebraicModel.RebuildGlobalMatrixPartially(conductivity, 
			model.EnumerateElements, conductivityProvider, rebuildConductivityPredicate);

		private void BuildCapacity() => capacity = algebraicModel.BuildGlobalMatrix(capacityProvider);

		#region IAnalyzerProvider Members
		public void ClearMatrices()
		{
			capacity = null;
			conductivity = null;
		}

		public void Reset()
		{
			// TODO: Check if we should clear material state - (goat) removed that, seemed erroneous
			//foreach (ISubdomain subdomain in model.Subdomains)
			//	foreach (IElement element in subdomain.Elements)
			//		element.ElementType.ClearMaterialState();

			conductivity = null;
			capacity = null;
		}

		public void GetProblemDofTypes()
		{
			//model.AllDofs.AddDof(ThermalDof.Temperature);
		}
		#endregion

		#region IImplicitIntegrationProvider Members

		public void LinearCombinationOfMatricesIntoEffectiveMatrix(TransientAnalysisCoefficients coefficients)
		{
			//TODO: when the matrix is mutated, the solver must be informed via observers (or just flags).
			IGlobalMatrix matrix = Conductivity;
			if (coefficients[DifferentiationOrder.First] != 0)
			{
				matrix.LinearCombinationIntoThis(coefficients[DifferentiationOrder.Zero], Capacity, coefficients[DifferentiationOrder.First]);
			}
			else
			{
				matrix.ScaleIntoThis(coefficients[DifferentiationOrder.Zero]);
			}

			solver.LinearSystem.Matrix = matrix;
			shouldRebuildConductivityMatrixForZeroOrderDerivativeMatrixVectorProduct = true;
		}

		public void LinearCombinationOfMatricesIntoEffectiveMatrixNoOverwrite(TransientAnalysisCoefficients coefficients)
		{
			IGlobalMatrix matrix = Conductivity.Copy();
			matrix.AxpyIntoThis(Capacity, coefficients[DifferentiationOrder.First]);
			solver.LinearSystem.Matrix = matrix;
		}

		public void ProcessRhs(TransientAnalysisCoefficients coefficients, IGlobalVector rhs)
		{
			// Method intentionally left empty.
		}

		private IGlobalVector GetSecondOrderDerivativeVectorFromBoundaryConditions(double time) => algebraicModel.CreateZeroVector();

		private IGlobalVector GetFirstOrderDerivativeVectorFromBoundaryConditions(double time)
		{
			var boundaryConditions = model.EnumerateBoundaryConditions(model.EnumerateSubdomains().First().ID).ToArray();
			foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IThermalDofType>>())
			{
				boundaryCondition.CurrentTime = time;
			}

			IGlobalVector velocities = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(id =>
			{
				var boundaryConditions = model.EnumerateBoundaryConditions(id).ToArray();
				foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IThermalDofType>>())
				{
					boundaryCondition.CurrentTime = time;
				}

				return boundaryConditions
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalCapacityBoundaryCondition>();
			},
				velocities);
			algebraicModel.AddToGlobalVector(boundaryConditions
					.SelectMany(x => x.EnumerateDomainBoundaryConditions())
					.OfType<IDomainCapacityBoundaryCondition>(),
				velocities);

			return velocities;
		}

		private IGlobalVector GetZeroOrderDerivativeVectorFromInitialConditions(double time)
		{
			IGlobalVector temperatures = algebraicModel.CreateZeroVector();

			if (time == 0)
			{
				algebraicModel.AddToGlobalVector(
				id =>
				{
					var initialConditions = model.EnumerateInitialConditions(id).ToArray();
					return initialConditions
						.SelectMany(x => x.EnumerateNodalInitialConditions())
						.OfType<INodalTemperatureInitialCondition>();
				},
				temperatures);

				algebraicModel.AddToGlobalVector(model.EnumerateInitialConditions(model.EnumerateSubdomains().First().ID)
					.SelectMany(x => x.EnumerateDomainInitialConditions())
					.OfType<IDomainTemperatureInitialCondition>(),
				temperatures);
			}

			return temperatures;
		}

		public IGlobalVector GetVectorFromModelConditions(DifferentiationOrder differentiationOrder, double time) => differentiationOrder switch
		{
			DifferentiationOrder.Zero => GetZeroOrderDerivativeVectorFromInitialConditions(time),
			DifferentiationOrder.First => GetFirstOrderDerivativeVectorFromBoundaryConditions(time),
			_ => algebraicModel.CreateZeroVector(),
		};

		public IGlobalVector GetRhs(double time)
		{
			solver.LinearSystem.RhsVector.Clear(); //TODO: this is also done by model.AssignLoads()
			AssignRhs();

			algebraicModel.AddToGlobalVector(id =>
			{
				var transientBCs = model.EnumerateBoundaryConditions(id).OfType<ITransientBoundaryConditionSet<IThermalDofType>>().ToArray();
				foreach (var boundaryCondition in transientBCs)
				{
					boundaryCondition.CurrentTime = time;
				}

				return transientBCs
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalHeatFluxBoundaryCondition>()
					.Where(x => model.EnumerateBoundaryConditions(id)
						.SelectMany(x => x.EnumerateNodalBoundaryConditions())
						.OfType<INodalTemperatureBoundaryCondition>()
						.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false);
			},
				solver.LinearSystem.RhsVector);

			IGlobalVector result = solver.LinearSystem.RhsVector.Copy();
			return result;
		}

		private IGlobalVector SecondOrderDerivativeMatrixVectorProduct(IGlobalVector vector) => algebraicModel.CreateZeroVector();

		private IGlobalVector FirstOrderDerivativeMatrixVectorProduct(IGlobalVector vector)
		{
			IGlobalVector result = algebraicModel.CreateZeroVector();
			Capacity.MultiplyVector(vector, result);
			return result;
		}

		private IGlobalVector ZeroOrderDerivativeMatrixVectorProduct(IGlobalVector vector)
		{
			if (shouldRebuildConductivityMatrixForZeroOrderDerivativeMatrixVectorProduct)
			{
				BuildConductivity();
				shouldRebuildConductivityMatrixForZeroOrderDerivativeMatrixVectorProduct = false;
			}

			IGlobalVector result = algebraicModel.CreateZeroVector();
			Conductivity.MultiplyVector(vector, result);
			return result;
		}

		public IGlobalVector MatrixVectorProduct(DifferentiationOrder differentiationOrder, IGlobalVector vector) => differentiationOrder switch
		{
			DifferentiationOrder.Zero => ZeroOrderDerivativeMatrixVectorProduct(vector),
			DifferentiationOrder.First => FirstOrderDerivativeMatrixVectorProduct(vector),
			_ => algebraicModel.CreateZeroVector(),
		};

		#endregion

		#region IStaticProvider Members

		public void CalculateMatrix()
		{
			if (conductivity == null) BuildConductivity();
			solver.LinearSystem.Matrix = conductivity;
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
					.OfType<INodalHeatFluxBoundaryCondition>()
					.Where(x => model.EnumerateBoundaryConditions(id)
						.SelectMany(x => x.EnumerateNodalBoundaryConditions())
						.OfType<INodalTemperatureBoundaryCondition>()
						.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false), 
					solver.LinearSystem.RhsVector);
		}

		public IEnumerable<INodalNeumannBoundaryCondition<IDofType>> EnumerateEquivalentNeumannBoundaryConditions(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateEquivalentNodalNeumannBoundaryConditions(model.EnumerateElements(subdomainID)))
				.OfType<INodalHeatFluxBoundaryCondition>()
				.Where(x => model.EnumerateBoundaryConditions(subdomainID)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalTemperatureBoundaryCondition>()
					.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false);

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering() =>
			model.EnumerateSubdomains()
				.SelectMany(x => model.EnumerateBoundaryConditions(x.ID)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions()).OfType<INodalTemperatureBoundaryCondition>()
					.OrderBy(x => x.Node.ID)
					.GroupBy(x => (x.Node.ID, x.DOF))
					.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount))))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateNodalBoundaryConditions()).OfType<INodalTemperatureBoundaryCondition>()
				.OrderBy(x => x.Node.ID)
				.GroupBy(x => (x.Node.ID, x.DOF))
				.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount)))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));
	}
}
