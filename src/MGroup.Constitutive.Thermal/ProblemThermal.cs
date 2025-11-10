using System;
using System.Collections.Generic;
using System.Linq;

using MGroup.Constitutive.Thermal.BoundaryConditions;
using MGroup.Constitutive.Thermal.InitialConditions;
using MGroup.Constitutive.Thermal.Providers;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution.AlgebraicModel;

//TODO: I am not too fond of the provider storing global sized matrices.
namespace MGroup.Constitutive.Thermal
{
	public class ProblemThermal : IAlgebraicModelInterpreter, ITransientAnalysisProvider, INonTransientAnalysisProvider, INonLinearProvider
	{
		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ElementConductivityProvider conductivityProvider = new ElementConductivityProvider();
		private readonly ElementCapacityProvider capacityProvider = new ElementCapacityProvider();
		private readonly IElementMatrixPredicate rebuildConductivityPredicate = new MaterialModifiedElementMarixPredicate();
		private IMatrix capacity, conductivity;
		private TransientAnalysisPhase analysisPhase = TransientAnalysisPhase.SteadyStateSolution;

		public ProblemThermal(IModel model, IAlgebraicModel algebraicModel)
		{
			this.model = model;
			this.algebraicModel = algebraicModel;
			algebraicModel.BoundaryConditionsInterpreter = this;

			ActiveDofs.AddDof(ThermalDof.Temperature);
		}

		private IMatrix Capacity
		{
			get
			{
				if (capacity == null) BuildCapacity();
				return capacity;
			}
		}

		private IMatrix Conductivity
		{
			get
			{
				if (conductivity == null) BuildConductivity();
				//else RebuildConductivityMatrices();
				return conductivity;
			}
		}

		public ActiveDofs ActiveDofs { get; } = new ActiveDofs();

		public void SetTransientAnalysisPhase(TransientAnalysisPhase phase) => analysisPhase = phase;

		public DifferentiationOrder ProblemOrder => DifferentiationOrder.First;

		private void BuildConductivity() => conductivity = algebraicModel.BuildGlobalMatrix(conductivityProvider);
		
		private void RebuildConductivity() => algebraicModel.RebuildGlobalMatrixPartially(conductivity, 
			model.EnumerateElements, conductivityProvider, rebuildConductivityPredicate);

		private void BuildCapacity() => capacity = algebraicModel.BuildGlobalMatrix(capacityProvider);

		public void Reset()
		{
			conductivity = null;
			capacity = null;
		}

		public IMatrix GetMatrix(DifferentiationOrder differentiationOrder) => differentiationOrder switch
		{
			DifferentiationOrder.Zero => Conductivity,
			DifferentiationOrder.First => Capacity,
			_ => algebraicModel.CreateEmptyMatrix(),
		};

		private IVector GetFirstOrderDerivativeVectorFromBoundaryConditions(double time)
		{
			var boundaryConditions = model.EnumerateBoundaryConditions(model.EnumerateSubdomains().First().ID).ToArray();
			foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IThermalDofType>>())
			{
				boundaryCondition.CurrentTime = time;
			}

			IVector temperatureDerivatives = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(id =>
			{
				var boundaryConditions = model.EnumerateBoundaryConditions(id).ToArray();
				foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IThermalDofType>>())
				{
					boundaryCondition.CurrentTime = time;
				}

				return boundaryConditions
					.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id)))
					.OfType<INodalCapacityBoundaryCondition>();
			},
				temperatureDerivatives);

			//algebraicModel.AddToGlobalVector(boundaryConditions
			//		.SelectMany(x => x.EnumerateDomainBoundaryConditions())
			//		.OfType<IDomainCapacityBoundaryCondition>(),
			//	temperatureDerivatives);

			return temperatureDerivatives;
		}

		private IVector GetZeroOrderDerivativeVectorFromInitialConditions(double time)
		{
			IVector temperatures = algebraicModel.CreateZeroVector();

			if (time == 0)
			{
				algebraicModel.AddToGlobalVector(
				id =>
				{
					var initialConditions = model.EnumerateInitialConditions(id).ToArray();
					return initialConditions
						.SelectMany(x => x.EnumerateNodalInitialConditions(model.EnumerateElements(id)))
						.OfType<INodalTemperatureInitialCondition>();
				},
				temperatures);

				//algebraicModel.AddToGlobalVector(model.EnumerateInitialConditions(model.EnumerateSubdomains().First().ID)
				//	.SelectMany(x => x.EnumerateDomainInitialConditions())
				//	.OfType<IDomainTemperatureInitialCondition>(),
				//temperatures);
			}

			return temperatures;
		}

		public IVector GetVectorFromModelConditions(DifferentiationOrder differentiationOrder, double time) => differentiationOrder switch
		{
			DifferentiationOrder.Zero => GetZeroOrderDerivativeVectorFromInitialConditions(time),
			DifferentiationOrder.First => GetFirstOrderDerivativeVectorFromBoundaryConditions(time),
			_ => algebraicModel.CreateZeroVector(),
		};

		public IVector GetRhs(double time)
		{
			IVector rhs = algebraicModel.CreateZeroVector();

			algebraicModel.AddToGlobalVector(id =>
				model.EnumerateBoundaryConditions(id)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id)))
					.OfType<INodalHeatFluxBoundaryCondition>()
					.Where(x => model.EnumerateBoundaryConditions(id)
						.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id)))
						.OfType<INodalTemperatureBoundaryCondition>()
						.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false),
					rhs);

			algebraicModel.AddToGlobalVector(id =>
			{
				var transientBCs = model.EnumerateBoundaryConditions(id).OfType<ITransientBoundaryConditionSet<IThermalDofType>>().ToArray();
				foreach (var boundaryCondition in transientBCs)
				{
					boundaryCondition.CurrentTime = time;
				}

				return transientBCs
					.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id)))
					.OfType<INodalHeatFluxBoundaryCondition>()
					.Where(x => model.EnumerateBoundaryConditions(id)
						.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id)))
						.OfType<INodalTemperatureBoundaryCondition>()
						.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false);
			},
				rhs);

			return rhs;
		}

		public double CalculateRhsNorm(IVector rhs) => rhs.Norm2();

		public void ProcessInternalRhs(IVector solution, IVector rhs) 
		{
			// Method intentionally left blank
		}

		public IVector CalculateResponseIntegralVector(IVector solution)
		{
			throw new NotImplementedException();
		}

		public void UpdateState(IHaveState externalState)
		{
			throw new NotImplementedException();
		}

		public IEnumerable<INodalNeumannBoundaryCondition<IDofType>> EnumerateEquivalentNeumannBoundaryConditions(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateEquivalentNodalNeumannBoundaryConditions(model.EnumerateElements(subdomainID), Array.Empty<(int NodeID, IDofType DOF)>()))
				.OfType<INodalHeatFluxBoundaryCondition>()
				.Where(x => model.EnumerateBoundaryConditions(subdomainID)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(subdomainID)))
					.OfType<INodalTemperatureBoundaryCondition>()
					.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false);

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering() =>
			model.EnumerateSubdomains()
				.Select(x => new Tuple<int, IEnumerable<IBoundaryConditionSet<IDofType>>>(x.ID, model.EnumerateBoundaryConditions(x.ID)))
					.Select(i => i.Item2
						.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(i.Item1)))
						.OfType<INodalTemperatureBoundaryCondition>())
					.SelectMany(x => x)
					.OrderBy(x => x.Node.ID)
					.GroupBy(x => (x.Node.ID, x.DOF))
					.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount)))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(subdomainID))).OfType<INodalTemperatureBoundaryCondition>()
				.OrderBy(x => x.Node.ID)
				.GroupBy(x => (x.Node.ID, x.DOF))
				.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount)))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));

		public IMatrix GetMatrix() => GetMatrix(DifferentiationOrder.Zero);

		public IVector GetRhs() => GetRhs(0);
	}
}
