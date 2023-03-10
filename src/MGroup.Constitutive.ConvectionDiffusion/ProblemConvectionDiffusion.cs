using System.Linq;
using System.Collections.Generic;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;
using MGroup.Constitutive.ConvectionDiffusion.Providers;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.DataStructures;
using System;

namespace MGroup.Constitutive.ConvectionDiffusion
{
	public class ProblemConvectionDiffusion : IAlgebraicModelInterpreter, ITransientAnalysisProvider, INonTransientAnalysisProvider, INonLinearProvider
	{
		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ElementProductionProvider productionProvider = new ElementProductionProvider();
		private readonly ElementIndependentProductionProvider independentProductionProvider = new ElementIndependentProductionProvider();
		private readonly ElementDiffusionProvider diffusionProvider = new ElementDiffusionProvider();
		private readonly ElementConvectionProvider convectionProvider = new ElementConvectionProvider();
		private readonly ElementCapacityMatrixProvider fistTimeDerivativeMatrixProvider = new ElementCapacityMatrixProvider();
		private readonly IElementMatrixPredicate rebuildDiffusionPredicate = new MaterialModifiedElementMarixPredicate();
		private IGlobalMatrix convection, diffusion, production, capacityMatrix;
		private TransientAnalysisPhase analysisPhase = TransientAnalysisPhase.SteadyStateSolution;

		public ProblemConvectionDiffusion(IModel model, IAlgebraicModel algebraicModel)
		{
			this.model = model;
			this.algebraicModel = algebraicModel;
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

		private IGlobalMatrix CapacityMatrix
		{
			get
			{
				if (capacityMatrix == null) BuildCapacityMatrix();
				return capacityMatrix;
			}
		}

		public ActiveDofs ActiveDofs { get; } = new ActiveDofs();

		public DifferentiationOrder ProblemOrder => DifferentiationOrder.First;

		public void SetTransientAnalysisPhase(TransientAnalysisPhase phase) => analysisPhase = phase;

		private void BuildConvection() => convection = algebraicModel.BuildGlobalMatrix(convectionProvider);

		private void BuildDiffusion() => diffusion = algebraicModel.BuildGlobalMatrix(diffusionProvider);

		private void BuildProduction() => production = algebraicModel.BuildGlobalMatrix(productionProvider);

		private void BuildCapacityMatrix() => capacityMatrix = algebraicModel.BuildGlobalMatrix(fistTimeDerivativeMatrixProvider);

		private void RebuildDiffusion() => algebraicModel.RebuildGlobalMatrixPartially(diffusion, 
			model.EnumerateElements, diffusionProvider, rebuildDiffusionPredicate);

		public void Reset()
		{
			convection = null;
			diffusion = null;
			production = null;
			capacityMatrix = null;
		}

		public IGlobalMatrix GetMatrix()
		{
			IGlobalMatrix matrix = Diffusion.Copy();
			matrix.AddIntoThis(Convection);
			matrix.AddIntoThis(Production);

			return matrix;
		}

		public IGlobalMatrix GetMatrix(DifferentiationOrder differentiationOrder) => differentiationOrder switch
		{
			DifferentiationOrder.Zero => GetMatrix(),
			DifferentiationOrder.First => CapacityMatrix,
			_ => algebraicModel.CreateEmptyMatrix(),
		};

		private IGlobalVector GetFirstOrderDerivativeVectorFromBoundaryConditions(double time)
		{
			var boundaryConditions = model.EnumerateBoundaryConditions(model.EnumerateSubdomains().First().ID).ToArray();
			foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IConvectionDiffusionDofType>>())
			{
				boundaryCondition.CurrentTime = time;
			}

			IGlobalVector capacities = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(
			id =>
			{
				var boundaryConditions = model.EnumerateBoundaryConditions(id).ToArray();
				foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IConvectionDiffusionDofType>>())
				{
					boundaryCondition.CurrentTime = time;
				}

				return boundaryConditions
					.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id)))
					.OfType<INodalCapacityBoundaryCondition>();
			},
			capacities);

			//algebraicModel.AddToGlobalVector(boundaryConditions
			//	.SelectMany(x => x.EnumerateDomainBoundaryConditions())
			//	.OfType<IDomainCapacityBoundaryCondition>(),
			//capacities);

			return capacities;
		}

		private IGlobalVector GetZeroOrderDerivativeVectorFromInitialConditions(double time)
		{
			IGlobalVector unknownVariables = algebraicModel.CreateZeroVector();

			if (time == 0)
			{
				algebraicModel.AddToGlobalVector(
				id =>
				{
					var initialConditions = model.EnumerateInitialConditions(id).ToArray();
					return initialConditions
						.SelectMany(x => x.EnumerateNodalInitialConditions(model.EnumerateElements(id)))
						.OfType<INodalUnknownVariableInitialCondition>();
				},
				unknownVariables);

				//algebraicModel.AddToGlobalVector(model.EnumerateInitialConditions(model.EnumerateSubdomains().First().ID)
				//	.SelectMany(x => x.EnumerateDomainInitialConditions())
				//	.OfType<IDomainUnknownVariableInitialCondition>(),
				//unknownVariables);
			}

			return unknownVariables;
		}

		public IGlobalVector GetVectorFromModelConditions(DifferentiationOrder differentiationOrder, double time) => differentiationOrder switch
		{
			DifferentiationOrder.Zero => GetZeroOrderDerivativeVectorFromInitialConditions(time),
			DifferentiationOrder.First => GetFirstOrderDerivativeVectorFromBoundaryConditions(time),
			_ => algebraicModel.CreateZeroVector(),
		};

		public IGlobalVector GetRhs(double time)
		{
			var rhs = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(rhs, independentProductionProvider);

			algebraicModel.AddToGlobalVector(id =>
				model.EnumerateBoundaryConditions(id)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id)))
					.OfType<INodalUnknownVariableFluxBoundaryCondition>()
					.Where(x => model.EnumerateBoundaryConditions(id)
						.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id)))
						.OfType<INodalUnknownVariableBoundaryCondition>()
						.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false),
					rhs);

			algebraicModel.AddToGlobalVector(id =>
			{
				var transientBCs = model.EnumerateBoundaryConditions(id).OfType<ITransientBoundaryConditionSet<IConvectionDiffusionDofType>>().ToArray();
				foreach (var boundaryCondition in transientBCs)
				{
					boundaryCondition.CurrentTime = time;
				}

				return transientBCs
					.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id)))
					.OfType<INodalUnknownVariableFluxBoundaryCondition>()
					.Where(x => model.EnumerateBoundaryConditions(id)
						.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id)))
						.OfType<INodalUnknownVariableBoundaryCondition>()
						.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false);
			},
				rhs);
			
			return rhs;
		}

		public double CalculateRhsNorm(IGlobalVector rhs) => rhs.Norm2();

		public void ProcessInternalRhs(IGlobalVector solution, IGlobalVector rhs) 
		{
			// Method intentionally left blank
		}

		public IGlobalVector CalculateResponseIntegralVector(IGlobalVector solution)
		{
			throw new NotImplementedException();
		}

		public void UpdateState(IHaveState externalState)
		{
			throw new NotImplementedException();
		}

		public IEnumerable<INodalNeumannBoundaryCondition<IDofType>> EnumerateEquivalentNeumannBoundaryConditions(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateEquivalentNodalNeumannBoundaryConditions(model.EnumerateElements(subdomainID)))
				.OfType<INodalUnknownVariableFluxBoundaryCondition>()
				.Where(x => model.EnumerateBoundaryConditions(subdomainID)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(subdomainID)))
					.OfType<INodalUnknownVariableBoundaryCondition>()
					.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false);

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering() =>
			model.EnumerateSubdomains()
				.Select(x => new Tuple<int, IEnumerable<IBoundaryConditionSet<IDofType>>>(x.ID, model.EnumerateBoundaryConditions(x.ID)))
					.Select(i => i.Item2
						.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(i.Item1)))
						.OfType<INodalUnknownVariableBoundaryCondition>())
					.SelectMany(x => x)
					.OrderBy(x => x.Node.ID)
					.GroupBy(x => (x.Node.ID, x.DOF))
					.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount)))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(subdomainID))).OfType<INodalUnknownVariableBoundaryCondition>()
				.OrderBy(x => x.Node.ID)
				.GroupBy(x => (x.Node.ID, x.DOF))
				.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount)))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));

		public IGlobalVector GetRhs() => GetRhs(0);
	}
}
