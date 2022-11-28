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
using MGroup.MSolve.AnalysisWorkflow.Transient;
using System;
using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;

namespace MGroup.Constitutive.ConvectionDiffusion
{
	public class ProblemConvectionDiffusion : IAlgebraicModelInterpreter, ITransientAnalysisProvider, INonTransientAnalysisProvider, INonLinearProvider
	{
		private bool shouldRebuildDiffusionMatrixForZeroOrderDerivativeMatrixVectorProduct = false;
		private IGlobalMatrix convection, diffusion, production, capacityMatrix;
		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ISolver solver;
		private ElementProductionProvider productionProvider = new ElementProductionProvider();
		private ElementIndependentProductionProvider independentProductionProvider = new ElementIndependentProductionProvider();
		private ElementDiffusionProvider diffusionProvider = new ElementDiffusionProvider();
		private ElementConvectionProvider convectionProvider = new ElementConvectionProvider();
		private ElementCapacityMatrixProvider fistTimeDerivativeMatrixProvider = new ElementCapacityMatrixProvider();

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

		private IGlobalMatrix CapacityMatrix
		{
			get
			{
				if (capacityMatrix == null) BuildCapacityMatrix();
				return capacityMatrix;
			}
		}

		public ActiveDofs ActiveDofs { get; } = new ActiveDofs();

		private void BuildConvection() => convection = algebraicModel.BuildGlobalMatrix(convectionProvider);

		private void BuildDiffusion() => diffusion = algebraicModel.BuildGlobalMatrix(diffusionProvider);

		private void BuildProduction() => production = algebraicModel.BuildGlobalMatrix(productionProvider);

		private void BuildCapacityMatrix() => capacityMatrix = algebraicModel.BuildGlobalMatrix(fistTimeDerivativeMatrixProvider);

		private void RebuildDiffusion() => algebraicModel.RebuildGlobalMatrixPartially(diffusion, 
			model.EnumerateElements, diffusionProvider, rebuildDiffusionPredicate);

		#region IAnalyzerProvider Members
		public void ClearMatrices()
		{
			convection = null;
			diffusion = null;
			production = null;
			capacityMatrix = null;
		}

		public void Reset()
		{
			ClearMatrices();
		}

		public void GetProblemDofTypes()
		{
			// function no longer used
		}
		#endregion

		#region IImplicitIntegrationProvider Members

		public void LinearCombinationOfMatricesIntoEffectiveMatrix(TransientAnalysisCoefficients coefficients)
		{
			IGlobalMatrix matrix = Diffusion;
			matrix.AddIntoThis(Convection);
			matrix.AddIntoThis(Production);
			if (coefficients[DifferentiationOrder.First] != 0)
			{
				matrix.LinearCombinationIntoThis(coefficients[DifferentiationOrder.Zero], CapacityMatrix, coefficients[DifferentiationOrder.First]);
			}
			else
			{
				matrix.ScaleIntoThis(coefficients[DifferentiationOrder.Zero]);
			}

			solver.LinearSystem.Matrix = matrix;
			shouldRebuildDiffusionMatrixForZeroOrderDerivativeMatrixVectorProduct = true;
		}

		public void LinearCombinationOfMatricesIntoEffectiveMatrixNoOverwrite(TransientAnalysisCoefficients coefficients)
		{
			IGlobalMatrix matrix = Diffusion.Copy();
			matrix.AddIntoThis(Convection);
			matrix.AddIntoThis(Production);
			matrix.AxpyIntoThis(CapacityMatrix, coefficients[DifferentiationOrder.First]);
			solver.LinearSystem.Matrix = matrix;
		}

		public void ProcessRhs(TransientAnalysisCoefficients coefficients, IGlobalVector rhs)
		{
			// Method intentionally left empty.
		}

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
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalCapacityBoundaryCondition>();
			},
			capacities);

			algebraicModel.AddToGlobalVector(boundaryConditions
				.SelectMany(x => x.EnumerateDomainBoundaryConditions())
				.OfType<IDomainCapacityBoundaryCondition>(),
			capacities);

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
						.SelectMany(x => x.EnumerateNodalInitialConditions())
						.OfType<INodalUnknownVariableInitialCondition>();
				},
				unknownVariables);

				algebraicModel.AddToGlobalVector(model.EnumerateInitialConditions(model.EnumerateSubdomains().First().ID)
					.SelectMany(x => x.EnumerateDomainInitialConditions())
					.OfType<IDomainUnknownVariableInitialCondition>(),
				unknownVariables);
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
			solver.LinearSystem.RhsVector.Clear();
			AssignRhs();

			algebraicModel.AddToGlobalVector(id =>
			{
				var transientBCs = model.EnumerateBoundaryConditions(id).OfType<ITransientBoundaryConditionSet<IConvectionDiffusionDofType>>().ToArray();
				foreach (var boundaryCondition in transientBCs)
				{
					boundaryCondition.CurrentTime = time;
				}

				return transientBCs
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalUnknownVariableFluxBoundaryCondition>()
					.Where(x => model.EnumerateBoundaryConditions(id)
						.SelectMany(x => x.EnumerateNodalBoundaryConditions())
						.OfType<INodalUnknownVariableBoundaryCondition>()
						.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false);
			},
				solver.LinearSystem.RhsVector);
			algebraicModel.AddToGlobalVector(EnumerateEquivalentNeumannBoundaryConditions, solver.LinearSystem.RhsVector); 
			
			IGlobalVector result = solver.LinearSystem.RhsVector.Copy();
			return result;
		}

		private IGlobalVector SecondOrderDerivativeMatrixVectorProduct(IGlobalVector vector) => algebraicModel.CreateZeroVector();

		private IGlobalVector FirstOrderDerivativeMatrixVectorProduct(IGlobalVector vector)
		{
			IGlobalVector result = algebraicModel.CreateZeroVector();
			CapacityMatrix.MultiplyVector(vector, result);
			return result;
		}

		private IGlobalVector ZeroOrderDerivativeMatrixVectorProduct(IGlobalVector vector)
		{
			if (shouldRebuildDiffusionMatrixForZeroOrderDerivativeMatrixVectorProduct)
			{
				BuildDiffusion();
				shouldRebuildDiffusionMatrixForZeroOrderDerivativeMatrixVectorProduct = false;
			}

			IGlobalVector result = algebraicModel.CreateZeroVector();
			diffusion.MultiplyVector(vector, result);
			
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
			if (diffusion == null) BuildDiffusion();
			if (convection == null) BuildConvection();
			if (production == null) BuildProduction();

			var effectiveMatrix = diffusion.Copy();
			effectiveMatrix.AddIntoThis(convection);
			effectiveMatrix.AddIntoThis(production);
			solver.LinearSystem.Matrix = effectiveMatrix;
		}
		#endregion

		#region INonLinearProvider Members

		public double CalculateRhsNorm(IGlobalVector rhs) => rhs.Norm2();

		public void ProcessInternalRhs(IGlobalVector solution, IGlobalVector rhs) { }

		#endregion

		public void AssignRhs()
		{
			solver.LinearSystem.RhsVector.Clear();

			algebraicModel.AddToGlobalVector(solver.LinearSystem.RhsVector, independentProductionProvider);

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
