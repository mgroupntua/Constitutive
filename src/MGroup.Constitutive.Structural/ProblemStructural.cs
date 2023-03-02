using System.Collections.Generic;
using System.Linq;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.Constitutive.Structural.Providers;
using MGroup.Constitutive.Structural.InitialConditions;
using MGroup.MSolve.AnalysisWorkflow.Providers;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.Structural
{
	public class ProblemStructural : IAlgebraicModelInterpreter, ITransientAnalysisProvider, INonTransientAnalysisProvider, INonLinearProvider
	{
		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ElementStructuralStiffnessProvider stiffnessProvider = new ElementStructuralStiffnessProvider();
		private readonly ElementStructuralMassProvider massProvider = new ElementStructuralMassProvider();
		private readonly ElementStructuralDampingProvider dampingProvider = new ElementStructuralDampingProvider();
		private readonly ElementStructuralInternalForcesProvider rhsProvider = new ElementStructuralInternalForcesProvider();
		private readonly ElementStructuralVolumeLoadsProvider integratedLoadsVolumeProvider = new ElementStructuralVolumeLoadsProvider();
		private readonly IElementMatrixPredicate rebuildStiffnessPredicate = new MaterialModifiedElementMarixPredicate();
		private IGlobalMatrix mass, damping, stiffness;
		private TransientAnalysisPhase analysisPhase = TransientAnalysisPhase.SteadyStateSolution;

		public ProblemStructural(IModel model, IAlgebraicModel algebraicModel)
		{
			this.model = model;
			this.algebraicModel = algebraicModel;
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

		public DifferentiationOrder ProblemOrder => DifferentiationOrder.Second;

		public void SetTransientAnalysisPhase(TransientAnalysisPhase phase) => analysisPhase = phase;

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

		public IGlobalMatrix GetMatrix(DifferentiationOrder differentiationOrder) => differentiationOrder switch
		{
			DifferentiationOrder.Zero => Stiffness,
			DifferentiationOrder.First => Damping,
			DifferentiationOrder.Second => Mass,
			_ => algebraicModel.CreateEmptyMatrix(),
		};

		private IGlobalVector GetSecondOrderDerivativeVectorFromBoundaryConditions(double time)
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

		private IGlobalVector GetFirstOrderDerivativeVectorFromBoundaryConditions(double time)
		{
			var boundaryConditions = model.EnumerateBoundaryConditions(model.EnumerateSubdomains().First().ID).ToArray();
			foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IStructuralDofType>>())
			{
				boundaryCondition.CurrentTime = time;
			}

			IGlobalVector velocities = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(id =>
			{
				var boundaryConditions = model.EnumerateBoundaryConditions(id).ToArray();
				foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IStructuralDofType>>())
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

			if (time == 0)
			{
				algebraicModel.AddToGlobalVector(id =>
				{
					var initialConditions = model.EnumerateInitialConditions(id).ToArray();
					return initialConditions
						.SelectMany(x => x.EnumerateNodalInitialConditions())
						.OfType<INodalVelocityInitialCondition>();
				},
					velocities);
			}

			return velocities;
		}

		private IGlobalVector GetZeroOrderDerivativeVectorFromInitialConditions(double time)
		{
			IGlobalVector displacements = algebraicModel.CreateZeroVector();
			if (time == 0)
			{
				algebraicModel.AddToGlobalVector(id =>
				{
					var initialConditions = model.EnumerateInitialConditions(id).ToArray();
					return initialConditions
						.SelectMany(x => x.EnumerateNodalInitialConditions())
						.OfType<INodalDisplacementInitialCondition>();
				},
					displacements);

				algebraicModel.AddToGlobalVector(model.EnumerateInitialConditions(model.EnumerateSubdomains().First().ID)
					.SelectMany(x => x.EnumerateDomainInitialConditions())
					.OfType<IDomainDisplacementInitialCondition>(),
				displacements);
			}

			return displacements;
		}

		public IGlobalVector GetVectorFromModelConditions(DifferentiationOrder differentiationOrder, double time) => differentiationOrder switch
		{
			DifferentiationOrder.Zero => GetZeroOrderDerivativeVectorFromInitialConditions(time),
			DifferentiationOrder.First => GetFirstOrderDerivativeVectorFromBoundaryConditions(time),
			DifferentiationOrder.Second => GetSecondOrderDerivativeVectorFromBoundaryConditions(time),
			_ => algebraicModel.CreateZeroVector(),
		};

		public IGlobalVector GetRhs(double time)
		{
			IGlobalVector rhs = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(rhs, integratedLoadsVolumeProvider);


			algebraicModel.AddToGlobalVector(id =>
				model.EnumerateBoundaryConditions(id)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalLoadBoundaryCondition>()
					.Where(x => model.EnumerateBoundaryConditions(id)
						.SelectMany(x => x.EnumerateNodalBoundaryConditions())
						.OfType<INodalDisplacementBoundaryCondition>()
						.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false),
				rhs);

			algebraicModel.AddToGlobalVector(id =>
			{
				var transientBCs = model.EnumerateBoundaryConditions(id).OfType<ITransientBoundaryConditionSet<IStructuralDofType>>().ToArray();
				foreach (var boundaryCondition in transientBCs)
				{
					boundaryCondition.CurrentTime = time;
				}

				return transientBCs
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalLoadBoundaryCondition>()
					.Where(x => model.EnumerateBoundaryConditions(id)
						.SelectMany(x => x.EnumerateNodalBoundaryConditions())
						.OfType<INodalDisplacementBoundaryCondition>()
						.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false);
			},
				rhs);
			//algebraicModel.AddToGlobalVector(EnumerateEquivalentNeumannBoundaryConditions, rhs);

			return rhs;
		}

		//TODO: I suggest splitting this into 2 methods. One for updating the elements/materials and one for calculating the internal rhs
		public IGlobalVector CalculateResponseIntegralVector(IGlobalVector solution)
		{
			IGlobalVector internalRhs = algebraicModel.CreateZeroVector();
			if (analysisPhase != TransientAnalysisPhase.InitialConditionEvaluation)
			{
				var dirichletBoundaryConditions = algebraicModel.BoundaryConditionsInterpreter.GetDirichletBoundaryConditionsWithNumbering()
				.Select(x => new NodalBoundaryCondition(x.Value.Node, x.Key.DOF, x.Value.Amount));
				// First update the state of the elements
				algebraicModel.DoPerElement<IStructuralElementType>(element =>
				{
					double[] elementDisplacements = algebraicModel.ExtractElementVector(solution, element);
					element.MapNodalBoundaryConditionsToElementVector(dirichletBoundaryConditions, elementDisplacements);
					element.CalculateResponse(elementDisplacements);
				});

				// Then calculate the internal rhs vector
				algebraicModel.AddToGlobalVector(internalRhs, rhsProvider);
			}
			else
			{
				Mass.MultiplyVector(solution, internalRhs);
			}

			return internalRhs;
		}

		public void UpdateState(IHaveState externalState)
		{
			if (analysisPhase != TransientAnalysisPhase.InitialConditionEvaluation)
			{
				algebraicModel.DoPerElement<IStructuralElementType>(element =>
				{
					element.SaveConstitutiveLawState(externalState);
				});
			}
		}

		public IGlobalVector GetRHSFromSolutionWithInitialDisplacemntsEffect(IGlobalVector solution, Dictionary<int, INode> boundaryNodes,
			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements)
		{
			// First update the state of the elements
			algebraicModel.DoPerElement<IElementType>(element =>
			{
				double[] localSolution = algebraicModel.ExtractElementVector(solution, element);
				ImposePrescribedDisplacementsWithInitialConditionSEffect(element, localSolution, boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
				element.CalculateResponse(localSolution);
			});

			// Then calculate the internal rhs vector
			IGlobalVector internalRhs = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(internalRhs, rhsProvider);
			return internalRhs;
		}

		public void ImposePrescribedDisplacementsWithInitialConditionSEffect(IElementType element, double[] localSolution, Dictionary<int, INode> boundaryNodes,
			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements)
		{
			var elementDOFTypes = element.DofEnumerator.GetDofTypesForMatrixAssembly(element);
			var matrixAssemblyNodes = element.DofEnumerator.GetNodesForMatrixAssembly(element);
			int iElementMatrixColumn = 0;
			for (int j = 0; j < elementDOFTypes.Count; j++)
			{
				INode nodeColumn = matrixAssemblyNodes[j];
				int nodalDofsNumber = elementDOFTypes[j].Count;
				if (boundaryNodes.ContainsKey(nodeColumn.ID))
				{
					Dictionary<IDofType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
					Dictionary<IDofType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
					int positionOfDofInNode = 0;
					foreach (IDofType doftype1 in elementDOFTypes[j])
					{
						if (nodalConvergedDisplacements.ContainsKey(doftype1))
						{
							localSolution[iElementMatrixColumn + positionOfDofInNode] = nodalConvergedDisplacements[doftype1] + (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
							// TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
						}
						positionOfDofInNode += 1;
					}
				}
				iElementMatrixColumn += nodalDofsNumber;
			}

		}

		public void ImposePrescribed_d_DisplacementsWithInitialConditionSEffect(IElementType element, double[] localSolution, Dictionary<int, INode> boundaryNodes,
			Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
			int nIncrement, int totalIncrements)
		{

			var elementDOFTypes = element.DofEnumerator.GetDofTypesForMatrixAssembly(element);
			var matrixAssemblyNodes = element.DofEnumerator.GetNodesForMatrixAssembly(element);
			int iElementMatrixColumn = 0;
			for (int j = 0; j < elementDOFTypes.Count; j++)
			{
				INode nodeColumn = matrixAssemblyNodes[j];
				int nodalDofsNumber = elementDOFTypes[j].Count;
				if (boundaryNodes.ContainsKey(nodeColumn.ID))
				{
					Dictionary<IDofType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
					Dictionary<IDofType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
					int positionOfDofInNode = 0;
					foreach (IDofType doftype1 in elementDOFTypes[j])
					{
						if (nodalConvergedDisplacements.ContainsKey(doftype1))
						{
							localSolution[iElementMatrixColumn + positionOfDofInNode] = (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
							// 1) den vazoume mono (1/increments) alla (nIncrement/increments) dioti metaxu aftwn twn nIncrements den exei mesolavhsei save sta material ths mikroklimakas
							// TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
						}
						positionOfDofInNode += 1;
					}
				}
				iElementMatrixColumn += nodalDofsNumber;
			}
		}

		public IEnumerable<INodalNeumannBoundaryCondition<IDofType>> EnumerateEquivalentNeumannBoundaryConditions(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateEquivalentNodalNeumannBoundaryConditions(model.EnumerateElements(subdomainID)))
				.OfType<INodalLoadBoundaryCondition>()
				.Where(x => model.EnumerateBoundaryConditions(subdomainID)
					.SelectMany(x => x.EnumerateNodalBoundaryConditions())
					.OfType<INodalDisplacementBoundaryCondition>()
					.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false);

		public double CalculateRhsNorm(IGlobalVector rhs) => rhs.Norm2();

		public void ProcessInternalRhs(IGlobalVector solution, IGlobalVector rhs) 
		{
			// Method intentionally left blank
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

		public IGlobalVector GetRhs() => GetRhs(0);

		public IGlobalMatrix GetMatrix() => GetMatrix(DifferentiationOrder.Zero);
	}
}
