#pragma warning disable SA1116 // Split parameters should start on line after declaration

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
using System.Net.Sockets;
using System;
using System.Diagnostics;

namespace MGroup.Constitutive.Structural
{
	public enum BoundaryConditionPrecedence
	{
		/// <summary>
		/// If a node has both a load and a displacement boundary condition, an exception is raised
		/// </summary>
		NoPrecedence = 0,
		/// <summary>
		/// If a node has both a load and a displacement boundary condition, the displacement boundary condition is IGNORED
		/// </summary>
		LoadPrecedence,
		/// <summary>
		/// If a node has both a load and a displacement boundary condition, the load boundary condition is IGNORED
		/// </summary>
		DisplacementPrecedence
	}

	public class ProblemStructural : IAlgebraicModelInterpreter, ITransientAnalysisProvider, INonTransientAnalysisProvider, INonLinearProvider
	{
		private readonly IModel model;
		private readonly IAlgebraicModel algebraicModel;
		private readonly ElementStructuralStiffnessProvider stiffnessProvider = new ElementStructuralStiffnessProvider();
		private readonly ElementStructuralMassProvider massProvider = new ElementStructuralMassProvider();
		private readonly ElementStructuralDampingProvider dampingProvider = new ElementStructuralDampingProvider();
		private readonly ElementStructuralInternalForcesProvider rhsProvider = new ElementStructuralInternalForcesProvider();
		private readonly IElementMatrixPredicate rebuildStiffnessPredicate = new MaterialModifiedElementMarixPredicate();
		private readonly BoundaryConditionPrecedence precedence;
		private IGlobalMatrix mass, damping, stiffness;
		private TransientAnalysisPhase analysisPhase = TransientAnalysisPhase.SteadyStateSolution;

		public ProblemStructural(IModel model, IAlgebraicModel algebraicModel, BoundaryConditionPrecedence precedence)
		{
			this.precedence = precedence;
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

		public ProblemStructural(IModel model, IAlgebraicModel algebraicModel)
			: this(model, algebraicModel, BoundaryConditionPrecedence.NoPrecedence)
		{
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

		private static IEnumerable<T> GetNeumannBoundaryConditionsAccordingToPrecedence<T>(INodalBoundaryCondition<IDofType>[] allBCs, BoundaryConditionPrecedence precedence)
			where T : INodalBoundaryCondition<IStructuralDofType>
		{
			var nodalNeumann = allBCs.OfType<T>().ToArray();
			var nodalDirichlet = allBCs.OfType<INodalDisplacementBoundaryCondition>();
			var validNeumann = nodalNeumann.Where(x => nodalDirichlet.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false).ToArray();
			if (precedence == BoundaryConditionPrecedence.NoPrecedence && nodalNeumann.Length != validNeumann.Length)
			{
				var duplicateLoadNodes = nodalNeumann.Except(validNeumann).Select(x => $"({x.Node.ID}, {x.DOF})").Aggregate(String.Empty, (s, n) => s + n + ", ");
				throw new ArgumentException($"Boundary conditions for both loads/velocities/accelerations and displacements at nodes {duplicateLoadNodes}");
			}

			return precedence == BoundaryConditionPrecedence.LoadPrecedence ? nodalNeumann : validNeumann;
		}

		private IEnumerable<T> GetNeumannBoundaryConditionsAccordingToPrecedence<T>(INodalBoundaryCondition<IDofType>[] allBCs)
			where T : INodalBoundaryCondition<IStructuralDofType> => GetNeumannBoundaryConditionsAccordingToPrecedence<T>(allBCs, precedence);

		private IEnumerable<T> GetNeumannInitialConditionsAccordingToPrecedence<T>(INodalInitialCondition<IDofType>[] allBCs)
			where T : INodalInitialCondition<IStructuralDofType>
		{
			var nodalNeumann = allBCs.OfType<T>().ToArray();
			var nodalDirichlet = allBCs.OfType<INodalDisplacementBoundaryCondition>();
			var validNeumann = nodalNeumann.Where(x => nodalDirichlet.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false).ToArray();
			if (precedence == BoundaryConditionPrecedence.NoPrecedence && nodalNeumann.Length != validNeumann.Length)
			{
				var duplicateLoadNodes = nodalNeumann.Except(validNeumann).Select(x => $"({x.Node.ID}, {x.DOF})").Aggregate(String.Empty, (s, n) => s + n + ", ");
				throw new ArgumentException($"Boundary conditions for both loads/velocities/accelerations and displacements at nodes {duplicateLoadNodes}");
			}

			return precedence == BoundaryConditionPrecedence.LoadPrecedence ? nodalNeumann : validNeumann;
		}

		private IEnumerable<T> GetDirichletBoundaryConditionsAccordingToPrecedence<T>(INodalBoundaryCondition<IDofType>[] allBCs)
			where T : INodalBoundaryCondition<IStructuralDofType>
		{
			var nodalDirichlet = allBCs.OfType<T>().ToArray();
			var nodalNeumann = allBCs.OfType<INodalLoadBoundaryCondition>();
			var validDirichlet = nodalDirichlet.Where(x => nodalNeumann.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false).ToArray();
			if (precedence == BoundaryConditionPrecedence.NoPrecedence && nodalDirichlet.Length != validDirichlet.Length)
			{
				var duplicateLoadNodes = nodalDirichlet.Except(validDirichlet).Select(x => $"({x.Node.ID}, {x.DOF})").Aggregate(String.Empty, (s, n) => s + n + ", ");
				throw new ArgumentException($"Boundary conditions for both loads/velocities/accelerations and displacements at nodes {duplicateLoadNodes}");
			}

			return precedence == BoundaryConditionPrecedence.DisplacementPrecedence ? nodalDirichlet : validDirichlet;
		}

		private IEnumerable<T> GetDirichletInitialConditionsAccordingToPrecedence<T>(INodalInitialCondition<IDofType>[] allBCs)
			where T : INodalInitialCondition<IStructuralDofType>
		{
			var nodalDirichlet = allBCs.OfType<T>().ToArray();
			var nodalNeumann = allBCs.OfType<INodalLoadBoundaryCondition>();
			var validDirichlet = nodalDirichlet.Where(x => nodalNeumann.Any(d => d.Node.ID == x.Node.ID && d.DOF == x.DOF) == false).ToArray();
			if (precedence == BoundaryConditionPrecedence.NoPrecedence && nodalDirichlet.Length != validDirichlet.Length)
			{
				var duplicateLoadNodes = nodalDirichlet.Except(validDirichlet).Select(x => $"({x.Node.ID}, {x.DOF})").Aggregate(String.Empty, (s, n) => s + n + ", ");
				throw new ArgumentException($"Boundary conditions for both loads/velocities/accelerations and displacements at nodes {duplicateLoadNodes}");
			}

			return precedence == BoundaryConditionPrecedence.DisplacementPrecedence ? nodalDirichlet : validDirichlet;
		}

		private IGlobalVector GetSecondOrderDerivativeVectorFromBoundaryConditions(double time)
		{
			IGlobalVector accelerations = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(id =>
			{
				var boundaryConditions = model.EnumerateBoundaryConditions(id).ToArray();
				foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IStructuralDofType>>())
				{
					boundaryCondition.CurrentTime = time;
				}

				return GetNeumannBoundaryConditionsAccordingToPrecedence<INodalAccelerationBoundaryCondition>(boundaryConditions.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id))).ToArray());
			},
				accelerations);

			return accelerations;
		}

		private IGlobalVector GetFirstOrderDerivativeVectorFromBoundaryConditions(double time)
		{
			IGlobalVector velocities = algebraicModel.CreateZeroVector();
			algebraicModel.AddToGlobalVector(id =>
			{
				var boundaryConditions = model.EnumerateBoundaryConditions(id).ToArray();
				foreach (var boundaryCondition in boundaryConditions.OfType<ITransientBoundaryConditionSet<IStructuralDofType>>())
				{
					boundaryCondition.CurrentTime = time;
				}

				return GetNeumannBoundaryConditionsAccordingToPrecedence<INodalVelocityBoundaryCondition>(boundaryConditions.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id))).ToArray());
			},
				velocities);

			if (time == 0)
			{
				algebraicModel.AddToGlobalVector(id =>
				{
					var initialConditions = model.EnumerateInitialConditions(id).ToArray();
					return GetNeumannInitialConditionsAccordingToPrecedence<INodalVelocityInitialCondition>(initialConditions.SelectMany(x => x.EnumerateNodalInitialConditions(model.EnumerateElements(id))).ToArray());
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
					return GetDirichletInitialConditionsAccordingToPrecedence<INodalDisplacementInitialCondition>(initialConditions.SelectMany(x => x.EnumerateNodalInitialConditions(model.EnumerateElements(id))).ToArray());
				},
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

			algebraicModel.AddToGlobalVector(id => GetNeumannBoundaryConditionsAccordingToPrecedence<INodalLoadBoundaryCondition>(model.EnumerateBoundaryConditions(id)
				.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id)))
				.ToArray()),
				rhs);
			algebraicModel.AddToGlobalVector(id =>
			{
				var transientBCs = model.EnumerateBoundaryConditions(id).OfType<ITransientBoundaryConditionSet<IStructuralDofType>>().ToArray();
				foreach (var boundaryCondition in transientBCs)
				{
					boundaryCondition.CurrentTime = time;
				}

				return GetNeumannBoundaryConditionsAccordingToPrecedence<INodalLoadBoundaryCondition>(transientBCs.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(id))).ToArray());
			},
				rhs);

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
				algebraicModel.DoPerElement<IStructuralElementType>(element =>
				{
					double[] elementDisplacements = algebraicModel.ExtractElementVector(solution, element);
					element.MapNodalBoundaryConditionsToElementVector(dirichletBoundaryConditions, elementDisplacements);
					element.CalculateResponse(elementDisplacements);
				});

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

		public IEnumerable<INodalNeumannBoundaryCondition<IDofType>> EnumerateEquivalentNeumannBoundaryConditions(int subdomainID)
		{
			var dirichletBoundaryConditions = model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(subdomainID)))
				.OfType<INodalDisplacementBoundaryCondition>()
				.ToArray();

			var dofsToExclude = precedence != BoundaryConditionPrecedence.LoadPrecedence
				? Array.Empty<(int NodeID, IDofType DOF)>()
				: dirichletBoundaryConditions
					.Select(x => (x.Node.ID, (IDofType)x.DOF))
					.Distinct()
					.ToArray();

			var nodalNeumannBoundaryConditions = model.EnumerateBoundaryConditions(subdomainID)
					.SelectMany(x => x.EnumerateEquivalentNodalNeumannBoundaryConditions(model.EnumerateElements(subdomainID), dofsToExclude))
					.Select(x => (INodalBoundaryCondition<IStructuralDofType>)x)
					.Concat(dirichletBoundaryConditions)
					.ToArray();

			return GetNeumannBoundaryConditionsAccordingToPrecedence<INodalLoadBoundaryCondition>(nodalNeumannBoundaryConditions, BoundaryConditionPrecedence.DisplacementPrecedence);
		}

		public double CalculateRhsNorm(IGlobalVector rhs) => rhs.Norm2();

		public void ProcessInternalRhs(IGlobalVector solution, IGlobalVector rhs) 
		{
			// Method intentionally left blank
		}

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering() =>
			model.EnumerateSubdomains()
				.Select(x => new Tuple<int, IEnumerable<IBoundaryConditionSet<IDofType>>>(x.ID, model.EnumerateBoundaryConditions(x.ID)))
					.Select(i => i.Item2
						.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(i.Item1)))
						.OfType<INodalDisplacementBoundaryCondition>())
					.SelectMany(x => x)
					.OrderBy(x => x.Node.ID)
					.GroupBy(x => (x.Node.ID, x.DOF))
					.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount)))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(subdomainID))).OfType<INodalDisplacementBoundaryCondition>()
				.OrderBy(x => x.Node.ID)
				.GroupBy(x => (x.Node.ID, x.DOF))
				.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount)))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));

		public IGlobalVector GetRhs() => GetRhs(0);

		public IGlobalMatrix GetMatrix() => GetMatrix(DifferentiationOrder.Zero);
	}
}
