using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions
{
	public class ConvectionDiffusionBoundaryConditionSet : IBoundaryConditionSet<IConvectionDiffusionDofType>
	{
		private readonly IEnumerable<INodalConvectionDiffusionDirichletBoundaryCondition> nodalDirichlet;
		private readonly IEnumerable<INodalConvectionDiffusionNeumannBoundaryCondition> nodalNeumann;
		private readonly IEnumerable<IDomainConvectionDiffusionDirichletBoundaryCondition> domainDirichlet;
		private readonly IEnumerable<IDomainConvectionDiffusionNeumannBoundaryCondition> domainNeumann;

		public ConvectionDiffusionBoundaryConditionSet(IEnumerable<INodalConvectionDiffusionDirichletBoundaryCondition> nodalDirichlet, 
			IEnumerable<INodalConvectionDiffusionNeumannBoundaryCondition> nodalNeumann, IEnumerable<IDomainConvectionDiffusionDirichletBoundaryCondition> domainDirichlet,
			IEnumerable<IDomainConvectionDiffusionNeumannBoundaryCondition> domainNeumann)
		{
			this.nodalDirichlet = nodalDirichlet;
			this.nodalNeumann = nodalNeumann;
			this.domainDirichlet = domainDirichlet;
			this.domainNeumann = domainNeumann;
		}

		public ConvectionDiffusionBoundaryConditionSet(IEnumerable<INodalConvectionDiffusionDirichletBoundaryCondition> nodalDirichlet,
			IEnumerable<INodalConvectionDiffusionNeumannBoundaryCondition> nodalNeumann) : this(nodalDirichlet, nodalNeumann, null, null)
		{
		}

		public ConvectionDiffusionBoundaryConditionSet(IEnumerable<IDomainConvectionDiffusionDirichletBoundaryCondition> domainDirichlet,
			IEnumerable<IDomainConvectionDiffusionNeumannBoundaryCondition> domainNeumann) : this(null, null, domainDirichlet, domainNeumann)
		{
		}

		public IBoundaryConditionSet<IConvectionDiffusionDofType> CreateBoundaryConditionSetOfSubdomain(ISubdomain subdomain) =>
			new ConvectionDiffusionBoundaryConditionSet(nodalDirichlet?.Where(x => x.Node.Subdomains.Contains(subdomain.ID)),
				nodalNeumann?.Where(x => x.Node.Subdomains.Contains(subdomain.ID)), domainDirichlet, domainNeumann);

		public IEnumerable<INodalBoundaryCondition<IConvectionDiffusionDofType>> EnumerateNodalBoundaryConditions() =>
			nodalDirichlet != null
			? nodalDirichlet
				.ToArray<INodalBoundaryCondition<IConvectionDiffusionDofType>>()
				.Concat(nodalNeumann != null ? nodalNeumann.ToArray<INodalBoundaryCondition<IConvectionDiffusionDofType>>() : Enumerable.Empty<INodalBoundaryCondition<IConvectionDiffusionDofType>>())
				.ToArray()
			: (nodalNeumann != null ? nodalNeumann.ToArray() : Enumerable.Empty<INodalBoundaryCondition<IConvectionDiffusionDofType>>());

		public IEnumerable<IDomainBoundaryCondition<IConvectionDiffusionDofType>> EnumerateDomainBoundaryConditions() =>
			domainDirichlet != null
			? domainDirichlet
				.ToArray<IDomainBoundaryCondition<IConvectionDiffusionDofType>>()
				.Concat(domainNeumann != null ? domainNeumann.ToArray<IDomainBoundaryCondition<IConvectionDiffusionDofType>>() : Enumerable.Empty<IDomainBoundaryCondition<IConvectionDiffusionDofType>>())
				.ToArray()
			: (domainNeumann != null ? domainNeumann.ToArray() : Enumerable.Empty<IDomainBoundaryCondition<IConvectionDiffusionDofType>>());

		public IEnumerable<INodalNeumannBoundaryCondition<IConvectionDiffusionDofType>> EnumerateEquivalentNodalNeumannBoundaryConditions(IEnumerable<IElementType> elements)
		{
			if (nodalDirichlet == null)
			{
				return Enumerable.Empty<INodalNeumannBoundaryCondition<IConvectionDiffusionDofType>>();
			}

			var equivalentNodalNeummannBoundaryConditions = new List<INodalUnknownVariableFluxBoundaryCondition>();
			var unknownVariableDirichletBoundaryConditions = nodalDirichlet.OfType<INodalUnknownVariableBoundaryCondition>();

			foreach (var element in elements.OfType<IConvectionDiffusionElementType>())
			{
				IReadOnlyList<INode> nodes = element.DofEnumerator.GetNodesForMatrixAssembly(element);
				if (nodes.Any(x => unknownVariableDirichletBoundaryConditions.Count(d => d.Node.ID == x.ID) > 0))
				{//TODO: Matrix should be Keffective? 
					var elementMatrix = element.DiffusionMatrix().Copy();
					elementMatrix.AddIntoThis(element.ConvectionMatrix());
					elementMatrix.AddIntoThis(element.ProductionMatrix());
					var nodalForces = GetNodalForcesFromMatrix(element, elementMatrix, unknownVariableDirichletBoundaryConditions);
					equivalentNodalNeummannBoundaryConditions.AddRange(nodalForces.Select(x => new NodalUnknownVariableFlux(x.Key.Node, x.Key.DOF, x.Value)));
				}
			}

			return equivalentNodalNeummannBoundaryConditions;
		}

		private IDictionary<(INode Node, IConvectionDiffusionDofType DOF), double> GetNodalForcesFromMatrix(IElementType element, IMatrix elementMatrix, IEnumerable<INodalConvectionDiffusionDirichletBoundaryCondition> boundaryConditions)
		{
			var elementNodalValues = new double[elementMatrix.NumColumns];
			element.MapNodalBoundaryConditionsToElementVector(boundaryConditions, elementNodalValues);
			var elementEquivalentForces = elementMatrix.Multiply(elementNodalValues);
			return element.MapElementVectorToNodalValues(elementEquivalentForces)
				.Where(x => x.Key.DOF is IConvectionDiffusionDofType)
				.ToDictionary(x => (x.Key.Node, (IConvectionDiffusionDofType)x.Key.DOF), x => x.Value);
		}

	}
}
