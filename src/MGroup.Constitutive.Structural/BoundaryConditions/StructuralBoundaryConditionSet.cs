using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public class StructuralBoundaryConditionSet : IBoundaryConditionSet<IStructuralDofType>
	{
		private readonly IEnumerable<INodalStructuralDirichletBoundaryCondition> nodalDirichlet;
		private readonly IEnumerable<INodalStructuralNeumannBoundaryCondition> nodalNeumann;
		private readonly IEnumerable<IDomainStructuralDirichletBoundaryCondition> domainDirichlet;
		private readonly IEnumerable<IDomainStructuralNeumannBoundaryCondition> domainNeumann;

		public StructuralBoundaryConditionSet(IEnumerable<INodalStructuralDirichletBoundaryCondition> nodalDirichlet, 
			IEnumerable<INodalStructuralNeumannBoundaryCondition> nodalNeumann, IEnumerable<IDomainStructuralDirichletBoundaryCondition> domainDirichlet,
			IEnumerable<IDomainStructuralNeumannBoundaryCondition> domainNeumann)
		{
			this.nodalDirichlet = nodalDirichlet;
			this.nodalNeumann = nodalNeumann;
			this.domainNeumann = domainNeumann;
			this.domainDirichlet = domainDirichlet;
		}

		public StructuralBoundaryConditionSet(IEnumerable<INodalStructuralDirichletBoundaryCondition> nodalDirichlet, 
			IEnumerable<INodalStructuralNeumannBoundaryCondition> nodalNeumann)
			: this(nodalDirichlet, nodalNeumann, null, null)
		{
		}

		public StructuralBoundaryConditionSet(IEnumerable<IDomainStructuralDirichletBoundaryCondition> domainDirichlet,
			IEnumerable<IDomainStructuralNeumannBoundaryCondition> domainNeumann)
			: this(null, null, domainDirichlet, domainNeumann)
		{
		}

		public IBoundaryConditionSet<IStructuralDofType> CreateBoundaryConditionSetOfSubdomain(ISubdomain subdomain) => 
			new StructuralBoundaryConditionSet(nodalDirichlet?.Where(x => x.Node.Subdomains.Contains(subdomain.ID)), 
				nodalNeumann?.Where(x => x.Node.Subdomains.Contains(subdomain.ID)), domainDirichlet, domainNeumann);

		public IEnumerable<INodalBoundaryCondition<IStructuralDofType>> EnumerateNodalBoundaryConditions() =>
			nodalDirichlet != null
			? nodalDirichlet
				.ToArray<INodalBoundaryCondition<IStructuralDofType>>()
				.Concat(nodalNeumann != null ? nodalNeumann.ToArray<INodalBoundaryCondition<IStructuralDofType>>() : Enumerable.Empty<INodalBoundaryCondition<IStructuralDofType>>())
				.ToArray()
			: (nodalNeumann != null ? nodalNeumann.ToArray() : Enumerable.Empty<INodalBoundaryCondition<IStructuralDofType>>());

		public IEnumerable<IDomainBoundaryCondition<IStructuralDofType>> EnumerateDomainBoundaryConditions() =>
			domainDirichlet != null
			? domainDirichlet
				.ToArray<IDomainBoundaryCondition<IStructuralDofType>>()
				.Concat(domainNeumann != null ? domainNeumann.ToArray<IDomainBoundaryCondition<IStructuralDofType>>() : Enumerable.Empty<IDomainBoundaryCondition<IStructuralDofType>>())
				.ToArray()
			: (domainNeumann != null ? domainNeumann.ToArray() : Enumerable.Empty<IDomainBoundaryCondition<IStructuralDofType>>());

		public IEnumerable<INodalNeumannBoundaryCondition<IStructuralDofType>> EnumerateEquivalentNodalNeumannBoundaryConditions(IEnumerable<IElementType> elements)
		{
			if (nodalDirichlet == null)
			{
				return Enumerable.Empty<INodalNeumannBoundaryCondition<IStructuralDofType>>();
			}

			var equivalentNodalNeummannBoundaryConditions = new List<INodalLoadBoundaryCondition>();
			var displacementBoundaryConditions = nodalDirichlet.OfType<INodalDisplacementBoundaryCondition>();

			foreach (var element in elements.OfType<IStructuralElementType>())
			{
				IReadOnlyList<INode> nodes = element.DofEnumerator.GetNodesForMatrixAssembly(element);
				if (nodes.Any(x => displacementBoundaryConditions.Count(d => d.Node.ID == x.ID) > 0))
				{
					var elementMatrix = element.StiffnessMatrix();
					var nodalForces = GetNodalForcesFromMatrix(element, elementMatrix, displacementBoundaryConditions);
					equivalentNodalNeummannBoundaryConditions.AddRange(nodalForces.Select(x => new NodalLoad(x.Key.Node, x.Key.DOF, x.Value)));
				}
			}

			return equivalentNodalNeummannBoundaryConditions;
		}

		private IDictionary<(INode Node, IStructuralDofType DOF), double> GetNodalForcesFromMatrix(IElementType element, IMatrix elementMatrix, IEnumerable<INodalStructuralDirichletBoundaryCondition> boundaryConditions)
		{
			var elementNodalValues = new double[elementMatrix.NumColumns];
			element.MapNodalBoundaryConditionsToElementVector(boundaryConditions, elementNodalValues);
			var elementEquivalentForces = elementMatrix.Multiply(elementNodalValues);
			return element.MapElementVectorToNodalValues(elementEquivalentForces)
				.Where(x => x.Key.DOF is IStructuralDofType)
				.ToDictionary(x => (x.Key.Node, (IStructuralDofType)x.Key.DOF), x => x.Value);
		}

	}
}
