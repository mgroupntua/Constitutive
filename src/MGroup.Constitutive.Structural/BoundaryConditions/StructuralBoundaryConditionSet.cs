using System;
using System.Collections.Generic;
using System.Linq;
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
		private readonly IEnumerable<IElementStructuralDirichletBoundaryCondition> elementDirichlet;
		private readonly IEnumerable<IElementStructuralNeumannBoundaryCondition> elementNeumann;
		private readonly IEnumerable<IDomainStructuralDirichletBoundaryCondition> domainDirichlet;
		private readonly IEnumerable<IDomainStructuralNeumannBoundaryCondition> domainNeumann;

		public StructuralBoundaryConditionSet(IEnumerable<INodalStructuralDirichletBoundaryCondition> nodalDirichlet, 
			IEnumerable<INodalStructuralNeumannBoundaryCondition> nodalNeumann, IEnumerable<IElementStructuralDirichletBoundaryCondition> elementDirichlet,
			IEnumerable<IElementStructuralNeumannBoundaryCondition> elementNeumann, IEnumerable<IDomainStructuralDirichletBoundaryCondition> domainDirichlet,
			IEnumerable<IDomainStructuralNeumannBoundaryCondition> domainNeumann)
		{
			this.nodalDirichlet = nodalDirichlet;
			this.nodalNeumann = nodalNeumann;
			this.elementDirichlet = elementDirichlet;
			this.elementNeumann = elementNeumann;
			this.domainNeumann = domainNeumann;
			this.domainDirichlet = domainDirichlet;
		}

		public StructuralBoundaryConditionSet(IEnumerable<INodalStructuralDirichletBoundaryCondition> nodalDirichlet, 
			IEnumerable<INodalStructuralNeumannBoundaryCondition> nodalNeumann)
			: this(nodalDirichlet, nodalNeumann, null, null, null, null)
		{
		}

		public StructuralBoundaryConditionSet(IEnumerable<IDomainStructuralDirichletBoundaryCondition> domainDirichlet,
			IEnumerable<IDomainStructuralNeumannBoundaryCondition> domainNeumann)
			: this(null, null, null, null, domainDirichlet, domainNeumann)
		{
		}

		public IBoundaryConditionSet<IStructuralDofType> CreateBoundaryConditionSetOfSubdomain(ISubdomain subdomain) => 
			new StructuralBoundaryConditionSet(nodalDirichlet?.Where(x => x.Node.Subdomains.Contains(subdomain.ID)), 
				nodalNeumann?.Where(x => x.Node.Subdomains.Contains(subdomain.ID)), elementDirichlet?.Where(x => x.Element.SubdomainID.Equals(subdomain.ID)),
				elementNeumann?.Where(x => x.Element.SubdomainID.Equals(subdomain.ID)), domainDirichlet, domainNeumann);

		public IEnumerable<INodalBoundaryCondition<IStructuralDofType>> EnumerateNodalBoundaryConditions(IEnumerable<IElementType> elements)
		{
			if (domainDirichlet != null)
			{
				throw new NotImplementedException("Dirichlet boundary conditions for domain are not implemented");
			}

			if (domainNeumann != null)
			{
				throw new NotImplementedException("Neumann boundary conditions for domain are not implemented");
			}

			if (elementDirichlet != null) 
			{
				throw new NotImplementedException("Dirichlet boundary conditions for elements are not implemented");
			}

			if (elementNeumann != null && elementNeumann.Any(x => x is IElementDistributedLoadBoundaryCondition == false))
			{
				throw new NotImplementedException("Neumann boundary conditions for elements that are NOT distributed loads, are not implemented");
			}

			var nodalQuantities = nodalDirichlet != null
			? nodalDirichlet
				.ToArray<INodalBoundaryCondition<IStructuralDofType>>()
				.Concat(nodalNeumann != null ? nodalNeumann.ToArray<INodalBoundaryCondition<IStructuralDofType>>() : Enumerable.Empty<INodalBoundaryCondition<IStructuralDofType>>())
				.ToArray()
			: (nodalNeumann != null ? nodalNeumann.ToArray() : Enumerable.Empty<INodalBoundaryCondition<IStructuralDofType>>());

			var allQuantities = new List<INodalBoundaryCondition<IStructuralDofType>>(nodalQuantities);
			if (elementNeumann != null)
			{
				var elementNeumanDictionary = elementNeumann.GroupBy(x => x.Element.ID).ToDictionary(x => x.Key, x => x.ToArray());
				foreach (var element in elements.OfType<IStructuralElementType>().Where(x => elementNeumanDictionary.Keys.Any(e => e == x.ID)))
				{
					IReadOnlyList<INode> nodes = element.DofEnumerator.GetNodesForMatrixAssembly(element);
					var elementValues = element.IntegrateElementModelQuantities(elementNeumanDictionary[element.ID])
						.Aggregate(new double[element.GetElementDofTypes().Aggregate(0, (a, x) => x.Count + a)], (a, x) => Enumerable.Zip(a, x, (aE, xE) => aE + xE).ToArray());
					allQuantities.AddRange(element.MapElementVectorToNodalValues(elementValues)
						.Where(x => x.Key.DOF is IStructuralDofType)
						.Select(x => new NodalLoad(x.Key.Node, (IStructuralDofType)x.Key.DOF, x.Value)));
				}
			}

			return allQuantities;
		}

		//public IEnumerable<IDomainBoundaryCondition<IStructuralDofType>> EnumerateDomainBoundaryConditions() =>
		//	domainDirichlet != null
		//	? domainDirichlet
		//		.ToArray<IDomainBoundaryCondition<IStructuralDofType>>()
		//		.Concat(domainNeumann != null ? domainNeumann.ToArray<IDomainBoundaryCondition<IStructuralDofType>>() : Enumerable.Empty<IDomainBoundaryCondition<IStructuralDofType>>())
		//		.ToArray()
		//	: (domainNeumann != null ? domainNeumann.ToArray() : Enumerable.Empty<IDomainBoundaryCondition<IStructuralDofType>>());

		public IEnumerable<INodalNeumannBoundaryCondition<IStructuralDofType>> EnumerateEquivalentNodalNeumannBoundaryConditions(IEnumerable<IElementType> elements, 
			IEnumerable<(int NodeID, IDofType DOF)> dofsToExclude)
		{
			if (domainDirichlet != null)
			{
				throw new NotImplementedException("Dirichlet boundary conditions for domain are not implemented");
			}

			if (elementDirichlet != null)
			{
				throw new NotImplementedException("Dirichlet boundary conditions for elements are not implemented");
			}

			if (nodalDirichlet == null)
			{
				return Enumerable.Empty<INodalNeumannBoundaryCondition<IStructuralDofType>>();
			}

			var equivalentNodalNeummannBoundaryConditions = new List<INodalLoadBoundaryCondition>();
			var displacementBoundaryConditions = nodalDirichlet.OfType<INodalDisplacementBoundaryCondition>()
				.Where(x => dofsToExclude.Any(d => d.NodeID == x.Node.ID && d.DOF == x.DOF) == false);

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
