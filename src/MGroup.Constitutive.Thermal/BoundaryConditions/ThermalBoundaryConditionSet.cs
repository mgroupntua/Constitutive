using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Thermal.BoundaryConditions
{
	public class ThermalBoundaryConditionSet : IBoundaryConditionSet<IThermalDofType>
	{
		private readonly IEnumerable<INodalThermalDirichletBoundaryCondition> nodalDirichlet;
		private readonly IEnumerable<INodalThermalNeumannBoundaryCondition> nodalNeumann;
		private readonly IEnumerable<IDomainThermalDirichletBoundaryCondition> domainDirichlet;
		private readonly IEnumerable<IDomainThermalNeumannBoundaryCondition> domainNeumann;

		public ThermalBoundaryConditionSet(IEnumerable<INodalThermalDirichletBoundaryCondition> nodalDirichlet, 
			IEnumerable<INodalThermalNeumannBoundaryCondition> nodalNeumann, IEnumerable<IDomainThermalDirichletBoundaryCondition> domainDirichlet,
			IEnumerable<IDomainThermalNeumannBoundaryCondition> domainNeumann)
		{
			this.nodalDirichlet = nodalDirichlet;
			this.nodalNeumann = nodalNeumann;
			this.domainDirichlet = domainDirichlet;
			this.domainNeumann = domainNeumann;
		}

		public ThermalBoundaryConditionSet(IEnumerable<INodalThermalDirichletBoundaryCondition> nodalDirichlet,
			IEnumerable<INodalThermalNeumannBoundaryCondition> nodalNeumann) : this(nodalDirichlet, nodalNeumann, null, null)
		{
		}

		public ThermalBoundaryConditionSet(IEnumerable<IDomainThermalDirichletBoundaryCondition> domainDirichlet,
			IEnumerable<IDomainThermalNeumannBoundaryCondition> domainNeumann) : this(null, null, domainDirichlet, domainNeumann)
		{
		}

		public IBoundaryConditionSet<IThermalDofType> CreateBoundaryConditionSetOfSubdomain(ISubdomain subdomain) =>
			new ThermalBoundaryConditionSet(nodalDirichlet?.Where(x => x.Node.Subdomains.Contains(subdomain.ID)),
				nodalNeumann?.Where(x => x.Node.Subdomains.Contains(subdomain.ID)), domainDirichlet, domainNeumann);

		public IEnumerable<INodalBoundaryCondition<IThermalDofType>> EnumerateNodalBoundaryConditions() =>
			nodalDirichlet != null
			? nodalDirichlet
				.ToArray<INodalBoundaryCondition<IThermalDofType>>()
				.Concat(nodalNeumann != null ? nodalNeumann.ToArray<INodalBoundaryCondition<IThermalDofType>>() : Enumerable.Empty<INodalBoundaryCondition<IThermalDofType>>())
				.ToArray()
			: (nodalNeumann != null ? nodalNeumann.ToArray() : Enumerable.Empty<INodalBoundaryCondition<IThermalDofType>>());

		public IEnumerable<IDomainBoundaryCondition<IThermalDofType>> EnumerateDomainBoundaryConditions() =>
			domainDirichlet != null
			? domainDirichlet
				.ToArray<IDomainBoundaryCondition<IThermalDofType>>()
				.Concat(domainNeumann != null ? domainNeumann.ToArray<IDomainBoundaryCondition<IThermalDofType>>() : Enumerable.Empty<IDomainBoundaryCondition<IThermalDofType>>())
				.ToArray()
			: (domainNeumann != null ? domainNeumann.ToArray() : Enumerable.Empty<IDomainBoundaryCondition<IThermalDofType>>());

		public IEnumerable<INodalNeumannBoundaryCondition<IThermalDofType>> EnumerateEquivalentNodalNeumannBoundaryConditions(IEnumerable<IElementType> elements)
		{
			if (nodalDirichlet == null)
			{
				return Enumerable.Empty<INodalNeumannBoundaryCondition<IThermalDofType>>();
			}

			var equivalentNodalNeummannBoundaryConditions = new List<INodalHeatFluxBoundaryCondition>();
			var temperatureBoundaryConditions = nodalDirichlet.OfType<INodalTemperatureBoundaryCondition>();

			foreach (var element in elements.OfType<IThermalElementType>())
			{
				IReadOnlyList<INode> nodes = element.DofEnumerator.GetNodesForMatrixAssembly(element);
				if (nodes.Any(x => temperatureBoundaryConditions.Count(d => d.Node.ID == x.ID) > 0))
				{
					var elementMatrix = element.ConductivityMatrix();
					var nodalForces = GetNodalForcesFromMatrix(element, elementMatrix, temperatureBoundaryConditions);
					equivalentNodalNeummannBoundaryConditions.AddRange(nodalForces.Select(x => new NodalHeatFlux(x.Key.Node, x.Key.DOF, x.Value)));
				}
			}

			return equivalentNodalNeummannBoundaryConditions;
		}

		private IDictionary<(INode Node, IThermalDofType DOF), double> GetNodalForcesFromMatrix(IElementType element, IMatrix elementMatrix, IEnumerable<INodalThermalDirichletBoundaryCondition> boundaryConditions)
		{
			var elementNodalValues = new double[elementMatrix.NumColumns];
			element.MapNodalBoundaryConditionsToElementVector(boundaryConditions, elementNodalValues);
			var elementEquivalentForces = elementMatrix.Multiply(elementNodalValues);
			return element.MapElementVectorToNodalValues(elementEquivalentForces)
				.Where(x => x.Key.DOF is IThermalDofType)
				.ToDictionary(x => (x.Key.Node, (IThermalDofType)x.Key.DOF), x => x.Value);
		}

	}
}
