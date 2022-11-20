using System.Collections.Generic;
using System.Linq;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Thermal.InitialConditions
{
	public class ThermalInitialConditionSet : IInitialConditionSet<IThermalDofType>
	{
		private readonly IEnumerable<INodalThermalInitialCondition> nodalDirichlet;
		private readonly IEnumerable<IDomainThermalInitialCondition> domainDirichlet;

		public ThermalInitialConditionSet(IEnumerable<INodalThermalInitialCondition> nodalDirichlet,
			IEnumerable<IDomainThermalInitialCondition> domainDirichlet)
		{
			this.nodalDirichlet = nodalDirichlet;
			this.domainDirichlet = domainDirichlet;
		}

		public IInitialConditionSet<IThermalDofType> CreateInitialConditionSetOfSubdomain(ISubdomain subdomain) =>
			new ThermalInitialConditionSet(nodalDirichlet?.Where(x => x.Node.Subdomains.Contains(subdomain.ID)), domainDirichlet);

		public IEnumerable<INodalInitialCondition<IThermalDofType>> EnumerateNodalInitialConditions() =>
			nodalDirichlet != null
			? nodalDirichlet.ToArray<INodalInitialCondition<IThermalDofType>>()
			: Enumerable.Empty<INodalInitialCondition<IThermalDofType>>();

		public IEnumerable<IDomainInitialCondition<IThermalDofType>> EnumerateDomainInitialConditions() =>
			domainDirichlet != null
			? domainDirichlet.ToArray<IDomainInitialCondition<IThermalDofType>>()
			: Enumerable.Empty<IDomainInitialCondition<IThermalDofType>>();
	}
}
