using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Structural.InitialConditions
{
	public class StructuralInitialConditionSet : IInitialConditionSet<IStructuralDofType>
	{
		private readonly IEnumerable<INodalStructuralInitialCondition> nodalDirichlet;
		private readonly IEnumerable<IDomainStructuralInitialCondition> domainDirichlet;

		public StructuralInitialConditionSet(IEnumerable<INodalStructuralInitialCondition> nodalDirichlet, 
			IEnumerable<IDomainStructuralInitialCondition> domainDirichlet)
		{
			this.nodalDirichlet = nodalDirichlet;
			this.domainDirichlet = domainDirichlet;
		}

		public IInitialConditionSet<IStructuralDofType> CreateInitialConditionSetOfSubdomain(ISubdomain subdomain) => 
			new StructuralInitialConditionSet(nodalDirichlet?.Where(x => x.Node.Subdomains.Contains(subdomain.ID)), domainDirichlet);

		public IEnumerable<INodalInitialCondition<IStructuralDofType>> EnumerateNodalInitialConditions() =>
			nodalDirichlet != null 
				? nodalDirichlet.ToArray<INodalInitialCondition<IStructuralDofType>>() 
				: Enumerable.Empty<INodalInitialCondition<IStructuralDofType>>();

		public IEnumerable<IDomainInitialCondition<IStructuralDofType>> EnumerateDomainInitialConditions() =>
			domainDirichlet != null
			? domainDirichlet.ToArray<IDomainInitialCondition<IStructuralDofType>>()
			: Enumerable.Empty<IDomainInitialCondition<IStructuralDofType>>();
	}
}
