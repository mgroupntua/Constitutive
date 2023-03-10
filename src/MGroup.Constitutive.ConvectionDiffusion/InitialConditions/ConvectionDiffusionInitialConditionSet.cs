using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.ConvectionDiffusion.InitialConditions
{
	public class ConvectionDiffusionInitialConditionSet : IInitialConditionSet<IConvectionDiffusionDofType>
	{
		private readonly IEnumerable<INodalConvectionDiffusionInitialCondition> nodalDirichlet;
		private readonly IEnumerable<IDomainConvectionDiffusionInitialCondition> domainDirichlet;

		public ConvectionDiffusionInitialConditionSet(IEnumerable<INodalConvectionDiffusionInitialCondition> nodalDirichlet, 
			IEnumerable<IDomainConvectionDiffusionInitialCondition> domainDirichlet)
		{
			this.nodalDirichlet = nodalDirichlet;
			this.domainDirichlet = domainDirichlet;
		}

		public IInitialConditionSet<IConvectionDiffusionDofType> CreateInitialConditionSetOfSubdomain(ISubdomain subdomain) =>
			new ConvectionDiffusionInitialConditionSet(nodalDirichlet?.Where(x => x.Node.Subdomains.Contains(subdomain.ID)), domainDirichlet);

		public IEnumerable<INodalInitialCondition<IConvectionDiffusionDofType>> EnumerateNodalInitialConditions(IEnumerable<IElementType> elements) =>
			nodalDirichlet != null
			? nodalDirichlet
				.ToArray<INodalInitialCondition<IConvectionDiffusionDofType>>().ToArray()
			: Enumerable.Empty<INodalInitialCondition<IConvectionDiffusionDofType>>();

		//public IEnumerable<IDomainInitialCondition<IConvectionDiffusionDofType>> EnumerateDomainInitialConditions() =>
		//	domainDirichlet != null
		//	? domainDirichlet
		//		.ToArray<IDomainInitialCondition<IConvectionDiffusionDofType>>()
		//		.ToArray()
		//	: Enumerable.Empty<IDomainInitialCondition<IConvectionDiffusionDofType>>();
	}
}
