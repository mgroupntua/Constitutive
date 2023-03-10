using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;

namespace MGroup.Constitutive.ConvectionDiffusion.InitialConditions
{
	public class DomainInitialUnknownVariable : IDomainUnknownVariableInitialCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IConvectionDiffusionDofType DOF { get; }

		public double Multiplier { get; }

		public DifferentiationOrder DifferentiationOrder => DifferentiationOrder.Zero;

		public DomainInitialUnknownVariable(IConvectionDiffusionDofType dofType, double multiplier, DomainFunction domainFunction)
		{
			this.DOF = dofType;
			this.Multiplier = multiplier;
			this.DomainFunction = domainFunction;
		}

		IDomainModelQuantity<IConvectionDiffusionDofType> IDomainModelQuantity<IConvectionDiffusionDofType>.WithMultiplier(double multiplier) => new DomainInitialUnknownVariable(DOF, multiplier, DomainFunction);
		IDomainInitialCondition<IConvectionDiffusionDofType> IDomainInitialCondition<IConvectionDiffusionDofType>.WithMultiplier(double multiplier) => new DomainInitialUnknownVariable(DOF, multiplier, DomainFunction);
	}
}
