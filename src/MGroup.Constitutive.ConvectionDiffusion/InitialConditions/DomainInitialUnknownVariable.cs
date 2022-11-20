using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization;

namespace MGroup.Constitutive.ConvectionDiffusion.InitialConditions
{
	public class DomainInitialUnknownVariable : IDomainUnknownVariableInitialCondition
	{
		public IConvectionDiffusionDofType DOF { get; }

		public double Amount { get; }

		public DifferentiationOrder DifferentiationOrder => DifferentiationOrder.Zero;

		public DomainInitialUnknownVariable(IConvectionDiffusionDofType dofType, double amount)
		{
			this.DOF = dofType;
			this.Amount = amount;
		}

		IDomainModelQuantity<IConvectionDiffusionDofType> IDomainModelQuantity<IConvectionDiffusionDofType>.WithAmount(double amount) => new DomainInitialUnknownVariable(DOF, amount);
		IDomainInitialCondition<IConvectionDiffusionDofType> IDomainInitialCondition<IConvectionDiffusionDofType>.WithAmount(double amount) => new DomainInitialUnknownVariable(DOF, amount);
	}
}
