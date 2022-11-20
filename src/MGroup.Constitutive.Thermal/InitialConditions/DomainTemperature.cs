using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization;

namespace MGroup.Constitutive.Thermal.InitialConditions
{
	public class DomainInitialTemperature : IDomainTemperatureInitialCondition
	{
		public IThermalDofType DOF { get; }

		public double Amount { get; }

		public DifferentiationOrder DifferentiationOrder => DifferentiationOrder.Zero;

		public DomainInitialTemperature(IThermalDofType dofType, double amount)
		{
			this.DOF = dofType;
			this.Amount = amount;
		}

		IDomainModelQuantity<IThermalDofType> IDomainModelQuantity<IThermalDofType>.WithAmount(double amount) => new DomainInitialTemperature(DOF, amount);
		IDomainInitialCondition<IThermalDofType> IDomainInitialCondition<IThermalDofType>.WithAmount(double amount) => new DomainInitialTemperature(DOF, amount);
	}
}
