using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;

namespace MGroup.Constitutive.Thermal.InitialConditions
{
	public class DomainInitialTemperature : IDomainTemperatureInitialCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IThermalDofType DOF { get; }

		public double Multiplier { get; }

		public DifferentiationOrder DifferentiationOrder => DifferentiationOrder.Zero;

		public DomainInitialTemperature(IThermalDofType dofType, double multiplier, DomainFunction domainFunction)
		{
			this.DOF = dofType;
			this.Multiplier = multiplier;
		}

		public DomainInitialTemperature(IThermalDofType dofType, double multiplier)
		{
		}

		IDomainModelQuantity<IThermalDofType> IDomainModelQuantity<IThermalDofType>.WithMultiplier(double multiplier) => new DomainInitialTemperature(DOF, multiplier, DomainFunction);
		IDomainInitialCondition<IThermalDofType> IDomainInitialCondition<IThermalDofType>.WithMultiplier(double multiplier) => new DomainInitialTemperature(DOF, multiplier, DomainFunction);
	}
}
