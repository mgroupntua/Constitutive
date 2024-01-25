using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Thermal.BoundaryConditions
{
	public class DomainTemperature : IDomainTemperatureBoundaryCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IThermalDofType DOF { get; }

		public double Multiplier { get; }

		public DomainTemperature(IThermalDofType dofType, double multiplier, DomainFunction domainFunction)
		{
			this.DOF = dofType;
			this.Multiplier = multiplier;
		}

		public DomainTemperature(IThermalDofType dofType, double multiplier)
			: this(dofType, multiplier, c => multiplier)
		{
		}

		public IDomainBoundaryCondition<IThermalDofType> WithMultiplier(double multiplier) => new DomainTemperature(DOF, multiplier, DomainFunction);
		IDomainModelQuantity<IThermalDofType> IDomainModelQuantity<IThermalDofType>.WithMultiplier(double multiplier) => new DomainTemperature(DOF, multiplier, DomainFunction);
	}
}
