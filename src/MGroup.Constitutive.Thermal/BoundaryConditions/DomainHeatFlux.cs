using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Thermal.BoundaryConditions
{
	public class DomainHeatFlux : IDomainHeatFluxBoundaryCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IThermalDofType DOF { get; }

		public double Multiplier { get; }

		public DomainHeatFlux(IThermalDofType dof, double multiplier, DomainFunction domainFunction)
		{
			DOF = dof;
			Multiplier = multiplier;
			DomainFunction = domainFunction;
		}

		public DomainHeatFlux(IThermalDofType dof, double multiplier)
			: this(dof, multiplier, c => multiplier)
		{
		}

		public IDomainBoundaryCondition<IThermalDofType> WithMultiplier(double multiplier) => new DomainHeatFlux(DOF, multiplier, DomainFunction);
		IDomainModelQuantity<IThermalDofType> IDomainModelQuantity<IThermalDofType>.WithMultiplier(double multiplier) => new DomainHeatFlux(DOF, multiplier, DomainFunction);
	}
}
