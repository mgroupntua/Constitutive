using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Thermal.BoundaryConditions
{
	public class DomainTemperature : IDomainTemperatureBoundaryCondition
	{
		public IThermalDofType DOF { get; }

		public double Amount { get; }

		public DomainTemperature(IThermalDofType dofType, double amount)
		{
			this.DOF = dofType;
			this.Amount = amount;
		}

		public IDomainBoundaryCondition<IThermalDofType> WithAmount(double amount) => new DomainTemperature(DOF, amount);
	}
}
