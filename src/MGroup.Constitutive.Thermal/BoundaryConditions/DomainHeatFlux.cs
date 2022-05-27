using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Thermal.BoundaryConditions
{
	public class DomainHeatFlux : IDomainHeatFluxBoundaryCondition
	{
		public IThermalDofType DOF { get; }

		public double Amount { get; }

		public DomainHeatFlux(IThermalDofType dof, double amount)
		{
			DOF = dof;
			Amount = amount;
		}

		public IDomainBoundaryCondition<IThermalDofType> WithAmount(double amount) => new DomainHeatFlux(DOF, amount);
	}
}
