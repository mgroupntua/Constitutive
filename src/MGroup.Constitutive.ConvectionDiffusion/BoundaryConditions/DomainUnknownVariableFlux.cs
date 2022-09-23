using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions
{
	public class DomainUnknownVariableFlux : IDomainUnknownVariableFluxBoundaryCondition
	{
		public IConvectionDiffusionDofType DOF { get; }

		public double Amount { get; }

		public DomainUnknownVariableFlux(IConvectionDiffusionDofType dof, double amount)
		{
			DOF = dof;
			Amount = amount;
		}

		public IDomainBoundaryCondition<IConvectionDiffusionDofType> WithAmount(double amount) => new DomainUnknownVariableFlux(DOF, amount);
	}
}
