using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public class DomainVelocity : IDomainVelocityBoundaryCondition
	{
		public IStructuralDofType DOF { get; }

		public double Amount { get; }

		public DomainVelocity(IStructuralDofType dofType, double amount)
		{
			this.DOF = dofType;
			this.Amount = amount;
		}

		public IDomainBoundaryCondition<IStructuralDofType> WithAmount(double amount) => new DomainVelocity(DOF, amount);
	}
}
