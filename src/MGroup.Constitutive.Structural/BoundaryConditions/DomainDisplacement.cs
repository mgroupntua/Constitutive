using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public class DomainDisplacement : IDomainDisplacementBoundaryCondition
	{
		public IStructuralDofType DOF { get; }

		public double Amount { get; }

		public DomainDisplacement(IStructuralDofType dofType, double amount)
		{
			this.DOF = dofType;
			this.Amount = amount;
		}

		public IDomainBoundaryCondition<IStructuralDofType> WithAmount(double amount) => new DomainDisplacement(DOF, amount);
	}
}
