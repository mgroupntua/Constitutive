using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public class DomainLoad : IDomainLoadBoundaryCondition
	{
		public IStructuralDofType DOF { get; }

		public double Amount { get; }

		public DomainLoad(IStructuralDofType dof, double amount)
		{
			DOF = dof;
			Amount = amount;
		}

		public IDomainBoundaryCondition<IStructuralDofType> WithAmount(double amount) => new DomainLoad(DOF, amount);
	}
}
