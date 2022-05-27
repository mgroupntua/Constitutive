using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public class NodalLoad : INodalLoadBoundaryCondition
	{
		public INode Node { get; }

		public IStructuralDofType DOF { get; }

		public double Amount { get; }

		public NodalLoad(INode node, IStructuralDofType dof, double amount)
		{
			Node = node;
			DOF = dof;
			Amount = amount;
		}

		public INodalBoundaryCondition<IStructuralDofType> WithAmount(double amount) => new NodalLoad(Node, DOF, amount);
	}
}
