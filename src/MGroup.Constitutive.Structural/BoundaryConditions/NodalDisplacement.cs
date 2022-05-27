using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public class NodalDisplacement : INodalDisplacementBoundaryCondition
	{
		public IStructuralDofType DOF { get; }

		public INode Node { get; }

		public double Amount { get; }

		public NodalDisplacement(INode node, IStructuralDofType dofType, double amount)
		{
			this.Node = node;
			this.DOF = dofType;
			this.Amount = amount;
		}

		public INodalBoundaryCondition<IStructuralDofType> WithAmount(double amount) => new NodalDisplacement(Node, DOF, amount);
	}
}
