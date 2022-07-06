using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions
{
	public class NodalUnknownVariable : INodalUnknownVariableBoundaryCondition
	{
		public IConvectionDiffusionDofType DOF { get; }

		public INode Node { get; }

		public double Amount { get; }

		public NodalUnknownVariable(INode node, IConvectionDiffusionDofType dofType, double amount)
		{
			this.Node = node;
			this.DOF = dofType;
			this.Amount = amount;
		}

		public INodalBoundaryCondition<IConvectionDiffusionDofType> WithAmount(double amount) => new NodalUnknownVariable(Node, DOF, amount);
	}
}
