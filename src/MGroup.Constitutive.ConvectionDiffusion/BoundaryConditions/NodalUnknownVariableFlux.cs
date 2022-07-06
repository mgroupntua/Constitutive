using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions
{
	public class NodalUnknownVariableFlux : INodalUnknownVariableFluxBoundaryCondition
	{
		public INode Node { get; }

		public IConvectionDiffusionDofType DOF { get; }

		public double Amount { get; }

		public NodalUnknownVariableFlux(INode node, IConvectionDiffusionDofType dof, double amount)
		{
			Node = node;
			DOF = dof;
			Amount = amount;
		}

		public INodalBoundaryCondition<IConvectionDiffusionDofType> WithAmount(double amount) => new NodalUnknownVariableFlux(Node, DOF, amount);
	}
}
