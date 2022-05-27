using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Thermal.BoundaryConditions
{
	public class NodalHeatFlux : INodalHeatFluxBoundaryCondition
	{
		public INode Node { get; }

		public IThermalDofType DOF { get; }

		public double Amount { get; }

		public NodalHeatFlux(INode node, IThermalDofType dof, double amount)
		{
			Node = node;
			DOF = dof;
			Amount = amount;
		}

		public INodalBoundaryCondition<IThermalDofType> WithAmount(double amount) => new NodalHeatFlux(Node, DOF, amount);
	}
}
