using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Thermal.BoundaryConditions
{
	public class NodalTemperature : INodalTemperatureBoundaryCondition
	{
		public IThermalDofType DOF { get; }

		public INode Node { get; }

		public double Amount { get; }

		public NodalTemperature(INode node, IThermalDofType dofType, double amount)
		{
			this.Node = node;
			this.DOF = dofType;
			this.Amount = amount;
		}

		public INodalBoundaryCondition<IThermalDofType> WithAmount(double amount) => new NodalTemperature(Node, DOF, amount);
	}
}
