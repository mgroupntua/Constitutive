using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Thermal.InitialConditions
{
	public class NodalInitialTemperature : INodalTemperatureInitialCondition
	{
		public IThermalDofType DOF { get; }

		public INode Node { get; }

		public double Amount { get; }

		public DifferentiationOrder DifferentiationOrder => DifferentiationOrder.Zero;

		public NodalInitialTemperature(INode node, IThermalDofType dofType, double amount)
		{
			this.Node = node;
			this.DOF = dofType;
			this.Amount = amount;
		}

		INodalModelQuantity<IThermalDofType> INodalModelQuantity<IThermalDofType>.WithAmount(double amount) => new NodalInitialTemperature(Node, DOF, amount);
		INodalInitialCondition<IThermalDofType> INodalInitialCondition<IThermalDofType>.WithAmount(double amount) => new NodalInitialTemperature(Node, DOF, amount);
	}
}
