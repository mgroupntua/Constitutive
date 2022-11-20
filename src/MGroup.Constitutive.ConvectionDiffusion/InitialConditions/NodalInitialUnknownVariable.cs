using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.ConvectionDiffusion.InitialConditions
{
	public class NodalInitialUnknownVariable : INodalUnknownVariableInitialCondition
	{
		public IConvectionDiffusionDofType DOF { get; }

		public INode Node { get; }

		public double Amount { get; }

		public DifferentiationOrder DifferentiationOrder => DifferentiationOrder.Zero;

		public NodalInitialUnknownVariable(INode node, IConvectionDiffusionDofType dofType, double amount)
		{
			this.Node = node;
			this.DOF = dofType;
			this.Amount = amount;
		}

		INodalModelQuantity<IConvectionDiffusionDofType> INodalModelQuantity<IConvectionDiffusionDofType>.WithAmount(double amount) => new NodalInitialUnknownVariable(Node, DOF, amount);
		INodalInitialCondition<IConvectionDiffusionDofType> INodalInitialCondition<IConvectionDiffusionDofType>.WithAmount(double amount) => new NodalInitialUnknownVariable(Node, DOF, amount);
	}
}
