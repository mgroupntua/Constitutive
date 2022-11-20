using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Structural.InitialConditions
{
	public class NodalInitialVelocity : INodalVelocityInitialCondition
	{
		public IStructuralDofType DOF { get; }

		public INode Node { get; }

		public double Amount { get; }

		public DifferentiationOrder DifferentiationOrder => DifferentiationOrder.First;

		public NodalInitialVelocity(INode node, IStructuralDofType dofType, double amount)
		{
			this.Node = node;
			this.DOF = dofType;
			this.Amount = amount;
		}

		INodalModelQuantity<IStructuralDofType> INodalModelQuantity<IStructuralDofType>.WithAmount(double amount) => new NodalInitialVelocity(Node, DOF, amount);
		public INodalInitialCondition<IStructuralDofType> WithAmount(double amount) => throw new System.NotImplementedException();
	}
}
