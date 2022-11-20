using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization;

namespace MGroup.Constitutive.Structural.InitialConditions
{
	public class DomainInitialVelocity : IDomainVelocityInitialCondition
	{
		public IStructuralDofType DOF { get; }

		public double Amount { get; }

		public DifferentiationOrder DifferentiationOrder => DifferentiationOrder.First;

		public DomainInitialVelocity(IStructuralDofType dofType, double amount)
		{
			this.DOF = dofType;
			this.Amount = amount;
		}

		IDomainModelQuantity<IStructuralDofType> IDomainModelQuantity<IStructuralDofType>.WithAmount(double amount) => new DomainInitialVelocity(DOF, amount);
		IDomainInitialCondition<IStructuralDofType> IDomainInitialCondition<IStructuralDofType>.WithAmount(double amount) => new DomainInitialVelocity(DOF, amount);
	}
}
