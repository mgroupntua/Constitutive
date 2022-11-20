using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization;

namespace MGroup.Constitutive.Structural.InitialConditions
{
	public class DomainInitialDisplacement : IDomainDisplacementInitialCondition
	{
		public IStructuralDofType DOF { get; }

		public double Amount { get; }

		public DifferentiationOrder DifferentiationOrder => DifferentiationOrder.Zero;

		public DomainInitialDisplacement(IStructuralDofType dofType, double amount)
		{
			this.DOF = dofType;
			this.Amount = amount;
		}

		IDomainInitialCondition<IStructuralDofType> IDomainInitialCondition<IStructuralDofType>.WithAmount(double amount) => new DomainInitialDisplacement(DOF, amount);
		IDomainModelQuantity<IStructuralDofType> IDomainModelQuantity<IStructuralDofType>.WithAmount(double amount) => new DomainInitialDisplacement(DOF, amount);
	}
}
