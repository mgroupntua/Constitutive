using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions
{
	public class DomainUnknownVariable : IDomainUnknownVariableBoundaryCondition
	{
		public IConvectionDiffusionDofType DOF { get; }

		public double Amount { get; }

		public DomainUnknownVariable(IConvectionDiffusionDofType dofType, double amount)
		{
			this.DOF = dofType;
			this.Amount = amount;
		}

		public IDomainBoundaryCondition<IConvectionDiffusionDofType> WithAmount(double amount) => new DomainUnknownVariable(DOF, amount);
		IDomainModelQuantity<IConvectionDiffusionDofType> IDomainModelQuantity<IConvectionDiffusionDofType>.WithAmount(double amount) => new DomainUnknownVariable(DOF, amount);
	}
}
