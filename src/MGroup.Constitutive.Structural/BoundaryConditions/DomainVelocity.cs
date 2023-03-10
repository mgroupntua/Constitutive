using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public class DomainVelocity : IDomainVelocityBoundaryCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IStructuralDofType DOF { get; }

		public double Multiplier { get; }

		public DomainVelocity(IStructuralDofType dofType, double multiplier, DomainFunction domainFunction)
		{
			this.DOF = dofType;
			this.Multiplier = multiplier;
			this.DomainFunction = domainFunction;
		}

		public DomainVelocity(IStructuralDofType dofType, double multiplier)
			: this(dofType, multiplier, c => multiplier)
		{
		}

		public IDomainBoundaryCondition<IStructuralDofType> WithMultiplier(double multiplier) => new DomainVelocity(DOF, multiplier, DomainFunction);
		IDomainModelQuantity<IStructuralDofType> IDomainModelQuantity<IStructuralDofType>.WithMultiplier(double multiplier) => new DomainVelocity(DOF, multiplier, DomainFunction);
	}
}
