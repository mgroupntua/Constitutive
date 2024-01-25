using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public class DomainAcceleration : IDomainAccelerationBoundaryCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IStructuralDofType DOF { get; }

		public double Multiplier { get; }

		public DomainAcceleration(IStructuralDofType dofType, double multiplier, DomainFunction domainFunction)
		{
			this.DOF = dofType;
			this.Multiplier = multiplier;
			this.DomainFunction = domainFunction;
		}

		public DomainAcceleration(IStructuralDofType dofType, double multiplier) 
			: this(dofType, multiplier, c => multiplier)
		{
		}

		public IDomainBoundaryCondition<IStructuralDofType> WithMultiplier(double multiplier) => new DomainAcceleration(DOF, multiplier, DomainFunction);
		IDomainModelQuantity<IStructuralDofType> IDomainModelQuantity<IStructuralDofType>.WithMultiplier(double multiplier) => new DomainAcceleration(DOF, multiplier, DomainFunction);
	}
}
