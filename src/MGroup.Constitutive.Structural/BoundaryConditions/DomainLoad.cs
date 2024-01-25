using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public class DomainLoad : IDomainLoadBoundaryCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IStructuralDofType DOF { get; }

		public double Multiplier { get; }

		public DomainLoad(IStructuralDofType dof, double multiplier, DomainFunction domainFunction)
		{
			DOF = dof;
			Multiplier = multiplier;
			DomainFunction = domainFunction;
		}

		public DomainLoad(IStructuralDofType dof, double multiplier)
			: this(dof, multiplier, c => multiplier)
		{
		}

		public IDomainBoundaryCondition<IStructuralDofType> WithMultiplier(double multiplier) => new DomainLoad(DOF, multiplier, DomainFunction);
		IDomainModelQuantity<IStructuralDofType> IDomainModelQuantity<IStructuralDofType>.WithMultiplier(double multiplier) => new DomainLoad(DOF, multiplier, DomainFunction);
	}
}
