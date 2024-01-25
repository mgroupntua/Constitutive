using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions
{
	public class DomainUnknownVariable : IDomainUnknownVariableBoundaryCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IConvectionDiffusionDofType DOF { get; }

		public double Multiplier { get; }

		public DomainUnknownVariable(IConvectionDiffusionDofType dofType, double multiplier, DomainFunction domainFunction)
		{
			this.DOF = dofType;
			this.Multiplier = multiplier;
			this.DomainFunction = domainFunction;
		}

		public DomainUnknownVariable(IConvectionDiffusionDofType dofType, double multiplier)
			: this(dofType, multiplier, c => multiplier)
		{
		}

		public IDomainBoundaryCondition<IConvectionDiffusionDofType> WithMultiplier(double multiplier) => new DomainUnknownVariable(DOF, multiplier, DomainFunction);
		IDomainModelQuantity<IConvectionDiffusionDofType> IDomainModelQuantity<IConvectionDiffusionDofType>.WithMultiplier(double multiplier) => new DomainUnknownVariable(DOF, multiplier, DomainFunction);
	}
}
