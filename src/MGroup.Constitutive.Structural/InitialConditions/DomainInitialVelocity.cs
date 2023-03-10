using MGroup.MSolve.DataStructures;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization;

namespace MGroup.Constitutive.Structural.InitialConditions
{
	public class DomainInitialVelocity : IDomainVelocityInitialCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IStructuralDofType DOF { get; }

		public double Multiplier { get; }

		public DifferentiationOrder DifferentiationOrder => DifferentiationOrder.First;

		public DomainInitialVelocity(IStructuralDofType dofType, double multiplier, DomainFunction domainFunction)
		{
			this.DOF = dofType;
			this.Multiplier = multiplier;
			this.DomainFunction = domainFunction;
		}

		public DomainInitialVelocity(IStructuralDofType dofType, double multiplier)
			: this(dofType, multiplier, c => multiplier)
		{
		}

		IDomainModelQuantity<IStructuralDofType> IDomainModelQuantity<IStructuralDofType>.WithMultiplier(double multiplier) => new DomainInitialVelocity(DOF, multiplier);
		IDomainInitialCondition<IStructuralDofType> IDomainInitialCondition<IStructuralDofType>.WithMultiplier(double multiplier) => new DomainInitialVelocity(DOF, multiplier);
	}
}
