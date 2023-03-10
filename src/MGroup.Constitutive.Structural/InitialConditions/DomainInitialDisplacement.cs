using MGroup.MSolve.DataStructures;
using MGroup.MSolve.AnalysisWorkflow.Transient;
using MGroup.MSolve.Discretization;

namespace MGroup.Constitutive.Structural.InitialConditions
{
	public class DomainInitialDisplacement : IDomainDisplacementInitialCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IStructuralDofType DOF { get; }

		public double Multiplier { get; }

		public DifferentiationOrder DifferentiationOrder => DifferentiationOrder.Zero;

		public DomainInitialDisplacement(IStructuralDofType dofType, double multiplier, DomainFunction domainFunction)
		{
			this.DOF = dofType;
			this.Multiplier = multiplier;
			this.DomainFunction = domainFunction;
		}

		public DomainInitialDisplacement(IStructuralDofType dofType, double multiplier)
			: this(dofType, multiplier, c => multiplier)
		{
		}

		IDomainInitialCondition<IStructuralDofType> IDomainInitialCondition<IStructuralDofType>.WithMultiplier(double multiplier) => new DomainInitialDisplacement(DOF, multiplier, DomainFunction);
		IDomainModelQuantity<IStructuralDofType> IDomainModelQuantity<IStructuralDofType>.WithMultiplier(double multiplier) => new DomainInitialDisplacement(DOF, multiplier, DomainFunction);
	}
}
