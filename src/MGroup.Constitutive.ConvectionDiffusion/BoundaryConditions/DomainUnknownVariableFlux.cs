using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions
{
	public class DomainUnknownVariableFlux : IDomainUnknownVariableFluxBoundaryCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IConvectionDiffusionDofType DOF { get; }

		public double Multiplier { get; }

		public DomainUnknownVariableFlux(IConvectionDiffusionDofType dof, double multiplier, DomainFunction domainFunction)
		{
			DOF = dof;
			Multiplier = multiplier;
			DomainFunction = domainFunction;
		}

		public DomainUnknownVariableFlux(IConvectionDiffusionDofType dof, double multiplier)
			: this(dof, multiplier, c => multiplier)
		{
		}

		public IDomainBoundaryCondition<IConvectionDiffusionDofType> WithMultiplier(double multiplier) => new DomainUnknownVariableFlux(DOF, multiplier, DomainFunction);
		IDomainModelQuantity<IConvectionDiffusionDofType> IDomainModelQuantity<IConvectionDiffusionDofType>.WithMultiplier(double multiplier) => new DomainUnknownVariableFlux(DOF, multiplier, DomainFunction);
	}
}
