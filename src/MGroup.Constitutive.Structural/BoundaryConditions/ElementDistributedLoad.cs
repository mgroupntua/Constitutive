using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public class ElementDistributedLoad : IElementDistributedLoadBoundaryCondition
	{
		private DomainFunction DomainFunction;

		public double Amount(double[] coordinates) => DomainFunction(coordinates);

		public IStructuralDofType DOF { get; }

		public IStructuralElementType Element { get; }

		IElementType IElementModelQuantity<IStructuralDofType>.Element => Element;

		public double Multiplier { get; }

		public ElementDistributedLoad(IStructuralElementType element, IStructuralDofType dof, double multiplier, DomainFunction domainFunction)
		{
			Element = element;
			DOF = dof;
			Multiplier = multiplier;
			DomainFunction = domainFunction;
		}

		public ElementDistributedLoad(IStructuralElementType element, IStructuralDofType dof, double multiplier)
			: this(element, dof, multiplier, c => multiplier)
		{
		}

		public IElementBoundaryCondition<IStructuralDofType> WithMultiplier(double multiplier) => new ElementDistributedLoad(Element, DOF, multiplier, DomainFunction);

		IElementModelQuantity<IStructuralDofType> IElementModelQuantity<IStructuralDofType>.WithMultiplier(double multiplier) => WithMultiplier(multiplier);
	}
}
