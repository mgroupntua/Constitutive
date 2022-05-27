using System.Collections.Generic;
using MGroup.MSolve.Constitutive;

namespace MGroup.Constitutive.Structural.Fiber
{
	public interface IStructuralElementWithFibersType : IStructuralElementType
	{
		IConstitutiveLaw Material { get; }
		IList<IFiber> Fibers { get; }
	}
}
