using MGroup.MSolve.Constitutive;

namespace MGroup.Constitutive.Structural.Fiber
{ 
	public interface IFiber
	{
		IConstitutiveLaw Material { get; }
	}
}
