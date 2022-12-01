using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.Structural.Providers
{
	public class ElementStructuralInternalForcesProvider : IElementVectorProvider
	{
		public double[] CalcVector(IElementType element) => element.CalculateResponseIntegral();
	}
}
