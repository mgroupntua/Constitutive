using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.Structural.Providers
{
	public class ElementConvectionDiffusionInternalForcesProvider : IElementVectorProvider
	{
		public double[] CalcVector(IElementType element) => element.CalculateResponseIntegral();
	}
}
