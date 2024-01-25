using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.ConvectionDiffusion.Providers
{
	public class ElementConvectionDiffusionInternalForcesProvider : IElementVectorProvider
	{
		public double[] CalcVector(IElementType element) => element.CalculateResponseIntegral();
	}
}
