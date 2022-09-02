using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.ConvectionDiffusion.Providers
{
	public class ElementIndependentProductionProvider : IElementVectorProvider
	{
		public double[] CalcVector(IElementType element) =>
			element is IConvectionDiffusionElementType ?
				((IConvectionDiffusionElementType)element).ProductionVector() :
				new double[element.GetElementDofTypes().Count];
		public IMatrix Matrix(IElementType element) => throw new System.NotImplementedException();
	}
}
