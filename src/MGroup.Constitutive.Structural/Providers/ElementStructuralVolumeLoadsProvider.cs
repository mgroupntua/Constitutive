using MGroup.Constitutive.Structural;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.Structural.Providers
{
	public class ElementStructuralVolumeLoadsProvider : IElementVectorProvider
	{
		public double[] CalcVector(IElementType element) =>
			element is IStructuralElementType ?
				((IStructuralElementType)element).VolumeLoads() :
				new double[element.GetElementDofTypes().Count];
		//public IMatrix Matrix(IElementType element) => throw new System.NotImplementedException();
	}
}
