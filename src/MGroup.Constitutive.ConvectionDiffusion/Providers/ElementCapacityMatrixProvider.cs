using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.ConvectionDiffusion.Providers
{
	public class ElementCapacityMatrixProvider : IElementMatrixProvider
	{
		public IMatrix Matrix(IElementType element) =>
			element is IConvectionDiffusionElementType ?
				((IConvectionDiffusionElementType)element).CapacityMatrix() :
				LinearAlgebra.Matrices.Matrix.CreateZero(element.GetElementDofTypes().Count, element.GetElementDofTypes().Count);
	}
}
