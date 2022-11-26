using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.ConvectionDiffusion.Providers
{
	public class ElementProductionProvider : IElementMatrixProvider
	{
		private static Matrix GetSymmetricZero(int count)
		{
			var m = LinearAlgebra.Matrices.Matrix.CreateZero(count, count);
			m.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
			return m;
		}

		public IMatrix Matrix(IElementType element) =>
			element is IConvectionDiffusionElementType ?
				((IConvectionDiffusionElementType)element).ProductionMatrix() :
				GetSymmetricZero(element.GetElementDofTypes().Count);
	}
}
