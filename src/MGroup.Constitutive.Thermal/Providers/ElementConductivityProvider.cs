using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.Thermal.Providers
{
    public class ElementConductivityProvider : IElementMatrixProvider
    {
		private static Matrix GetSymmetricZero(int count)
		{
			var m = LinearAlgebra.Matrices.Matrix.CreateZero(count, count);
			m.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
			return m;
		}

		public IMatrix Matrix(IElementType element) =>
            element is IThermalElementType ?
                ((IThermalElementType)element).ConductivityMatrix() :
				GetSymmetricZero(element.GetElementDofTypes().Count);
	}
}
