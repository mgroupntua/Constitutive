using System.Xml.Linq;

using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.Structural.Providers
{
    public class ElementStructuralStiffnessProvider : IElementMatrixProvider
    {
		private static Matrix GetSymmetricZero(int count)
		{
			var m = LinearAlgebra.Matrices.Matrix.CreateZero(count, count);
			m.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
			return m;
		}

		public IMatrix Matrix(IElementType element) =>
            element is IStructuralElementType ?
                ((IStructuralElementType)element).StiffnessMatrix() :
                GetSymmetricZero(element.GetElementDofTypes().Count);
    }
}

