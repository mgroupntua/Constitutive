using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.Structural.Providers
{
    public class ElementStructuralDampingProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElementType element) =>
            element is IStructuralElementType ?
                ((IStructuralElementType)element).DampingMatrix() :
                LinearAlgebra.Matrices.Matrix.CreateZero(element.GetElementDofTypes().Count, element.GetElementDofTypes().Count);
    }
}
