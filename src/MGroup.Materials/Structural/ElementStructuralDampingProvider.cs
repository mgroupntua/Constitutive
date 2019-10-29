using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Interfaces;

//TODO: Should this be in the Problems project?
namespace MGroup.Constitutive.Structural
{
    public class ElementStructuralDampingProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElement element) => element.ElementType.DampingMatrix(element);
    }
}