using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;

//TODO: Should this be in the Problems project?
namespace MGroup.Constitutive.Structural
{
    public class ElementStructuralMassProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElement element) => element.ElementType.MassMatrix(element);
    }
}
