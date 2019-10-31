using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;

namespace MGroup.Constitutive.Thermal
{
    public class ElementStructuralCapacityProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElement element) => element.ElementType.MassMatrix(element);
    }
}
