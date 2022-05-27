using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.Thermal.Providers
{
    public class ElementCapacityProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElementType element) =>
            element is IThermalElementType ?
                ((IThermalElementType)element).CapacityMatrix() :
                LinearAlgebra.Matrices.Matrix.CreateZero(element.GetElementDofTypes().Count, element.GetElementDofTypes().Count);
    }
}
