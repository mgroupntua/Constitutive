using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;

namespace MGroup.Constitutive.Thermal
{
    public class ElementStructuralConductivityProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElement element) => element.ElementType.StiffnessMatrix(element);
    }
}
