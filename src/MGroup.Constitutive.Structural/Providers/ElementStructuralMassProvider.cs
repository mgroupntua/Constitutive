using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.Structural.Providers
{
    public class ElementStructuralMassProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElementType element) =>
            element is IStructuralElementType ?
                ((IStructuralElementType)element).MassMatrix() :
                LinearAlgebra.Matrices.Matrix.CreateZero(element.GetElementDofTypes().Count, element.GetElementDofTypes().Count);

        //double[] ElementalLoad(IElement eleme, IElementBoundaryCondition load) => //StiffnessMatrix * ElementLoad.Amount

        //double[] CalculateAccelerationResponseIntegral(IElement element, IList<MassAccelerationLoad> loads) { return  ElementalLoad(IElement eleme, IElementBoundaryCondition load } //TODO: go here?
    }
}
