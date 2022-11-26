using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.Structural.Providers
{
    public class ElementStructuralMassProvider : IElementMatrixProvider
    {
		private static Matrix GetSymmetricZero(int count)
		{
			var m = LinearAlgebra.Matrices.Matrix.CreateZero(count, count);
			m.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
			return m;
		}

		public IMatrix Matrix(IElementType element) =>
            element is IStructuralElementType ?
                ((IStructuralElementType)element).MassMatrix() :
				GetSymmetricZero(element.GetElementDofTypes().Count);

		//double[] ElementalLoad(IElement eleme, IElementBoundaryCondition load) => //StiffnessMatrix * ElementLoad.Amount

		//double[] CalculateAccelerationResponseIntegral(IElement element, IList<MassAccelerationLoad> loads) { return  ElementalLoad(IElement eleme, IElementBoundaryCondition load } //TODO: go here?
	}
}
