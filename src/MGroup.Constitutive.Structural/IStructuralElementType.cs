using MGroup.MSolve.Discretization;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Constitutive.Structural
{
	public interface IStructuralElementType : IElementType
	{
		IMatrix StiffnessMatrix();
		IMatrix MassMatrix();
		IMatrix DampingMatrix();
	}
}
