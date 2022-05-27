using MGroup.LinearAlgebra.Matrices;
using MGroup.Constitutive.Structural;
using MGroup.MSolve.Discretization;

namespace MGroup.Constitutive.PorousMedia
{
	public interface IPorousElementType : IStructuralElementType
	{
		IMatrix PermeabilityMatrix();
		IMatrix CouplingMatrix();
		IMatrix SaturationMatrix();
	}
}
