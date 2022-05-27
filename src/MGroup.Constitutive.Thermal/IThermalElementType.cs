using MGroup.MSolve.Discretization;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Constitutive.Thermal
{
	public interface IThermalElementType : IElementType
	{
		IMatrix ConductivityMatrix();

		IMatrix CapacityMatrix();
	}
}
