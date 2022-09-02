using MGroup.MSolve.Discretization;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Constitutive.ConvectionDiffusion
{
	public interface IConvectionDiffusionElementType : IElementType
	{
		IMatrix DiffusionMatrix();

		IMatrix ConvectionMatrix();

		IMatrix CapacityMatrix();

		IMatrix ProductionMatrix();

		double[] ProductionVector();
	}
}
