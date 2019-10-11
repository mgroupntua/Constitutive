using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Materials.Interfaces
{
	public interface IContinuumMaterial3D : IFiniteElementMaterial
	{
		/// <summary>
		/// Interface for materials laws implementations to be used in 3D finite elements
		/// </summary>
		double[] Stresses { get; }
		IMatrixView ConstitutiveMatrix { get; }
		void UpdateMaterial(double[] strains);
	}
}
