using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Materials.Interfaces
{
	/// <summary>
	/// Interface for materials laws implementations to be used in deformation gradient based 3D finite elements formulations
	/// </summary>
	public interface IContinuumMaterial3DDefGrad : IFiniteElementMaterial
	{
		double[] Stresses { get; }
		IMatrixView ConstitutiveMatrix { get; }
		void UpdateMaterial(double[] strains);
		void ClearState();
		void SaveState();
		void ClearStresses();
	}
}
