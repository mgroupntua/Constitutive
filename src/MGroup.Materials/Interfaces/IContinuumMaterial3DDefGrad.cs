using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Materials.Interfaces
{
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
