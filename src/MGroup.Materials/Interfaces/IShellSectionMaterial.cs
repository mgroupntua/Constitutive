using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Materials.Interfaces
{
    public interface IShellSectionMaterial:IFiniteElementMaterial
	{
		new IShellSectionMaterial Clone();
		double[] MembraneForces { get; }
		double[] Moments { get; }
		IMatrixView MembraneConstitutiveMatrix { get; }
		IMatrixView BendingConstitutiveMatrix { get; }
		IMatrixView CouplingConstitutiveMatrix { get; }
		void UpdateMaterial(double[] membraneStrains, double[] bendingStrains);
	}
}