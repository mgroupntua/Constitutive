using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Constitutive.Structural.Shells
{
	/// <summary>
	/// Interface for materials laws implementations to be used with shell finite elements with through section thickness integration 
	/// embedded in the material law code
	/// </summary>
	public interface IShellSectionMaterial : IStructuralMaterial
	{
		//new IShellSectionMaterial Clone();
		double[] MembraneForces { get; }
		double[] Moments { get; }
		IReadOnlyMatrix MembraneConstitutiveMatrix { get; }
		IReadOnlyMatrix BendingConstitutiveMatrix { get; }
		IReadOnlyMatrix CouplingConstitutiveMatrix { get; }
		void UpdateMaterial(double[] membraneStrains, double[] bendingStrains);
	}
}
