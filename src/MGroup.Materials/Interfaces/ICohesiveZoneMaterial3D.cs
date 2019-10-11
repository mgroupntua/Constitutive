using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Materials.Interfaces
{
	/// <summary>
	/// Interface that traction separation laws implementations should follow
	/// </summary>
	public interface ICohesiveZoneMaterial3D  
	{
		double[] Tractions { get; }
		IMatrixView ConstitutiveMatrix { get; }
		void UpdateMaterial(double[] strains);
		int ID { get; }
		bool Modified { get; }
		void ResetModified();
		void SaveState();
		void ClearState();
		void ClearTractions();
		ICohesiveZoneMaterial3D Clone();
	}
}
