using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Constitutive.Structural.Shells
{
	/// <summary>
	/// Interface for materials laws implementations to be used with shell finite elements
	/// </summary>
	public interface IShellMaterial : IStructuralMaterial
	{
		double[] NormalVectorV3 { set; }
		double[] TangentVectorV2 { set; }
		double[] TangentVectorV1 { set; }
	}
}
