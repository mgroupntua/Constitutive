namespace MGroup.Constitutive.Structural.Fiber
{
	/// <summary>
	/// Interface for materials laws implementations to be used in beam sections analysis 
	/// </summary>
	public interface IFiberMaterial : IStructuralMaterial
	{
		double Stress { get; }
		double Strain { get; }
		IFiberMaterial Clone(IFiberStructuralMaterial parent);
	}
}
