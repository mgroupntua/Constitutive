namespace MGroup.Materials.Interfaces
{
	/// <summary>
	/// Interface for materials laws implementations to be used in beam sections analysis 
	/// </summary>
	public interface IFiberMaterial : IFiniteElementMaterial
	{
		double Stress { get; }
		double Strain { get; }
		void UpdateMaterial(double dStrain);
		void SaveState();
		void ClearStresses();
		IFiberMaterial Clone(IFiberFiniteElementMaterial parent);
		double YoungModulus { get; set; }
		double PoissonRatio { get; set; } //It might be useless
		//double YoungModulusElastic { get; set; }
	}
}
