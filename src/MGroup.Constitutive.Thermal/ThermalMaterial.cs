using MGroup.MSolve.Constitutive;

namespace MGroup.Constitutive.Thermal
{
	public class ThermalMaterial : IThermalMaterial
	{
		/// <summary>
		/// Constructs a new material for the purposes of heat transfer applications. 
		/// This material is characterized by the following properties: specific heat, thermal conductivity and material density
		/// </summary>
		/// <param name="density">material density</param>
		/// <param name="specialHeatCoeff">specific heat,</param>
		/// <param name="thermalConductivity">thermal conductivity</param>
		public ThermalMaterial(double density, double specialHeatCoeff, double thermalConductivity)
		{
			this.Density = density;
			this.SpecialHeatCoeff = specialHeatCoeff;
			this.ThermalConductivity = thermalConductivity;
		}

		public double Density { get; }
		public double SpecialHeatCoeff { get; }
		public double ThermalConductivity { get; }

		public ThermalMaterial Clone() => new ThermalMaterial(Density, SpecialHeatCoeff, ThermalConductivity);
	}
}
