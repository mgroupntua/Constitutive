namespace MGroup.Constitutive.Thermal
{
	public class ThermalProperties : IThermalProperties
	{
		/// <summary>
		/// Constructs a new material for the purposes of heat transfer applications. 
		/// This material is characterized by the following properties: specific heat, thermal conductivity and material density
		/// </summary>
		/// <param name="density">material density</param>
		/// <param name="specialHeatCoeff">specific heat,</param>
		/// <param name="thermalConductivity">thermal conductivity</param>
		public ThermalProperties(double density, double specialHeatCoeff, double thermalConductivity)
		{
			this.Density = density;
			this.SpecialHeatCoeff = specialHeatCoeff;
			this.ThermalConductivity = thermalConductivity;
		}

		public double Density { get; }
		public double SpecialHeatCoeff { get; }
		public double ThermalConductivity { get; }

		public ThermalProperties Clone() => new ThermalProperties(Density, SpecialHeatCoeff, ThermalConductivity);
	}
}
