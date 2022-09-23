namespace MGroup.Constitutive.ConvectionDiffusion
{
	public class ConvectionDiffusionProperties : IConvectionDiffusionProperties
	{
		/// <summary>
		/// Constructs a new material for the purposes of heat transfer applications. 
		/// This material is characterized by the following properties: specific heat, thermal conductivity and material density
		/// </summary>
		/// <param name="capacityCoeff">coefficient of the first time derivative term.</param>
		/// <param name="diffusionCoeff">coefficient of the diffusion term.</param>
		/// <param name="convectionCoeff">coefficent of the convection term, multiplied by the velocity vector</param>
		/// <param name="dependentSourceCoeff">coefficient of the dependent term of the linear source.</param>
		/// <param name="independentSourceCoeff">coefficient of the independent term of the linear source.</param>
		public ConvectionDiffusionProperties(double capacityCoeff, double diffusionCoeff, double[] convectionCoeff, double dependentSourceCoeff, double independentSourceCoeff)
		{
			CapacityCoeff = capacityCoeff;
			DiffusionCoeff = diffusionCoeff;
			ConvectionCoeff = convectionCoeff;
			DependentSourceCoeff = dependentSourceCoeff;
			IndependentSourceCoeff = independentSourceCoeff;
		}

		public double CapacityCoeff { get; }

		public double DiffusionCoeff { get; }

		public double[] ConvectionCoeff { get; }

		public double DependentSourceCoeff { get; }

		public double IndependentSourceCoeff { get; }

		public ConvectionDiffusionProperties Clone() => new ConvectionDiffusionProperties(CapacityCoeff, DiffusionCoeff, ConvectionCoeff, DependentSourceCoeff, IndependentSourceCoeff);
	}
}
