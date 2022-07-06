namespace MGroup.Constitutive.ConvectionDiffusion
{
	public class ConvectionDiffusionProperties : IConvectionDiffusionProperties
	{
		/// <summary>
		/// Constructs a new material for the purposes of heat transfer applications. 
		/// This material is characterized by the following properties: specific heat, thermal conductivity and material density
		/// </summary>
		/// <param name="firstTimeDerivativeCoeff">coefficient of the first time derivative term.</param>
		/// <param name="diffusionCoeff">coefficient of the diffusion term.</param>
		/// <param name="convectionCoeff">coefficent of the convection term.</param>
		/// <param name="dependentSourceCoeff">coefficient of the dependent term of the linear source.</param>
		/// <param name="independentSourceCoeff">coefficient of the independent term of the linear source.</param>
		public ConvectionDiffusionProperties(double firstTimeDerivativeCoeff, double diffusionCoeff, double convectionCoeff, double dependentSourceCoeff, double independentSourceCoeff)
		{
			FirstTimeDerivativeCoeff = firstTimeDerivativeCoeff;
			DiffusionCoeff = diffusionCoeff;
			ConvectionCoeff = convectionCoeff;
			DependentSourceCoeff = dependentSourceCoeff;
			IndependentSourceCoeff = independentSourceCoeff;
		}

		public double FirstTimeDerivativeCoeff { get; }

		public double DiffusionCoeff { get; }

		public double ConvectionCoeff { get; }

		public double DependentSourceCoeff { get; }

		public double IndependentSourceCoeff { get; }

		public ConvectionDiffusionProperties Clone() => new ConvectionDiffusionProperties(FirstTimeDerivativeCoeff, DiffusionCoeff, ConvectionCoeff, DependentSourceCoeff, IndependentSourceCoeff);
	}
}
