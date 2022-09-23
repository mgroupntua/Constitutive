namespace MGroup.Constitutive.ConvectionDiffusion
{
	/// <summary>
	/// Contains material properties used for  problems. These are uniform for the whole element and immutable.
	/// Authors: Orestis Papas, Christodoulou Theofilos.
	/// </summary>
	public interface IConvectionDiffusionProperties
	{
		double CapacityCoeff { get; }

		double DiffusionCoeff { get; }

		double[] ConvectionCoeff { get; }

		double DependentSourceCoeff { get; }

		double IndependentSourceCoeff { get; }
	}
}
