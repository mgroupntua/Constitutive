//TODO: Density is also used in body forces, so it is not only for dynamic/modal problems.
//TODO: These cannot vary throughout the element. The solution is having a MaterialField that returns properties at specific 
//      natural points.
using MGroup.MSolve.Constitutive;

namespace MGroup.Constitutive.Structural
{
	/// <summary>
	/// Contains material properties used for dynamic or modal analysis. These are uniform for the whole element and immutable.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public class DynamicMaterial : IDynamicMaterial
	{
		/// <summary>
		/// Creates an object of class <see cref="DynamicMaterial"/> containing dynamic material properties.
		/// </summary>
		/// <param name="density"> Material density</param>
		/// <param name="rayleighCoeffMass"> </param>
		/// <param name="rayleighCoeffStiffness"></param>
		public DynamicMaterial(double density, double rayleighCoeffMass, double rayleighCoeffStiffness)
		{
			this.Density = density;
			this.RayleighCoeffMass = rayleighCoeffMass;
			this.RayleighCoeffStiffness = rayleighCoeffStiffness;
		}

		public double Density { get; }
		public double RayleighCoeffMass { get; }
		public double RayleighCoeffStiffness { get; }
	}
}
