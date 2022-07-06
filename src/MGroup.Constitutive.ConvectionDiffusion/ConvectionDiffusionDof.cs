using MGroup.MSolve.Discretization.Dofs;

namespace MGroup.Constitutive.ConvectionDiffusion
{
	/// <summary>
	/// Degree of freedom corresponding to the uknown coefficient at a single point. Implements enum pattern.
	/// Authors: Orestis Papas, Christodoulou Theofilos.
	/// </summary>
	public class ConvectionDiffusionDof : IConvectionDiffusionDofType
	{
		/// <summary>
		/// This dof corresponds to the unknown coefficient.
		/// </summary>
		public static readonly ConvectionDiffusionDof UnknownVariable = new ConvectionDiffusionDof("UnknownVariable");

		private readonly string name;

		private ConvectionDiffusionDof(string name) => this.name = name;

		public override string ToString() => name;
	}
}
