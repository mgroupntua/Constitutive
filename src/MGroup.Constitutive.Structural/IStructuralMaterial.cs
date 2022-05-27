using MGroup.MSolve.Constitutive;

namespace MGroup.Constitutive.Structural
{
	/// <summary>
	/// Interface for material law implementations to be used in structural mechanics problems
	/// </summary>
	public interface IStructuralMaterial : IConstitutiveLawWithGenericState
	{
	}
}
