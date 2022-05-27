using System.Collections.Generic;

namespace MGroup.Constitutive.Structural.Fiber
{
	/// <summary>
	/// Interface for materials laws implementations to be used in beam sections analysis 
	/// </summary>
	public interface IFiberStructuralMaterial : IStructuralMaterial
	{
		IList<IFiberMaterial> FiberMaterials { get; }
	}
}
