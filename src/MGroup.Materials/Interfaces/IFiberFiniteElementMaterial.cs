using System.Collections.Generic;

namespace MGroup.Materials.Interfaces
{
	/// <summary>
	/// Interface for materials laws implementations to be used in beam sections analysis 
	/// </summary>
	public interface IFiberFiniteElementMaterial : IFiniteElementMaterial
	{
		IList<IFiberMaterial> FiberMaterials { get; }
	}
}
