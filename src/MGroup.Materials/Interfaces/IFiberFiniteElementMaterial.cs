using System.Collections.Generic;

namespace MGroup.Materials.Interfaces
{
    public interface IFiberFiniteElementMaterial : IFiniteElementMaterial
    {
        IList<IFiberMaterial> FiberMaterials { get; }
    }
}
