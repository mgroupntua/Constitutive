using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public interface IElementStructuralNeumannBoundaryCondition : IElementNeumannBoundaryCondition<IStructuralDofType>
	{
		IStructuralElementType Element { get; }
	}
}
