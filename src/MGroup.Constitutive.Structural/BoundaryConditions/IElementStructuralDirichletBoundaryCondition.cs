using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public interface IElementStructuralDirichletBoundaryCondition : IElementDirichletBoundaryCondition<IStructuralDofType>
	{
		IStructuralElementType Element { get; }
	}
}
