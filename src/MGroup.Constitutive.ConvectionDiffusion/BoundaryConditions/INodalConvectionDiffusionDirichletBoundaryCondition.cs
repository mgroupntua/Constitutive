using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions
{
	public interface INodalConvectionDiffusionDirichletBoundaryCondition : INodalDirichletBoundaryCondition<IConvectionDiffusionDofType>
	{
	}
}
