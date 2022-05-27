using MGroup.MSolve.Discretization.BoundaryConditions;

namespace MGroup.Constitutive.Thermal.BoundaryConditions
{
	public interface INodalThermalNeumannBoundaryCondition : INodalNeumannBoundaryCondition<IThermalDofType>
	{
	}
}
