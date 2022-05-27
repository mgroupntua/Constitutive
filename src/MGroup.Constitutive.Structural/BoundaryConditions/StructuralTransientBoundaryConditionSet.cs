using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.Constitutive.Structural.BoundaryConditions
{
	public class StructuralTransientBoundaryConditionSet : TransientBoundaryConditionSet<IStructuralDofType>
	{
		public StructuralTransientBoundaryConditionSet(IEnumerable<IBoundaryConditionSet<IStructuralDofType>> boundaryConditionSets, Func<double, double, double> timeFunc) 
			: base(boundaryConditionSets, timeFunc)
		{
		}

		public override IBoundaryConditionSet<IStructuralDofType> CreateBoundaryConditionSetOfSubdomain(ISubdomain subdomain) =>
			new StructuralTransientBoundaryConditionSet(boundaryConditionSets.Select(x => x.CreateBoundaryConditionSetOfSubdomain(subdomain)), timeFunc);

		public override IEnumerable<INodalNeumannBoundaryCondition<IStructuralDofType>> EnumerateEquivalentNodalNeumannBoundaryConditions(IEnumerable<IElementType> elements) => 
			boundaryConditionSets.SelectMany(x => x.EnumerateEquivalentNodalNeumannBoundaryConditions(elements));
	}
}
