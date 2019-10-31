using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Loads;

//TODO: time logging must be refactored
//TODO: perhaps this belongs to Solvers.Assemblers, since the vector type depends on the solver. In that case, the 
//      elementMatrixProvider should be injected by the problem/provider.
namespace MGroup.Constitutive.Thermal
{
	/// <summary>
	/// Calculates the equivalent nodal forces (at the subdomain level) due to Dirichlet boundary conditions.
	/// Authors: Maria Tavlaki
	/// </summary>
	public class DirichletEquivalentLoadsThermal : IDirichletEquivalentLoadsAssembler
	{
		private IElementMatrixProvider elementProvider; //TODO: not sure if df = K * du is the best way to calcuate df.

		public DirichletEquivalentLoadsThermal(IElementMatrixProvider elementProvider)
		{
			this.elementProvider = elementProvider;
		}

		public IVector GetEquivalentNodalLoads(ISubdomain subdomain, IVectorView solution, double constraintScalingFactor) 
		{
			//var times = new Dictionary<string, TimeSpan>();
			//var totalStart = DateTime.Now;
			//times.Add("rowIndexCalculation", DateTime.Now - totalStart);
			//times.Add("element", TimeSpan.Zero);
			//times.Add("addition", TimeSpan.Zero);

			var subdomainEquivalentForces = Vector.CreateZero(subdomain.FreeDofOrdering.NumFreeDofs);
			foreach (IElement element in subdomain.Elements) //TODO: why go through all the elements? Most of them will not have Dirichlet bc.
			{
				//var elStart = DateTime.Now;
				IMatrix elementK = elementProvider.Matrix(element);

				//double[] localSolution = subdomain.CalculateElementNodalDisplacements(element, solution);
				//double[] localdSolution = subdomain.CalculateElementIcrementalConstraintDisplacements(element, constraintScalingFactor);
				double[] localdSolution = 
					subdomain.CalculateElementIncrementalConstraintDisplacements(element, constraintScalingFactor);

				var elementEquivalentForces = elementK.Multiply(localdSolution);

				subdomain.FreeDofOrdering.AddVectorElementToSubdomain(element, elementEquivalentForces, subdomainEquivalentForces);

				//times["addition"] += DateTime.Now - elStart;
			}

			//var totalTime = DateTime.Now - totalStart;

			return subdomainEquivalentForces;
		}        
	}
}
