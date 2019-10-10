using System;
using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization.Interfaces;

namespace MGroup.Multiscale.Interfaces
{
	/// <summary>
	/// Indicates the nesessary methods that should be implemented by builders of models intended to be used as rves in Multiscale Problems
	/// Authors: Gerasimos Sotiropoulos
	/// </summary>
	public interface IRVEbuilder
	{
		Tuple<Model, Dictionary<int, INode>,double> GetModelAndBoundaryNodes();
		IRVEbuilder Clone(int a);
	}
}
