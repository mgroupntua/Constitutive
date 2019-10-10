using System.Collections.Generic;
using MGroup.FEM.Entities;
using MGroup.MSolve.Discretization.FreedomDegrees;
using MGroup.MSolve.Discretization.Interfaces;

namespace MGroup.Multiscale.Interfaces
{
	/// <summary>
	/// Indicates the primary stress and strain tensor that will be used in the used multiscale homogenization scheme for the micro to macro transitions
	/// Authors: Gerasimos Sotiropoulos
	/// </summary>
	public interface IScaleTransitions
	{
		double[] MacroToMicroTransition(Node boundaryNode, double[] MacroScaleVariable);
		double[] MicroToMacroTransition(INode boundaryNode, double[] MicroScaleVariable);
		void ModifyMicrostructureTotalPrescribedBoundaryDisplacementsVectorForMacroStrainVariable(INode boundaryNode,
			double[] MacroScaleVariable, Dictionary<int, Dictionary<IDofType, double>> totalPrescribedBoundaryDisplacements);
		void ImposeAppropriateConstraintsPerBoundaryNode(Model model, Node boundaryNode);
		void ImposeAppropriateAndRigidBodyConstraintsPerBoundaryNode(Model model, Node boundaryNode, Dictionary<Node, IList<IDofType>> RigidBodyNodeConstraints); //TODO: enopoihsh twn duo duplicate
		int PrescribedDofsPerNode(); // TODO: pithanws epistrofh kai to poioi einai me input sugkekrimeno node
		int MacroscaleVariableDimension();
	}
}
