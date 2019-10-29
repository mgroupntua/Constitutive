using System;
using MGroup.MSolve.Constitutive;

namespace MGroup.Constitutive.Structural
{
	/// <summary>
	/// Elastic material properties to be assigned on beam finite elements.
	/// </summary>
	public class ElasticMaterial : IFiniteElementMaterial
	{
		private readonly double[] strains = new double[3];
		private readonly double[] incrementalStrains = new double[3];
		private readonly double[] stresses = new double[3];
		private double[] stressesNew = new double[3];
		private double[,] constitutiveMatrix = null;
		public double YoungModulus { get; set; }
		public double PoissonRatio { get; set; }
		public double[] Coordinates { get; set; }

		#region IFiniteElementMaterial Members

		/// <summary>
		/// Returns the ID of the material class indicating a specific material law implementation.
		/// </summary>
		public int ID
		{
			get { return 1; }
		}

		/// <summary>
		/// Returns a boolean indicating if the constitutive matrix of the material has changed for the current iteratively update 
		/// of the deformation state.
		/// </summary>
		public bool Modified
		{
			get { return false; }
		}

		/// <summary>
		/// Resets the boolean that indicates if the constitutive matrix of the material has changed for the current iteratively update 
		/// of the deformation state.
		/// </summary>
		public void ResetModified()
		{
		}

		#endregion

		#region ICloneable Members
		/// <summary>
		/// Creates a clone of material object with the same parameters.
		/// </summary>
		/// <returns>The created material clone</returns>
		object ICloneable.Clone() => Clone();

		/// <summary>
		/// Creates a clone of material object with the same parameters.
		/// </summary>
		/// <returns>The created material clone</returns>
		public ElasticMaterial Clone()
		{
			return new ElasticMaterial() { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
		}

		/// <summary>
		/// Saves the current stress strain state of the material (after convergence of the iterative solution process
		/// for a given loading step).
		/// </summary>
		public void SaveState()
		{
			Array.Copy(this.stressesNew, this.stresses, 3);
		}

		/// <summary>
		/// Clears the saved stress strain point connected to the last converged analysis step.
		/// </summary>
		public void ClearState()
		{
			if (constitutiveMatrix != null) Array.Clear(constitutiveMatrix, 0, constitutiveMatrix.Length);
			Array.Clear(incrementalStrains, 0, incrementalStrains.Length);
			Array.Clear(stresses, 0, stresses.Length);
			Array.Clear(stressesNew, 0, stressesNew.Length);
		}

		/// <summary>
		/// Clears stresses - Currently not a valid operation.
		/// </summary>
		public void ClearStresses()
		{
			throw new NotImplementedException();
		}

		#endregion
	}
}
