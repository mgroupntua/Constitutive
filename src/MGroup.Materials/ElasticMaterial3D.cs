using System;
using MGroup.LinearAlgebra;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Materials.Interfaces;

namespace MGroup.Materials
{
	/// <summary>
	/// It is an implementation of a linear elastic isotropic constitutive law for 3D finite elements.
	/// </summary>
	public class ElasticMaterial3D : IIsotropicContinuumMaterial3D
	{
		//private readonly double[] strains = new double[6];
		private readonly double[] stresses = new double[6];
		private Matrix constitutiveMatrix = null;
		public double YoungModulus { get; set; }
		public double PoissonRatio { get; set; }
		public double[] Coordinates { get; set; }
		private readonly double[] incrementalStrains = new double[6];
		private double[] stressesNew = new double[6];

		private Matrix GetConstitutiveMatrix()
		{
			double fE1 = YoungModulus / (double)(1 + PoissonRatio);
			double fE2 = fE1 * PoissonRatio / (double)(1 - 2 * PoissonRatio);
			double fE3 = fE1 + fE2;
			double fE4 = fE1 * 0.5;
			var afE = Matrix.CreateZero(6, 6);
			afE[0, 0] = fE3;
			afE[0, 1] = fE2;
			afE[0, 2] = fE2;
			afE[1, 0] = fE2;
			afE[1, 1] = fE3;
			afE[1, 2] = fE2;
			afE[2, 0] = fE2;
			afE[2, 1] = fE2;
			afE[2, 2] = fE3;
			afE[3, 3] = fE4;
			afE[4, 4] = fE4;
			afE[5, 5] = fE4;

			return afE;
		}

		private void CalculateNextStressStrainPoint()
		{
			var stressesElastic = new double[6];
			for (int i = 0; i < 6; i++)
			{
				stressesElastic[i] = this.stresses[i];
				for (int j = 0; j < 6; j++)
					stressesElastic[i] += this.constitutiveMatrix[i, j] * this.incrementalStrains[j];
			}

			this.stressesNew = stressesElastic;
		}

		#region IFiniteElementMaterial Members

		/// <summary>
		/// Returns the ID of the material class indicating a specific material law implementation.
		/// </summary>
		public int ID => 1;

		/// <summary>
		/// Returns a boolean indicating if the constitutive matrix of the material has changed for the current iteratively update 
		/// of the deformation state.
		/// </summary>
		public bool Modified => false;

		/// <summary>
		/// Resets the boolean that indicates if the constitutive matrix of the material has changed for the current iteratively update 
		/// of the deformation state.
		/// </summary>
		public void ResetModified() { }

		#endregion

		#region IFiniteElementMaterial3D Members
		/// <summary>
		/// Returns the stresses of this material for the current strain state.
		/// </summary>
		public double[] Stresses => stressesNew;

		/// <summary>
		/// Returns the constitutive matrix of the material for the current strain state
		/// </summary>
		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (constitutiveMatrix == null) UpdateMaterial(new double[6]);
				return constitutiveMatrix;
			}
		}

		/// <summary>
		/// Updates the material state for a given new strain point.
		/// </summary>
		/// <param name="Delta">The given strain point.</param>
		public void UpdateMaterial(double[] strainsIncrement)
		{
			//throw new NotImplementedException();
			this.incrementalStrains.CopyFrom(strainsIncrement);
			constitutiveMatrix = GetConstitutiveMatrix();
			this.CalculateNextStressStrainPoint();

		}

		/// <summary>
		/// Clears the saved stress strain point connected to the last converged analysis step.
		/// Currently not a valid operation.
		/// </summary>
		public void ClearState()
		{
			//constitutiveMatrix.Clear();
			incrementalStrains.Clear();
			stresses.Clear();
			stressesNew.Clear();
		}

		/// <summary>
		/// Saves the current stress strain state of the material (after convergence of the iterative solution process
		/// for a given loading step).
		/// </summary>
		public void SaveState() => stresses.CopyFrom(stressesNew);

		/// <summary>
		/// Clears stresses.
		/// </summary>
		public void ClearStresses()
		{
			stresses.Clear();
			stressesNew.Clear();
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
		public ElasticMaterial3D Clone()
		{
			return new ElasticMaterial3D() { YoungModulus = this.YoungModulus, PoissonRatio = this.PoissonRatio };
		}

		#endregion

	}

}
