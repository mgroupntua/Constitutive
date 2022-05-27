using System;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;

namespace MGroup.Constitutive.Structural.Continuum
{
	/// <summary>
	/// It is an implementation of a linear elastic isotropic constitutive law for 3D finite elements.
	/// </summary>
	public class ElasticMaterial3D : IIsotropicContinuumMaterial3D
	{
		private const string STRESS_X = "Stress X";
		private const string STRESS_Y = "Stress Y";
		private const string STRESS_Z = "Stress Z";
		private const string STRESS_XY = "Stress XY";
		private const string STRESS_XZ = "Stress XZ";
		private const string STRESS_YZ = "Stress YZ";

		//private readonly double[] strains = new double[6];
		private readonly double[] stresses = new double[6];
		private Matrix constitutiveMatrix = null;
		public double YoungModulus { get; }
		public double PoissonRatio { get; }
		//public double[] Coordinates { get; set; }
		private readonly double[] incrementalStrains = new double[6];
		private double[] stressesNew = new double[6];
		private GenericConstitutiveLawState currentState;

		public ElasticMaterial3D(double youngModulus, double poissonRatio)
		{
			YoungModulus = youngModulus;
			PoissonRatio = poissonRatio;
			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(STRESS_X, 0d),
				(STRESS_Y, 0d),
				(STRESS_Z, 0d),
				(STRESS_XY, 0d),
				(STRESS_XZ, 0d),
				(STRESS_YZ, 0d),
			});
		}

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
		public bool IsCurrentStateDifferent() => false;

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
				if (constitutiveMatrix == null) UpdateConstitutiveMatrixAndEvaluateResponse(new double[6]);
				return constitutiveMatrix;
			}
		}

		/// <summary>
		/// Updates the material state for a given new strain point.
		/// </summary>
		/// <param name="Delta">The given strain point.</param>
		public double[] UpdateConstitutiveMatrixAndEvaluateResponse(double[] strainsIncrement)
		{
			//throw new NotImplementedException();
			this.incrementalStrains.CopyFrom(strainsIncrement);
			constitutiveMatrix = GetConstitutiveMatrix();
			this.CalculateNextStressStrainPoint();

			return stressesNew;
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
		public GenericConstitutiveLawState CreateState()
		{
			stresses.CopyFrom(stressesNew);
			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(STRESS_X, stresses[0]),
				(STRESS_Y, stresses[1]),
				(STRESS_Z, stresses[2]),
				(STRESS_XY, stresses[3]),
				(STRESS_XZ, stresses[4]),
				(STRESS_YZ, stresses[5]),
			});

			return currentState;
		}
		IHaveState ICreateState.CreateState() => CreateState();
		public GenericConstitutiveLawState CurrentState
		{
			get => currentState;
			set
			{
				currentState = value;
				stresses[0] = currentState.StateValues[STRESS_X];
				stresses[1] = currentState.StateValues[STRESS_Y];
				stresses[2] = currentState.StateValues[STRESS_Z];
				stresses[3] = currentState.StateValues[STRESS_XY];
				stresses[4] = currentState.StateValues[STRESS_XZ];
				stresses[5] = currentState.StateValues[STRESS_YZ];
			}
		}

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
		public ElasticMaterial3D Clone() => new ElasticMaterial3D(YoungModulus, PoissonRatio);

		#endregion

	}

}
