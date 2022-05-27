using System;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;

namespace MGroup.Constitutive.Structural.Planar
{
	/// <summary>
	/// It is an implementation of a linear elastic isotropic constitutive law for 2D finite elements.
	/// </summary>
	public class ElasticMaterial2D : IIsotropicContinuumMaterial2D
	{
		private const string STRESS_X = "Stress X";
		private const string STRESS_Y = "Stress Y";
		private const string STRESS_XY = "Stress XY";
		
		private readonly double[] strains = new double[3];
		private double[] stresses = new double[3];
		private double[] stressesNew = new double[3];
		private Matrix constitutiveMatrix = null;
		private GenericConstitutiveLawState currentState;

		public double PoissonRatio { get; }
		public StressState2D StressState { get; }
		public double YoungModulus { get; }

		/// <summary>
		/// Creastes a new object of <see cref="ElasticMaterial2D"/> class.
		/// </summary>
		/// <param name="stressState"> The stress strain state of the 2D problem.</param>
		public ElasticMaterial2D(double youngModulus, double poissonRatio, StressState2D stressState)
		{
			this.StressState = stressState;
			this.YoungModulus = youngModulus;
			this.PoissonRatio = poissonRatio;
			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(STRESS_X, 0d),
				(STRESS_Y, 0d),
				(STRESS_XY, 0d),
			});
		}

		#region IFiniteElementMaterial3D

		/// <summary>
		/// Returns the constitutive matrix of the material for the current strain state
		/// </summary>
		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (constitutiveMatrix == null) UpdateConstitutiveMatrixAndEvaluateResponse(new double[3]);
				return constitutiveMatrix;
			}
		}

		/// <summary>
		/// Returns the stresses of this material for the current strain state.
		/// </summary>
		public double[] Stresses => stresses;

		/// <summary>
		/// Clears the saved stress strain point connected to the last converged analysis step.
		/// Currently not a valid operation.
		/// </summary>
		public void ClearState()
		{
			strains.Clear();
			stresses.Clear();
		}

		/// <summary>
		/// Clears stresses - Currently not a valid operation.
		/// </summary>
		public void ClearStresses()
		{
			throw new NotImplementedException();
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
				(STRESS_XY, stresses[2]),
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
				stresses[2] = currentState.StateValues[STRESS_XY];
			}
		}

		/// <summary>
		/// Updates the material state for a given new strain point.
		/// </summary>
		/// <param name="Delta">The given strain point.</param>
		public double[] UpdateConstitutiveMatrixAndEvaluateResponse(double[] strains)
		{
			this.strains.CopyFrom(strains);
			constitutiveMatrix = Matrix.CreateZero(3, 3); //TODO: This should be cached in the constitutive matrix property and used here.
			if (StressState == StressState2D.PlaneStress)
			{
				double aux = YoungModulus / (1 - PoissonRatio * PoissonRatio);
				constitutiveMatrix[0, 0] = aux;
				constitutiveMatrix[1, 1] = aux;
				constitutiveMatrix[0, 1] = PoissonRatio * aux;
				constitutiveMatrix[1, 0] = PoissonRatio * aux;
				constitutiveMatrix[2, 2] = (1 - PoissonRatio) / 2 * aux;
			}
			else
			{
				double aux = YoungModulus / (1 + PoissonRatio) / (1 - 2 * PoissonRatio);
				constitutiveMatrix[0, 0] = aux * (1 - PoissonRatio);
				constitutiveMatrix[1, 1] = aux * (1 - PoissonRatio);
				constitutiveMatrix[0, 1] = PoissonRatio * aux;
				constitutiveMatrix[1, 0] = PoissonRatio * aux;
				constitutiveMatrix[2, 2] = (1 - 2 * PoissonRatio) / 2 * aux;
			}

			var stressesElastic = new double[6];
			for (int i = 0; i < 3; i++)
			{
				stressesElastic[i] = this.stresses[i];
				for (int j = 0; j < 3; j++)
					stressesElastic[i] += this.constitutiveMatrix[i, j] * this.strains[j];
			}

			this.stressesNew = stressesElastic;
			return stressesNew;
			//// correction for stress calculation 
			//stresses = (constitutiveMatrix * Vector.CreateFromArray( strains)).CopyToArray();
			//// TODO correct strains add update strains depends of how we define update strain tou ulikou

			//return stresses;
		}

		#endregion

		#region IFiniteElementMaterial

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
		public ElasticMaterial2D Clone() => new ElasticMaterial2D(YoungModulus, PoissonRatio, StressState);

		#endregion
		/// <summary>
		/// Returns coordinates. Currently not a valid operation.
		/// </summary>
		public double[] Coordinates { get; set; }
	}
}
