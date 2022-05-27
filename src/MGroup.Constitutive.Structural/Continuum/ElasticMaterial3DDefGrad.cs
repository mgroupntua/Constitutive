using System;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;

namespace MGroup.Constitutive.Structural.Continuum
{
	/// <summary>
	/// Deformation Gradient based implementation of a linear elastic isotropic material
	/// Authors Gerasimos Sotiropoulos
	/// </summary>
	public class ElasticMaterial3DDefGrad : IContinuumMaterial3DDefGrad
	{
		private readonly double[] strains = new double[6];
		private double[] stresses = new double[6];
		private double[,] constitutiveMatrix = null;
		public double YoungModulus { get; }
		public double PoissonRatio { get; }
		//public double[] Coordinates { get; set; }
		private readonly GenericConstitutiveLawState currentState;

		public ElasticMaterial3DDefGrad(double youngModulus, double poissonRatio)
		{
			YoungModulus = youngModulus;
			PoissonRatio = poissonRatio;
			currentState = new GenericConstitutiveLawState(this, new (string, double)[0]);
		}

		private double[,] GetConstitutiveMatrix()
		{
			double fE1 = YoungModulus / (double)(1 + PoissonRatio);
			double fE2 = fE1 * PoissonRatio / (double)(1 - 2 * PoissonRatio);
			double fE3 = fE1 + fE2;
			double fE4 = fE1 * 0.5;
			double[,] afE = new double[6, 6];
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

			Vector s = (Matrix.CreateFromArray(afE)) * (Vector.CreateFromArray(strains));
			stresses = s.CopyToArray();

			return afE;
		}

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
		public bool IsCurrentStateDifferent() => false;

		/// <summary>
		/// Resets the boolean that indicates if the constitutive matrix of the material has changed for the current iteratively update 
		/// of the deformation state.
		/// </summary>
		public void ResetModified()
		{
		}

		#endregion

		#region IFiniteElementMaterial3D Members
		/// <summary>
		/// Returns the stresses of this material for the current strain state.
		/// </summary>
		public double[] Stresses { get => stresses; }

		/// <summary>
		/// Returns the constitutive matrix of the material for the current strain state
		/// </summary>
		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (constitutiveMatrix == null) UpdateConstitutiveMatrixAndEvaluateResponse(new double[9]);
				return Matrix.CreateFromArray(constitutiveMatrix);
			}
		}

		/// <summary>
		/// Updates the material state for a given new strain point.
		/// </summary>
		/// <param name="Delta">The given strain point.</param>
		public double[] UpdateConstitutiveMatrixAndEvaluateResponse(double[] DefGradVec )
		{
			//throw new NotImplementedException();

			double[,] DGtr = new double[3, 3] { { DefGradVec [0], DefGradVec[7], DefGradVec[5] },
												{ DefGradVec [3], DefGradVec[1], DefGradVec[8] },
												{ DefGradVec [6], DefGradVec[4], DefGradVec[2] }};

			double[,] GL = new double[3, 3];
			double[] GLvec = new double[6];
			for(int m = 0; m < 3; m++)
			{ 
				for (int n = 0; n < 3; n++)
				{
					GL[m, n] = 0;
					for (int p = 0; p < 3; p++)
					{
						GL[m, n] += DGtr[m, p] * DGtr[n, p];
					}
				}
			}
			for (int m = 0; m < 3; m++)
			{
				GL[m, m] += -1;
			}
			for (int m = 0; m < 3; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					GL[m, n] = 0.5 * GL[m, n];
				}
			}

			//
			for (int m = 0; m < 3; m++)
			{
				GLvec[m] = GL[m, m];
			}
			GLvec[3] = 2 * GL[0, 1];
			GLvec[4] = 2 * GL[1, 2];
			GLvec[5] = 2 * GL[2, 0];


			double[] strains = GLvec;
			strains.CopyTo(this.strains, 0);
			constitutiveMatrix = GetConstitutiveMatrix();

			return stresses;
		}

		/// <summary>
		/// Clears the saved stress strain point connected to the last converged analysis step.
		/// Currently not a valid operation.
		/// </summary>
		public void ClearState()
		{
			//throw new NotImplementedException();
		}

		/// <summary>
		/// Saves the current stress strain state of the material (after convergence of the iterative solution process
		/// for a given loading step).
		/// </summary>
		public GenericConstitutiveLawState CreateState() => currentState;
		IHaveState ICreateState.CreateState() => CreateState();
		public GenericConstitutiveLawState CurrentState
		{
			get => currentState;
			set { }
		}

		/// <summary>
		/// Clears stresses - Currently not a valid operation.
		/// </summary>
		public void ClearStresses()
		{
			//throw new NotImplementedException();
		}

		#endregion

		#region ICloneable Members
		/// <summary>
		/// Creates a clone of material object with the same parameters.
		/// </summary>
		/// <returns>The created material clone</returns>
		public object Clone() => new ElasticMaterial3DDefGrad(YoungModulus, PoissonRatio);

		/// <summary>
		/// Creates a clone of material object with the same parameters.
		/// </summary>
		/// <returns>The created material clone</returns>
		object ICloneable.Clone() => Clone();

		#endregion

	}
}
