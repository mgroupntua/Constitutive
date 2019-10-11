using System;
using MGroup.LinearAlgebra;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Materials.Interfaces;

namespace MGroup.Materials
{
	/// <summary>
	/// It is an implementation of a linear elastic isotropic constitutive law for 2D finite elements.
	/// </summary>
	public class ElasticMaterial2D : IIsotropicContinuumMaterial2D
	{
		private readonly double[] strains = new double[3];
		private double[] stresses = new double[3];
		private Matrix constitutiveMatrix = null;

		
		public double PoissonRatio { get; set; }
		public StressState2D StressState { get; }
		public double YoungModulus { get; set; }

		/// <summary>
		/// Creastes a new object of <see cref="ElasticMaterial2D"/> class.
		/// </summary>
		/// <param name="stressState"> The stress strain state of the 2D problem.</param>
		public ElasticMaterial2D(StressState2D stressState)
		{
			this.StressState = stressState;
		}

		#region IFiniteElementMaterial3D

		/// <summary>
		/// Returns the constitutive matrix of the material for the current strain state
		/// </summary>
		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (constitutiveMatrix == null) UpdateMaterial(new double[3]);
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
		public void SaveState()
		{
			throw new NotImplementedException();
		}

		/// <summary>
		/// Updates the material state for a given new strain point.
		/// </summary>
		/// <param name="Delta">The given strain point.</param>
		public void UpdateMaterial(double[] strains)
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
			// correction for stress calculation 
			stresses = (constitutiveMatrix * Vector.CreateFromArray( strains)).CopyToArray();
			// TODO correct strains add update strains depends of how we define update strain tou ulikou

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
		public bool Modified => false;

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
		public ElasticMaterial2D Clone()
		{
			return new ElasticMaterial2D(StressState)
			{
				PoissonRatio = this.PoissonRatio,
				YoungModulus = this.YoungModulus
			};
		}

		#endregion
		/// <summary>
		/// Returns coordinates. Currently not a valid operation.
		/// </summary>
		public double[] Coordinates { get; set; }
	}
}
