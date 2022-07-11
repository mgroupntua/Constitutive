using System;
using System.IO;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;

namespace MGroup.Constitutive.Structural.Continuum
{
	

	public class VonMises3DNonLinearHardening : IIsotropicContinuumMaterial3D
	{
		private const string PLASTIC_STRAIN = "Plastic strain";
		private const string STRESS_X = "Stress X";
		private const string STRESS_Y = "Stress Y";
		private const string STRESS_Z = "Stress Z";
		private const string STRESS_XY = "Stress XY";
		private const string STRESS_XZ = "Stress XZ";
		private const string STRESS_YZ = "Stress YZ";
		private GenericConstitutiveLawState currentState;
		/// <summary>
		///   The Poisson ratio value of an incompressible solid.
		/// </summary>
		private const double PoissonRatioForIncompressibleSolid = 0.5;

		/// <summary>
		///   The total number of strains.
		/// </summary>
		private const int TotalStrains = 6;

		/// <summary>
		///   The total number of stresses.
		/// </summary>
		private const int TotalStresses = TotalStrains;

		/// <summary>
		///   An array needed for the formulation of the consistent constitutive matrix.
		/// </summary>
		private static readonly double[,] SupportiveMatrixForConsistentConstitutiveMatrix = new[,]
			{
				{  2.0/3.0, -1.0/3.0, -1.0/3.0, 0,   0,   0 },
				{ -1.0/3.0, 2.0/3.0, -1.0/3.0, 0,   0,   0 },
				{ -1.0/3.0, -1.0/3.0, 2.0/3.0, 0,   0,   0  },
				{  0,  0,  0, 1.0, 0,   0   },
				{  0,  0,  0, 0,   1.0, 0   },
				{  0,  0,  0, 0,   0,   1.0 }
			};

		/// <summary>
		///   The constitutive matrix of the material while still in the elastic region.
		/// </summary>
		private readonly Matrix elasticConstitutiveMatrix;
		public static int npoi; // number of points in the isotropic and kinematic hardening curves
		/// <summary>
		///  Isotropic Hardening Curve 1 row : plastic strain 2 row: yield stress
		/// </summary>
		public double[,] IsotropicHardeningCurve;
		// Kinematic Hardening Curve 1 row : plastic strain 2 row: backstress
		/// </summary>
		public double[,] KinematicHardeningCurve;

		/// <summary>
		///   The Poisson ratio.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Poisson%27s_ratio">Wikipedia - Poisson's Ratio</a>
		/// </remarks>
		private double poissonRatio;

		/// <summary>
		///   The shear modulus.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Shear_modulus">Wikipedia - Shear Modulus</a>
		/// </remarks>
		private readonly double shearModulus;

		/// <summary>
		///   The yields stress.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Yield_%28engineering%29">Yield (engineering)</a>
		///   The current yield strength or yield point of a material is defined in engineering and materials science as the stress at which a material begins to deform plastically.
		/// </remarks>
		private double yieldStress;
		///   The initial yield strength or yield point of a material is defined in engineering and materials science as the stress at which a material begins to deform plastically.
		/// </remarks>
		private double yieldStressNew;
		/// <summary>
		///   The young modulus.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Young%27s_modulus">Wikipedia - Young's Modulus</a>
		/// </remarks>
		private double youngModulus;

		/// <summary>
		///   The constitutive matrix of the material.
		/// </summary>
		private Matrix constitutiveMatrix;

		/// <summary>
		///   The array of incremental strains.
		/// </summary>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Deformation_%28engineering%29">Deformation (engineering)</a>
		/// </remarks>
		private double[] incrementalStrains = new double[6];

		/// <summary>
		///   Indicates whether this <see cref = "IFiniteElementMaterial" /> is modified.
		/// </summary>
		private bool modified;

		/// <summary>
		///   The current plastic strain.
		/// </summary>
		private double plasticStrain;

		/// <summary>
		///   The new plastic strain.
		/// </summary>
		private double plasticStrainNew;
		/// <summary>
		///   The Scalar value of the kinematic hardening curve.
		/// </summary>
		private double kinematicscalarb;
		/// <summary>
		///   The array of stresses.
		/// </summary>
		private double[] stresses = new double[6];

		/// <summary>
		///   The array of backstresses-center of yield surface.
		/// </summary>
		private double[] backstresses = new double[6];
		/// <summary>
		///   The array of new backstresses-center of yield surface.
		/// </summary>
		private double[] backstressesnew = new double[6];
		/// <summary>
		///   The array of new stresses.
		/// </summary>
		private double[] stressesNew = new double[6];

		/// <summary>
		///   Initializes a new instance of the <see cref = "VonMisesMaterial3D" /> class.
		/// </summary>
		/// <param name = "youngModulus">
		///   The young modulus.
		/// </param>
		/// <param name = "poissonRatio">
		///   The Poisson ratio.
		/// </param>
		/// <param name = "yieldStress">
		///   The yield stress.
		/// </param>
		/// <param name = "hardeningRatio">
		///   The hardening ratio.
		/// </param>
		/// <exception cref = "ArgumentException"> When Poisson ratio is equal to 0.5.</exception>
		public VonMises3DNonLinearHardening(double youngModulus, double poissonRatio, double yieldStress)
		{
			this.youngModulus = youngModulus;

			if (poissonRatio == PoissonRatioForIncompressibleSolid)
			{
				throw new ArgumentException(
					"Poisson ratio cannot be" + PoissonRatioForIncompressibleSolid + "(incompressible solid)");
			}

			this.poissonRatio = poissonRatio;
			this.yieldStress = yieldStress;
			this.yieldStressNew = this.yieldStress;
			readMatrixData("isotropiccurve.txt", out this.IsotropicHardeningCurve);
			readMatrixData("kinematiccurve.txt", out this.KinematicHardeningCurve);
			npoi = this.IsotropicHardeningCurve.GetLength(0);
			backstresses[0] = 0.0;
			backstresses[1] = 0.0;
			backstresses[2] = 0.0;
			backstresses[3] = 0.0;
			backstresses[4] = 0.0;
			backstresses[5] = 0.0;

			this.shearModulus = this.YoungModulus / (2 * (1 + this.PoissonRatio));
			var value1 = this.youngModulus / ((1 + this.poissonRatio) * (1 - 2 * this.poissonRatio));
			var lamda = this.poissonRatio * value1;
			this.elasticConstitutiveMatrix = Matrix.CreateZero(6, 6);
			this.elasticConstitutiveMatrix[0, 0] = value1 * (1 - this.poissonRatio);
			this.elasticConstitutiveMatrix[0, 1] = lamda;
			this.elasticConstitutiveMatrix[0, 2] = lamda;
			this.elasticConstitutiveMatrix[1, 0] = lamda;
			this.elasticConstitutiveMatrix[1, 1] = value1 * (1 - this.poissonRatio);
			this.elasticConstitutiveMatrix[1, 2] = lamda;
			this.elasticConstitutiveMatrix[2, 0] = lamda;
			this.elasticConstitutiveMatrix[2, 1] = lamda;
			this.elasticConstitutiveMatrix[2, 2] = value1 * (1 - this.poissonRatio);
			this.elasticConstitutiveMatrix[3, 3] = this.shearModulus;
			this.elasticConstitutiveMatrix[4, 4] = this.shearModulus;
			this.elasticConstitutiveMatrix[5, 5] = this.shearModulus;
			this.constitutiveMatrix = this.elasticConstitutiveMatrix;
		}

		public double[] Coordinates { get; set; }

		/// <summary>
		///   Gets the constitutive matrix.
		/// </summary>
		/// <value>
		///   The constitutive matrix.
		/// </value>
		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (this.constitutiveMatrix == null) UpdateMaterial(new double[6]);
				return constitutiveMatrix;
			}
		}

		/// <summary>
		///   Gets the ID of the material.
		/// </summary>
		/// <value>
		///   The id.
		/// </value>
		public int ID => 1;

		/// <summary>
		///   Gets the incremental strains of the finite element's material.
		/// </summary>
		/// <value>
		///   The incremental strains.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Deformation_%28engineering%29">Deformation (engineering)</a>
		/// </remarks>
		public double[] IncrementalStrains => this.incrementalStrains;

		/// <summary>
		///   Gets a value indicating whether this <see cref = "IFiniteElementMaterial" /> is modified.
		/// </summary>
		/// <value>
		///   <c>true</c> if modified; otherwise, <c>false</c>.
		/// </value>
		public bool Modified => this.modified;

		/// <summary>
		///   Gets the plastic strain.
		/// </summary>
		/// <value>
		///   The plastic strain.
		/// </value>
		public double PlasticStrain => this.plasticStrain;

		/// <summary>
		///   Gets the Poisson ratio.
		/// </summary>
		/// <value>
		///   The Poisson ratio.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Poisson%27s_ratio">Wikipedia - Poisson's Ratio</a>
		/// </remarks>
		public double PoissonRatio
		{
			get
			{
				return this.poissonRatio;
			}
			set
			{
				this.poissonRatio = value;
			}
		}

		/// <summary>
		///   Gets the stresses of the finite element's material.
		/// </summary>
		/// <value>
		///   The stresses.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Stress_%28mechanics%29">Stress (mechanics)</a>
		/// </remarks>
		public double[] Stresses => this.stressesNew;

		/// <summary>
		///   Gets the Young's Modulus.
		/// </summary>
		/// <value>
		///   The young modulus.
		/// </value>
		/// <remarks>
		///   <a href = "http://en.wikipedia.org/wiki/Young%27s_modulus">Wikipedia - Young's Modulus</a>
		/// </remarks>
		public double YoungModulus
		{
			get => this.youngModulus;
			set => this.youngModulus = value;
		}

		/// <summary>
		///   Creates a new object that is a copy of the current instance.
		/// </summary>
		/// <returns>
		///   A new object that is a copy of this instance.
		/// </returns>
		public object Clone()
		{
			return new VonMises3DNonLinearHardening(this.youngModulus, this.poissonRatio, this.yieldStress)
			{
				modified = this.Modified,
				plasticStrain = this.plasticStrain,
				incrementalStrains = incrementalStrains.Copy()
			};
		}

		/// <summary>
		///   Resets the indicator of whether the material is modified.
		/// </summary>
		public void ResetModified() => this.modified = false;

		public GenericConstitutiveLawState CreateState()
		{
			this.plasticStrain = this.plasticStrainNew;
			stresses.CopyFrom(stressesNew);
			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(PLASTIC_STRAIN, plasticStrain),
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
				plasticStrain = currentState.StateValues[PLASTIC_STRAIN];
				stresses[0] = currentState.StateValues[STRESS_X];
				stresses[1] = currentState.StateValues[STRESS_Y];
				stresses[2] = currentState.StateValues[STRESS_Z];
				stresses[3] = currentState.StateValues[STRESS_XY];
				stresses[4] = currentState.StateValues[STRESS_XZ];
				stresses[5] = currentState.StateValues[STRESS_YZ];
			}
		}
		/// <summary>
		///   Updates the element's material with the provided incremental strains.
		/// </summary>
		/// <param name = "strainsIncrement">The incremental strains to use for the next step.</param>
		public double[] UpdateConstitutiveMatrixAndEvaluateResponse(double[] strainsIncrement)
		{
			incrementalStrains.CopyFrom(strainsIncrement);
			this.UpdateMaterial(strainsIncrement);

			return stressesNew;
		}

		/// <summary>
		///   Saves the state of the element's material.
		/// </summary>
		public void SaveState()
		{
			this.plasticStrain = this.plasticStrainNew;
			stresses.CopyFrom(stressesNew);
			this.backstresses = this.backstressesnew;
			this.yieldStress = this.yieldStressNew;
		}

		/// <summary>
		///   Updates the element's material with the provided incremental strains.
		/// </summary>
		/// <param name = "strainsIncrement">The incremental strains to use for the next step.</param>
		public void UpdateMaterial(double[] strainsIncrement)
		{
			incrementalStrains.CopyFrom(strainsIncrement);
			this.CalculateNextStressStrainPoint();
		}

		/// <summary>
		///   Builds the consistent tangential constitutive matrix.
		/// </summary>
		/// <param name = "vonMisesStress"> This is a constant already calculated in the calling method. </param>
		/// <remarks>
		///  
		///   
		///   
		/// </remarks>
		private Matrix BuildConsistentTangentialConstitutiveMatrix(double vonMisesStress, double[] unityvector)
		{
			Matrix consistenttangent = Matrix.CreateZero(TotalStresses, TotalStrains);
			double dgamma = this.plasticStrainNew - this.plasticStrain;
			double v1 = -dgamma * 6 * Math.Pow(this.shearModulus, 2) / vonMisesStress;
			double Hk = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.KinematicHardeningCurve);
			double Hi = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			double v2 = (dgamma / vonMisesStress - 1 / (3 * this.shearModulus + Hk + Hi)) * 6 * Math.Pow(this.shearModulus, 2);
			for (int i = 0; i < 6; i++)
			{
				for (int j = 0; j < 6; j++)
				{
					consistenttangent[i, j] = this.elasticConstitutiveMatrix[i, j] + v1 * SupportiveMatrixForConsistentConstitutiveMatrix[i, j] + v2 * unityvector[i] * unityvector[j];
				}
			}
			return consistenttangent;
		}

		/// <summary>
		///   Builds the tangential constitutive matrix.
		/// </summary>
		private void BuildTangentialConstitutiveMatrix()
		{
			this.constitutiveMatrix = Matrix.CreateZero(TotalStresses, TotalStrains);
			double invariantJ2New = this.GetDeviatorSecondStressInvariant(stressesNew);

			double value2 = (3 * this.shearModulus * this.shearModulus) /
							((GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve) + (3 * this.shearModulus)) * invariantJ2New);

			var stressDeviator = this.GetStressDeviator(stressesNew);
			for (int k1 = 0; k1 < TotalStresses; k1++)
			{
				for (int k2 = 0; k2 < TotalStresses; k2++)
				{
					this.constitutiveMatrix[k2, k1] = this.elasticConstitutiveMatrix[k2, k1] -
													  (value2 * stressDeviator[k2] * stressDeviator[k1]);
				}
			}
		}

		/// <summary>
		///   Calculates the next stress-strain point.
		/// </summary>
		/// <exception cref = "InvalidOperationException"> When the new plastic strain is less than the previous one.</exception>
		private void CalculateNextStressStrainPoint()
		{
			var stressesElastic = new double[6];
			for (int i = 0; i < 6; i++)
			{
				stressesElastic[i] = this.stresses[i] - this.backstresses[i]; //Considering the kinematic hardening back stresses
				for (int j = 0; j < 6; j++)
					stressesElastic[i] += this.elasticConstitutiveMatrix[i, j] * this.incrementalStrains[j];
			}
			double[] unityvector = GetStressDeviator(stressesElastic);
			unityvector.AxpyIntoThis(backstresses, -1.0);
			double invariantJ2Elastic = this.GetDeviatorSecondStressInvariant(stressesElastic);
			unityvector.ScaleIntoThis(Math.Sqrt(1 / (2 * invariantJ2Elastic))); //Note that the norm of the vector as considered in Von Mises Stress is the sqrt(J2*2). 
			double vonMisesStress = Math.Sqrt(3 * invariantJ2Elastic);
			if (this.plasticStrain != 0.0)
			{
				this.yieldStress = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			}
			double vonMisesStressMinusYieldStress = vonMisesStress - this.yieldStress;

			bool materialIsInElasticRegion = vonMisesStressMinusYieldStress <= 0;

			if (materialIsInElasticRegion)
			{
				this.stressesNew = stressesElastic.Axpy(backstresses, 1.0);
				this.constitutiveMatrix = this.elasticConstitutiveMatrix;
				this.plasticStrainNew = this.plasticStrain;
			}
			else
			{
				kinematicscalarb = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.KinematicHardeningCurve);
				double deltaPlasticStrain = ReturnMapping(vonMisesStress, this.shearModulus, kinematicscalarb, this.yieldStress, this.plasticStrain);
				this.plasticStrainNew = this.plasticStrain + deltaPlasticStrain;
				double kinematicscalarb2 = GetYieldBackStressFromPlasticStrain(this.plasticStrainNew, this.KinematicHardeningCurve);
				this.backstressesnew = backstresses.Axpy(unityvector, Math.Sqrt(2.0 / 3.0) * (kinematicscalarb2 - kinematicscalarb));
				stressesNew = stressesElastic.Axpy(backstresses, 1.0);
				unityvector.ScaleIntoThis(-2 * shearModulus * deltaPlasticStrain * Math.Sqrt(3.0 / 2.0));
				stressesNew.AddIntoThis(unityvector);
				unityvector = GetStressDeviator(stressesElastic);//reform the unity vector because it has been changed before.
				unityvector.AxpyIntoThis(backstresses, -1.0);
				unityvector.ScaleIntoThis(Math.Sqrt(1 / (2 * invariantJ2Elastic)));
				this.yieldStressNew = GetYieldBackStressFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
				this.constitutiveMatrix = this.BuildConsistentTangentialConstitutiveMatrix(vonMisesStress, unityvector);
			}

			if (Math.Abs(this.plasticStrainNew) < Math.Abs(this.plasticStrain))
			{
				throw new InvalidOperationException("Plastic strain cannot decrease.");
			}
			if (this.plasticStrainNew != this.plasticStrain)
			{
				this.modified = true;
			}
		}

		/// <summary>
		///   Calculates and returns the first stress invariant (I1).
		/// </summary>
		/// <returns> The first stress invariant (I1).</returns>
		public double GetFirstStressInvariant(double[] stresses) => stresses[0] + stresses[1] + stresses[2];

		/// <summary>
		///   Calculates and returns the mean hydrostatic stress.
		/// </summary>
		/// <returns> The mean hydrostatic stress.</returns>
		public double GetMeanStress(double[] stresses) => GetFirstStressInvariant(stresses) / 3.0;

		/// <summary>
		///   Calculates and returns the second stress invariant (I2).
		/// </summary>
		/// <returns> The second stress invariant (I2).</returns>
		public double GetSecondStressInvariant(double[] stresses)
			=> (stresses[0] * stresses[1]) + (stresses[1] * stresses[2]) + (stresses[0] * stresses[2])
			- Math.Pow(stresses[5], 2) - Math.Pow(stresses[3], 2) - Math.Pow(stresses[4], 2);

		/// <summary>
		///   Calculates and returns the stress deviator tensor in vector form.
		/// </summary>
		/// <returns> The stress deviator tensor in vector form.</returns>
		public double[] GetStressDeviator(double[] stresses)
		{
			var hydrostaticStress = this.GetMeanStress(stresses);
			var stressDeviator = new double[]
			{
				stresses[0] - hydrostaticStress,
				stresses[1] - hydrostaticStress,
				stresses[2] - hydrostaticStress,
				stresses[3],
				stresses[4],
				stresses[5]
			};

			return stressDeviator;
		}
		/// <summary>
		///   Calculates the yield/back stress from the plastic strain using the isotropic/kinematic hardening curve
		/// </summary>
		/// <returns> The yield stress.</returns>
		public double GetYieldBackStressFromPlasticStrain(double plasticstrain, double[,] curve)
		{
			double yieldstresstest = 0.0;

			if (plasticstrain != 0)
			{
				for (int i = 0; i < npoi - 2; i++)
				{
					if ((plasticstrain > curve[i, 0]) & (plasticstrain < curve[i + 1, 0]) & (yieldstresstest == 0.0))
					{
						double slope = (curve[i + 1, 1] - curve[i, 1]) / (curve[i + 1, 0] - curve[i, 0]);
						yieldstresstest = curve[i, 1] + slope * (plasticstrain - curve[i, 0]);
					}
				}
				if (plasticstrain > curve[npoi - 1, 0])
				{
					yieldstresstest = curve[npoi - 1, 1];
					//Console.WriteLine("Hey curve");
				}
			}
			else
			{
				yieldstresstest = curve[0, 1];
			}
			return yieldstresstest;
		}
		/// <summary>
		///   Calculates the yield/back stress slope of the isotropic/kinematic hardening curve
		/// </summary>
		/// <returns> The yield stress.</returns>
		public double GetYieldBackSlopeFromPlasticStrain(double plasticstrain, double[,] curve)
		{
			double slope = 0.0;
			if (plasticstrain != 0)
			{
				for (int i = 0; i < npoi - 2; i++)
				{
					if ((plasticstrain > curve[i, 0]) & (plasticstrain < curve[i + 1, 0]) & (slope == 0.0))
					{
						slope = (curve[i + 1, 1] - curve[i, 1]) / (curve[i + 1, 0] - curve[i, 0]);
					}
				}
				if (plasticstrain > curve[npoi - 1, 0])
				{
					slope = 0.0;
					//Console.WriteLine("Hey slope");
				}
			}
			else
			{
				slope = (curve[1, 1] - curve[0, 1]) / (curve[1, 0] - curve[0, 0]);
			}
			if (double.IsNaN(slope))
			{
				slope = 0.0;
			}
			return slope;
		}
		public double ReturnMapping(double vonmisesstress, double shearmodulus, double kinematicscalarb, double yieldstress, double plasticstrain)
		{
			double dgamma = 0.0;
			double phi = vonmisesstress - yieldstress;
			double dyieldstress = GetYieldBackSlopeFromPlasticStrain(plasticstrain, this.IsotropicHardeningCurve);
			double db = GetYieldBackSlopeFromPlasticStrain(plasticstrain, this.KinematicHardeningCurve);
			double dphi = -3 * shearmodulus - db - dyieldstress;
			double dgammaprev = dgamma;
			dgamma = dgamma - phi / dphi;
			double convergencerate = (dgamma - dgammaprev) / (dgamma);
			bool hasconverged = Math.Abs(convergencerate) - Math.Pow(10, -5) < 0;
			while (!hasconverged)
			{
				double kinematicsscalarb1 = GetYieldBackStressFromPlasticStrain(plasticstrain + dgamma, this.KinematicHardeningCurve);
				double yieldstress1 = GetYieldBackStressFromPlasticStrain(plasticstrain + dgamma, this.IsotropicHardeningCurve);
				phi = vonmisesstress - 3 * shearmodulus * dgamma + kinematicscalarb - kinematicsscalarb1 - yieldstress1;
				double slopeKinematic = GetYieldBackSlopeFromPlasticStrain(plasticstrain + dgamma, this.KinematicHardeningCurve);
				double slopeIsotropic = GetYieldBackSlopeFromPlasticStrain(plasticstrain + dgamma, this.IsotropicHardeningCurve);
				dphi = -3 * shearModulus - slopeKinematic - slopeIsotropic;
				dgammaprev = dgamma;
				dgamma = dgamma - phi / dphi;
				convergencerate = (dgamma - dgammaprev) / (dgamma);
				hasconverged = Math.Abs(convergencerate) - Math.Pow(10, -5) < 0;
			}
			return dgamma;
		}
		/// <summary>
		///   Calculates and returns the third stress invariant (I3).
		/// </summary>
		/// <returns> The third stress invariant (I3). </returns>
		public double GetThirdStressInvariant(double[] stresses)
			=> (stresses[0] * stresses[1] * stresses[2]) + (2 * stresses[5] * stresses[3] * stresses[4])
			- (Math.Pow(stresses[6], 2) * stresses[2]) - (Math.Pow(stresses[4], 2) * stresses[0])
			- (Math.Pow(stresses[5], 2) * stresses[1]); // Considering the classical FEM definition of the stress vector 

		/// <summary>
		///   Returns the first stress invariant of the stress deviator tensor (J1), which is zero.
		/// </summary>
		/// <returns> The first stress invariant of the stress deviator tensor (J1). </returns>
		public double GetDeviatorFirstStressInvariant(double[] stresses) => 0;

		/// <summary>
		///   Calculates and returns the second stress invariant of the stress deviator tensor (J2).
		/// </summary>
		/// <returns> The second stress invariant of the stress deviator tensor (J2). </returns>
		public double GetDeviatorSecondStressInvariant(double[] stresses)
		{
			double j2 = 0.0;
			j2 = 1.0 / 6.0 * (Math.Pow((stresses[0] - stresses[1]), 2) + Math.Pow((stresses[2] - stresses[1]), 2) + Math.Pow((stresses[0] - stresses[2]), 2));
			j2 += Math.Pow((stresses[3]), 2) + Math.Pow((stresses[4]), 2) + Math.Pow((stresses[5]), 2);
			return j2;
		}

		/// <summary>
		///   Calculates and returns the third stress invariant of the stress deviator tensor (J3).
		/// </summary>
		/// <returns> The third deviator stress invariant (J3). </returns>
		public double GetDeviatorThirdStressInvariant(double[] stresses)
		{
			double[] sd = GetStressDeviator(stresses);
			double j3 = GetThirdStressInvariant(sd);
			return j3;
		}
		/// <summary>
		///   Reads the data points of the curves. Stored in .txt files.
		/// </summary>
		/// <returns> The curve that is stored in the .txt file. </returns>
		public static void readMatrixData(string DataFileName, out double[,] array)
		{
			string dataLine;
			string[] dataFields;
			string[] numSeparators1 = { ":" };
			string[] numSeparators2 = { " " };
			StreamReader rStream;
			rStream = File.OpenText(DataFileName);
			int dim = 1;
			int dim1 = 1;
			dataLine = rStream.ReadLine();
			dataFields = dataLine.Split(numSeparators1, StringSplitOptions.RemoveEmptyEntries);
			dim = int.Parse(dataFields[0]);
			dataLine = rStream.ReadLine();
			dataFields = dataLine.Split(numSeparators1, StringSplitOptions.RemoveEmptyEntries);
			dim1 = int.Parse(dataFields[0]);
			double[,] array1 = new double[dim, dim1];
			for (int i = 0; i < dim; i++)
			{
				dataLine = rStream.ReadLine();
				dataFields = dataLine.Split(numSeparators1, StringSplitOptions.RemoveEmptyEntries);
				for (int j = 0; j < dim1; j++)
				{
					array1[i, j] = double.Parse(dataFields[j]);
				}
			}
			rStream.Close();
			array = array1;
		}
	}
}
