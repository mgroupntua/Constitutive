
using System;
using System.Linq;
using System.IO;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;
namespace MGroup.Constitutive.Structural.Continuum
{


	public class DruckerPrager3DNonLinearHardening : IIsotropicContinuumMaterial3D
	{
		private const string PLASTIC_STRAIN = "Plastic strain";
		private const string STRESS_X = "Stress X";
		private const string STRESS_Y = "Stress Y";
		private const string STRESS_Z = "Stress Z";
		private const string STRESS_XY = "Stress XY";
		private const string STRESS_XZ = "Stress XZ";
		private const string STRESS_YZ = "Stress YZ";
		private GenericConstitutiveLawState currentState;
		private const double PoissonRatioForIncompressibleSolid = 0.5;
		private bool modified;
		//private readonly double[,] elasticConstitutiveMatrix;
		private int plasticRegion;
		private readonly Matrix elasticConstitutiveMatrix;
		private Matrix constitutiveMatrix;
		private double[,] constitutiveMatrixNew = new double[6, 6];
		private Matrix PrConstMatr = Matrix.CreateFromArray(new double[6, 6]);
		private double[] incrementalStrains = new double[6];
		private double[] stresses = new double[6];
		private double[] stressesNew = new double[6];
		public double[,] IsotropicHardeningCurve;
		public double[,] KinematicHardeningCurve;
		public double[] Stressestrial;
		public double[] tempstresses;
		public double ni;
		public double ksi;
		public double nipaula;
		public double J2trial;
		public double ptrial;
		private static readonly double[,] SupportiveMatrixForConsistentConstitutiveMatrix = new[,]
			{
				{  2.0/3.0, -1.0/3.0, -1.0/3.0, 0,   0,   0 },
				{ -1.0/3.0, 2.0/3.0, -1.0/3.0, 0,   0,   0 },
				{ -1.0/3.0, -1.0/3.0, 2.0/3.0, 0,   0,   0  },
				{  0,  0,  0, 1.0, 0,   0   },
				{  0,  0,  0, 0,   1.0, 0   },
				{  0,  0,  0, 0,   0,   1.0 }
			};
		public double tolerance = Math.Pow(10, -8);
		public Vector PrincipalStresses { get; set; }
		public Matrix PrincipalVectors { get; set; }
		private double youngModulus, shearModulus, poissonRatio, cohesion, friction, dilation, cohesionNew, frictionNew, dilationNew;
		public int npoi;
		private double plasticStrain;
		private double plasticStrainNew;
		bool validitytosmoothportion = true;
		bool validitytoedges = true;
		public Matrix Rotation = Matrix.CreateFromArray(new double[6, 6]);
		public IMatrixView PrincipalStressConstitutiveMatrix { get; set; }
		public DruckerPrager3DNonLinearHardening(double youngModulus, double poissonRatio, double cohesion, double friction, double dilation, string type)
		{
			this.tempstresses = new double[6];
			this.youngModulus = youngModulus;
			this.cohesion = cohesion;
			this.friction = friction;
			this.dilation = dilation;
			if (type == "Outer Cone")
			{
				this.ni = 6 * Math.Sin(friction) / (Math.Sqrt(3) * (3 - Math.Sin(friction)));
				this.nipaula = 6 * Math.Sin(dilation) / (Math.Sqrt(3) * (3 - Math.Sin(dilation)));
				this.ksi = 6 * Math.Cos(friction) / (Math.Sqrt(3) * (3 - Math.Sin(friction)));
			}
			else if (type == "Inner Cone")
			{
				this.ni = 6 * Math.Sin(friction) / (Math.Sqrt(3) * (3 + Math.Sin(friction)));
				this.nipaula = 6 * Math.Sin(dilation) / (Math.Sqrt(3) * (3 + Math.Sin(dilation)));
				this.ksi = 6 * Math.Cos(friction) / (Math.Sqrt(3) * (3 + Math.Sin(friction)));
			}
			else if (type == "Plane Strain")
			{
				this.ni = 3 * Math.Tan(friction) / (Math.Sqrt(9 + 12 * Math.Pow(Math.Tan(friction), 2)));
				this.nipaula = 3 * Math.Tan(dilation) / (Math.Sqrt(9 + 12 * Math.Pow(Math.Tan(dilation), 2)));
				this.ksi = 3 / (Math.Sqrt(9 + 12 * Math.Pow(Math.Tan(friction), 2)));
			}
			else if (type == "Common Uniaxial Failure Loads")
			{
				this.ni = 3 * Math.Sin(friction) / Math.Sqrt(3);
				this.nipaula = 3 * Math.Sin(dilation) / Math.Sqrt(3);
				this.ksi = 2 * Math.Cos(friction) / Math.Sqrt(3);
			}
			else if (type == "Biaxial Fit")
			{
				this.ni = 3 * Math.Sin(friction) / (2 * Math.Sqrt(3));
				this.nipaula = 3 * Math.Sin(dilation) / (2 * Math.Sqrt(3));
				this.ksi = 2 * Math.Cos(friction) / Math.Sqrt(3);
			}
			else if (type == "Uniaxial Fit-Outer Cone")
			{
				double fcomp = cohesion;
				double ftens = friction;
				if (fcomp == ftens)
				{
					throw new ArgumentException("Compression and Tension Stresses are equal while given type is for different strengths. Please run 'Common Uniaxial Failure Loads'");
				}
				this.friction = Math.Asin((fcomp - ftens) / (fcomp + ftens));
				this.cohesion = Math.Tan(this.friction) * fcomp * ftens / (fcomp - ftens);
				this.ni = 6 * Math.Sin(this.friction) / (Math.Sqrt(3) * (3 - Math.Sin(this.friction)));
				this.nipaula = 6 * Math.Sin(this.dilation) / (Math.Sqrt(3) * (3 - Math.Sin(this.dilation)));
				this.ksi = 6 * Math.Cos(this.friction) / (Math.Sqrt(3) * (3 - Math.Sin(this.friction)));
			}
			else if (type == "Uniaxial Fit-Inner Cone")
			{
				double fcomp = cohesion;
				double ftens = friction;
				if (fcomp == ftens)
				{
					throw new ArgumentException("Compression and Tension Stresses are equal while given type is for different strengths. Please run 'Common Uniaxial Failure Loads'");
				}
				this.friction = Math.Asin((fcomp - ftens) / (fcomp + ftens));
				this.cohesion = Math.Tan(this.friction) * fcomp * ftens / (fcomp - ftens);
				this.ni = 6 * Math.Sin(this.friction) / (Math.Sqrt(3) * (3 + Math.Sin(this.friction)));
				this.nipaula = 6 * Math.Sin(this.dilation) / (Math.Sqrt(3) * (3 + Math.Sin(this.dilation)));
				this.ksi = 6 * Math.Cos(this.friction) / (Math.Sqrt(3) * (3 + Math.Sin(this.friction)));
			}
			else if (type == "Uniaxial Fit-Plane Strain")
			{
				double fcomp = cohesion;
				double ftens = friction;
				if (fcomp == ftens)
				{
					throw new ArgumentException("Compression and Tension Stresses are equal while given type is for different strengths. Please run 'Common Uniaxial Failure Loads'");
				}
				this.friction = Math.Asin((fcomp - ftens) / (fcomp + ftens));
				this.cohesion = Math.Tan(this.friction) * fcomp * ftens / (fcomp - ftens);
				this.ni = 3 * Math.Tan(this.friction) / (Math.Sqrt(9 + 12 * Math.Pow(Math.Tan(this.friction), 2)));
				this.nipaula = 3 * Math.Tan(this.dilation) / (Math.Sqrt(9 + 12 * Math.Pow(Math.Tan(this.dilation), 2)));
				this.ksi = 3 / (Math.Sqrt(9 + 12 * Math.Pow(Math.Tan(this.friction), 2)));
			}
			else
			{
				throw new ArgumentException("Your string indicating type is not correct.");
			}
			if (poissonRatio == PoissonRatioForIncompressibleSolid)
			{
				throw new ArgumentException(
					"Poisson ratio cannot be" + PoissonRatioForIncompressibleSolid + "(incompressible solid)");
			}

			this.poissonRatio = poissonRatio;
			readMatrixData("isotropiccurve.txt", out this.IsotropicHardeningCurve);
			readMatrixData("kinematiccurve.txt", out this.KinematicHardeningCurve);
			npoi = this.IsotropicHardeningCurve.GetLength(0);

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
			this.PrincipalStresses = Vector.CreateFromArray(new double[3]);
			this.PrincipalVectors = Matrix.CreateFromArray(new double[3, 3]);
		}

		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (this.constitutiveMatrix == null) UpdateMaterial(new double[6]);
				return constitutiveMatrix;
			}
		}

		public double[] Stresses => this.stressesNew;

		//public IMatrixView ConstitutiveMatrix => Matrix.CreateFromArray(constitutiveMatrix); //TODO: this copies stuff and is not efficient.

		public void UpdateMaterial(double[] strainsIncrement)
		{
			this.incrementalStrains.CopyFrom(strainsIncrement);
			for (int ii = 0; ii < 6; ii++)
			{
				this.stresses[ii] = this.tempstresses[ii];
			}
			this.CalculateNextStressStrainPoint(strainsIncrement, this.stresses, this.ConstitutiveMatrix);
		}

		public void ClearState()
		{
			//constitutiveMatrix = new double[6, 6];
			//constitutiveMatrixNew = new double[6, 6];
			//incrementalStrains = new double[6];
			//stresses = new double[6];
			//stressesNew = new double[6];

			//var Dinv = new double[6, 6];
			//DlinElas(youngModulus, poissonRatio, 6, constitutiveMatrix, Dinv);
		}

		public void SaveState()
		{
			//Array.Copy(this.constitutiveMatrixNew, this.constitutiveMatrix, 6 * 6);
			Array.Copy(this.stressesNew, this.stresses, 6);
			for (int ii = 0; ii < 6; ii++)
			{
				this.tempstresses[ii] = this.stressesNew[ii];
			}
			this.plasticStrain = this.plasticStrainNew;
			if (this.plasticStrain > 0)
			{
				this.cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			}
		}

		public void ClearStresses()
		{
			Array.Clear(this.stresses, 0, 6);
			Array.Clear(this.stressesNew, 0, 6);
		}

		public int ID => 998;

		public bool Modified => modified;

		public void ResetModified() => modified = false;
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
		public double YoungModulus
		{
			get { return youngModulus; }
			set { throw new InvalidOperationException(); }
		}

		public double PoissonRatio
		{
			get { return poissonRatio; }
			set { throw new InvalidOperationException(); }
		}

		public double[] Coordinates { get; set; }
		public double Cohesion { get { return cohesion; } }
		public double Friction { get { return friction; } }
		public double Dilation { get { return dilation; } }

		public bool hasfailed { get; set; }

		public object Clone()
		{
			//ΜΑΥΒΕ Α PROBLEM
			//var constitutiveMatrixCopy = Matrix.CreateZero(6,6);
			//Array.Copy(constitutiveMatrix, constitutiveMatrixCopy, 36);
			//var strainsCopy = new double[incrementalStrains.Length];
			//Array.Copy(incrementalStrains, strainsCopy, incrementalStrains.Length);
			//var stressesCopy = new double[stresses.Length];
			//Array.Copy(stresses, stressesCopy, stresses.Length);

			var m = new DruckerPrager3DNonLinearHardening(YoungModulus,PoissonRatio,Cohesion,Friction,Dilation,"Inner Cone")
			{
				//modified = this.Modified,
				//constitutiveMatrix = constitutiveMatrixCopy,
				//incrementalStrains = strainsCopy,
				//stresses = stressesCopy
			};
			return m;
		}
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
		#region calculations
		private Matrix BuildConsistentTangentialConstitutiveMatrixtoApex()
		{
			Matrix temp1 = PrincipalStressConstitutiveMatrix.MultiplyLeft(Rotation);
			Matrix final = temp1.MultiplyRight(Rotation.Transpose());
			Matrix check = final - final.Transpose();
			return final;
		}
		private Matrix BuildConsistentTangentialConstitutiveMatrix()
		{
			double dgamma = this.plasticStrainNew - this.plasticStrain;
			Vector deviator = Vector.CreateFromArray(GetStressDeviator(Stressestrial));
			double edevelastic = deviator.Norm2() / (2 * this.shearModulus);
			double H = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			double bulk = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			double A = 1 / (this.shearModulus + bulk * this.ni * this.nipaula + Math.Pow(this.ksi, 2) * H);
			Vector unityvector = Vector.CreateFromArray(Stressestrial);
			double norm = unityvector.Norm2();
			unityvector.ScaleIntoThis(1 / norm);
			Vector Identity = Vector.CreateFromArray(new double[6]);
			Identity[0] = 1;
			Identity[1] = 1;
			Identity[2] = 1;
			Matrix IdentityTensorIdentity = Matrix.CreateFromArray(new double[6, 6]);
			IdentityTensorIdentity[0, 0] = 1;
			IdentityTensorIdentity[0, 1] = 1;
			IdentityTensorIdentity[0, 2] = 1;
			IdentityTensorIdentity[1, 0] = 1;
			IdentityTensorIdentity[1, 1] = 1;
			IdentityTensorIdentity[1, 2] = 1;
			IdentityTensorIdentity[2, 0] = 1;
			IdentityTensorIdentity[2, 1] = 1;
			IdentityTensorIdentity[2, 2] = 1;
			Matrix temp1 = Matrix.CreateFromArray(SupportiveMatrixForConsistentConstitutiveMatrix);
			Matrix temp2 = unityvector.TensorProduct(unityvector);
			Matrix temp3 = unityvector.TensorProduct(Identity);
			Matrix temp4 = Identity.TensorProduct(unityvector);
			temp1.ScaleIntoThis(-dgamma / (Math.Sqrt(2) * edevelastic));
			temp2.ScaleIntoThis(2 * this.shearModulus * (dgamma / (Math.Sqrt(2) * edevelastic) - this.shearModulus * A));
			temp3.ScaleIntoThis(-Math.Sqrt(2) * this.shearModulus * A * bulk * this.ni);
			temp4.ScaleIntoThis(-Math.Sqrt(2) * this.shearModulus * A * bulk * this.nipaula);
			IdentityTensorIdentity.ScaleIntoThis(bulk * (-bulk * this.ni * this.nipaula * A));
			Matrix final = this.elasticConstitutiveMatrix + temp1 + temp2 + temp3 + temp4 + IdentityTensorIdentity;
			Matrix test = final - final.Transpose();
			return final;
		}
		public void CalculateNextStressStrainPoint(double[] de, double[] Stresses, IMatrixView ConstitutiveMatrix)
		{
			Stressestrial = compds(de, Stresses);
			this.ptrial = (Stressestrial[0] + Stressestrial[1] + Stressestrial[2]) / 3;
			this.J2trial = GetJ2(Stressestrial);
			var fvalue = compf(this.J2trial, this.ptrial, this.plasticStrain);
			if (fvalue <= 0)
			{
				for (int iii = 0; iii < 6; iii++)
				{
					stressesNew[iii] = Stressestrial[iii];
				}
				this.constitutiveMatrix = this.elasticConstitutiveMatrix;
			}
			else
			{
				stressesNew = ReturnMapping(Stresses, Stressestrial);
				for (int iii = 0; iii < 6; iii++)
				{
					Stresses[iii] = stressesNew[iii];
				}
			}
		}
		public double[] compds(double[] de, double[] Stresses)
		{
			double[] output = new double[6];
			for (int iii = 0; iii < 6; iii++)
			{
				var help = 0.0;
				for (int jjj = 0; jjj < 6; jjj++)
				{
					help += this.elasticConstitutiveMatrix[iii, jjj] * de[jjj];
				}
				output[iii] = help + Stresses[iii];
			}
			return output;
		}
		public double GetFirstStressInvariant(double[] stresses) => stresses[0] + stresses[1] + stresses[2];
		/// <summary>
		///   Calculates and returns the second stress invariant (I2).
		/// </summary>
		/// <returns> The second stress invariant (I2).</returns>
		public double GetSecondStressInvariant(double[] stresses)
			=> (stresses[0] * stresses[1]) + (stresses[1] * stresses[2]) + (stresses[0] * stresses[2])
			- Math.Pow(stresses[5], 2) - Math.Pow(stresses[3], 2) - Math.Pow(stresses[4], 2);
		/// <summary>
		///   Calculates and returns the third stress invariant (I3).
		/// </summary>
		/// <returns> The third stress invariant (I3). </returns>
		public double GetThirdStressInvariant(double[] stresses)
			=> (stresses[0] * stresses[1] * stresses[2]) + 2 * (stresses[3] * stresses[4] * stresses[5]) - stresses[0] * Math.Pow(stresses[4], 2) - stresses[1] * Math.Pow(stresses[5], 2) - stresses[2] * Math.Pow(stresses[3], 2); // Considering the classical FEM definition of the stress vector 

		public double GetJ2(double[] stresses)
		{
			double J2 = 0.0;
			double I1 = GetFirstStressInvariant(stresses);
			double I2 = GetSecondStressInvariant(stresses);
			J2 = Math.Pow(I1, 2) / 3 - I2;
			return J2;
		}
		public Vector compprincipalstresses(double[] stresses)
		{
			Matrix stressestensor = Matrix.CreateFromArray(new double[3, 3]);
			stressestensor[0, 0] = stresses[0];
			stressestensor[1, 1] = stresses[1];
			stressestensor[2, 2] = stresses[2];
			stressestensor[0, 1] = stresses[3];
			stressestensor[1, 0] = stresses[3];
			stressestensor[1, 2] = stresses[4];
			stressestensor[2, 1] = stresses[4];
			stressestensor[0, 2] = stresses[5];
			stressestensor[2, 0] = stresses[5];
			(this.PrincipalStresses, this.PrincipalVectors) = stressestensor.CalcEigensystemSymmetric();
			Matrix Rotation = Matrix.CreateFromArray(new double[3, 3]);
			Rotation[0, 2] = 1;
			Rotation[1, 1] = 1;
			Rotation[2, 0] = 1;
			Vector prstr = Rotation.Multiply(this.PrincipalStresses);
			this.PrincipalStresses = prstr;
			Vector temp1 = Vector.CreateFromArray(new double[3]);
			temp1[0] = this.PrincipalVectors[0, 0];
			temp1[1] = this.PrincipalVectors[1, 0];
			temp1[2] = this.PrincipalVectors[2, 0];
			Vector temp2 = Vector.CreateFromArray(new double[3]);
			temp2[0] = this.PrincipalVectors[0, 2];
			temp2[1] = this.PrincipalVectors[1, 2];
			temp2[2] = this.PrincipalVectors[2, 2];
			this.PrincipalVectors[0, 0] = temp2[0];
			this.PrincipalVectors[1, 0] = temp2[1];
			this.PrincipalVectors[2, 0] = temp2[2];
			this.PrincipalVectors[0, 2] = temp1[0];
			this.PrincipalVectors[1, 2] = temp1[1];
			this.PrincipalVectors[2, 2] = temp1[2];
			return this.PrincipalStresses;
		}
		public double GetMeanStress(double[] stresses) => GetFirstStressInvariant(stresses) / 3.0;
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
		public double compf(double J2, double p, double plasticstrain)
		{
			cohesion = GetYieldBackStressFromPlasticStrain(plasticstrain, this.IsotropicHardeningCurve);
			double fval = Math.Sqrt(J2) + this.ni * p - this.ksi * cohesion;
			return fval;
		}
		public double[] ReturnMapping(double[] Stresses, double[] StressesTrial)
		{
			Stresses = ReturntoSmoothPortion(Stresses, Stressestrial);
			this.constitutiveMatrix = BuildConsistentTangentialConstitutiveMatrix();
			if (validitytosmoothportion == false)
			{
				Stresses = ReturntoApex(Stresses, Stressestrial);
				this.PrincipalStresses = compprincipalstresses(StressesTrial);
				double[] help = RotatefromPrincipaltoCartesianStresses(StressesTrial, this.PrincipalStresses, this.PrincipalVectors);
				this.constitutiveMatrix = BuildConsistentTangentialConstitutiveMatrixtoApex();
			}
			return Stresses;
		}
		public double[] ReturntoApex(double[] Stresses, double[] StressesTrial)
		{
			double dgammavol = 0.0;
			double cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double alpha = this.ksi / this.ni;
			double beta = this.ksi / this.nipaula;
			double residual = cohesion * beta - this.ptrial;
			double slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double bulk = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			double d = alpha * beta * slope1 + bulk;
			double dgammavolprev = dgammavol;
			dgammavol = dgammavol - residual / d;
			double convergencerate = (dgammavol - dgammavolprev) / dgammavol;
			double p = 0.0;
			bool hasconverged = Math.Abs(convergencerate) - tolerance < 0;
			while (!hasconverged)
			{
				cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain + dgammavol * alpha, this.IsotropicHardeningCurve);
				p = this.ptrial - bulk * dgammavol;
				residual = beta * cohesion - p;
				slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain + dgammavol * alpha, this.IsotropicHardeningCurve);
				d = alpha * beta * slope1 + bulk;
				dgammavolprev = dgammavol;
				dgammavol = dgammavol - residual / d;
				convergencerate = (dgammavol - dgammavolprev) / dgammavol;
				hasconverged = Math.Abs(convergencerate) - tolerance < 0;
			}
			Stresses[0] = p;
			Stresses[1] = p;
			Stresses[2] = p;
			Stresses[3] = 0;
			Stresses[4] = 0;
			Stresses[5] = 0;
			this.plasticStrainNew = this.plasticStrain + dgammavol * alpha;
			this.cohesionNew = GetYieldBackStressFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			PrConstMatr.Clear();
			PrConstMatr[0, 0] = bulk * (1 - bulk / (bulk + alpha * beta * slope1));
			PrConstMatr[0, 1] = bulk * (1 - bulk / (bulk + alpha * beta * slope1));
			PrConstMatr[0, 2] = bulk * (1 - bulk / (bulk + alpha * beta * slope1));
			PrConstMatr[1, 0] = bulk * (1 - bulk / (bulk + alpha * beta * slope1));
			PrConstMatr[1, 1] = bulk * (1 - bulk / (bulk + alpha * beta * slope1));
			PrConstMatr[1, 2] = bulk * (1 - bulk / (bulk + alpha * beta * slope1));
			PrConstMatr[2, 0] = bulk * (1 - bulk / (bulk + alpha * beta * slope1));
			PrConstMatr[2, 1] = bulk * (1 - bulk / (bulk + alpha * beta * slope1));
			PrConstMatr[2, 2] = bulk * (1 - bulk / (bulk + alpha * beta * slope1));
			PrConstMatr[3, 3] = this.shearModulus;
			PrConstMatr[4, 4] = this.shearModulus;
			PrConstMatr[5, 5] = this.shearModulus;
			this.PrincipalStressConstitutiveMatrix = PrConstMatr;
			return Stresses;
		}
		public double[] ReturntoSmoothPortion(double[] Stresses, double[] StressesTrial)
		{
			double dgamma = 0;
			double cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double bulk = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			double phi = compf(this.J2trial, this.ptrial, this.plasticStrain);
			double slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double slope = -this.shearModulus - bulk * this.ni * this.nipaula - Math.Pow(this.ksi, 2) * slope1;
			double dgammaprev = dgamma;
			dgamma = dgamma - phi / slope;
			double convergencerate = (dgamma - dgammaprev) / dgamma;
			bool hasconverged = Math.Abs(convergencerate) - tolerance < 0;
			while (!hasconverged)
			{
				cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain + dgamma * this.ksi, this.IsotropicHardeningCurve);
				phi = Math.Sqrt(GetJ2(StressesTrial)) - this.shearModulus * dgamma + this.ni * (this.ptrial - bulk * this.nipaula * dgamma) - this.ksi * cohesion;
				slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain + dgamma * this.ksi, this.IsotropicHardeningCurve);
				slope = -this.shearModulus - bulk * this.ni * this.nipaula - Math.Pow(this.ksi, 2) * slope1;
				dgammaprev = dgamma;
				dgamma = dgamma - phi / slope;
				convergencerate = (dgamma - dgammaprev) / dgamma;
				hasconverged = Math.Abs(convergencerate) - tolerance < 0;
			}
			double[] deviatoricpart = GetStressDeviator(StressesTrial);
			double p = this.ptrial - bulk * this.nipaula * dgamma;
			Stresses[0] = deviatoricpart[0] * (1 - this.shearModulus * dgamma / Math.Sqrt(GetJ2(StressesTrial))) + p;
			Stresses[1] = deviatoricpart[1] * (1 - this.shearModulus * dgamma / Math.Sqrt(GetJ2(StressesTrial))) + p;
			Stresses[2] = deviatoricpart[2] * (1 - this.shearModulus * dgamma / Math.Sqrt(GetJ2(StressesTrial))) + p;
			Stresses[3] = deviatoricpart[3] * (1 - this.shearModulus * dgamma / Math.Sqrt(GetJ2(StressesTrial)));
			Stresses[4] = deviatoricpart[4] * (1 - this.shearModulus * dgamma / Math.Sqrt(GetJ2(StressesTrial)));
			Stresses[5] = deviatoricpart[5] * (1 - this.shearModulus * dgamma / Math.Sqrt(GetJ2(StressesTrial)));
			PrConstMatr.Clear();
			validitytosmoothportion = (Math.Sqrt(GetJ2(StressesTrial)) - this.shearModulus * dgamma >= 0);
			this.plasticStrainNew = this.plasticStrain + dgamma * this.ksi;
			this.cohesionNew = GetYieldBackStressFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			double check = compf(GetJ2(Stresses), p, this.plasticStrainNew);
			bool checkfyield = Math.Abs(check) <= tolerance;
			double H = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			//Constitutive Matrix needed
			return Stresses;
		}
		public double[] RotatefromPrincipaltoCartesianStresses(double[] Stresses, Vector PrincipalStresses, Matrix PrincipalVectors)
		{
			Vector vector1 = Vector.CreateFromArray(new double[3]);
			Vector vector2 = Vector.CreateFromArray(new double[3]);
			Vector vector3 = Vector.CreateFromArray(new double[3]);
			vector1[0] = this.PrincipalVectors[0, 0];
			vector1[1] = this.PrincipalVectors[1, 0];
			vector1[2] = this.PrincipalVectors[2, 0];
			vector2[0] = this.PrincipalVectors[0, 1];
			vector2[1] = this.PrincipalVectors[1, 1];
			vector2[2] = this.PrincipalVectors[2, 1];
			vector3[0] = this.PrincipalVectors[0, 2];
			vector3[1] = this.PrincipalVectors[1, 2];
			vector3[2] = this.PrincipalVectors[2, 2];
			//if (vector1[0] < 0) vector1.ScaleIntoThis(-1);
			//if (vector2[1] < 0) vector2.ScaleIntoThis(-1);
			//if (vector3[2] < 0) vector3.ScaleIntoThis(-1);
			//Matrix s1 = vector1.TensorProduct(vector1);
			//Matrix s2 = vector2.TensorProduct(vector2);
			//Matrix s3 = vector3.TensorProduct(vector3);
			//Matrix stresstensor = Matrix.CreateFromArray(new double[3, 3]);
			//stresstensor = Stresses[0] * s1 + Stresses[1] * s2 + Stresses[2] * s3;
			double l1 = vector1[0] / vector1.Norm2();
			double l2 = vector1[1] / vector1.Norm2();
			double l3 = vector1[2] / vector1.Norm2();
			double m1 = vector2[0] / vector2.Norm2();
			double m2 = vector2[1] / vector2.Norm2();
			double m3 = vector2[2] / vector2.Norm2();
			double n1 = vector3[0] / vector3.Norm2();
			double n2 = vector3[1] / vector3.Norm2();
			double n3 = vector3[2] / vector3.Norm2();
			Matrix Ms = Matrix.CreateFromArray(new double[6, 6]);
			//1st row
			Ms[0, 0] = Math.Pow(l1, 2);
			Ms[0, 1] = Math.Pow(m1, 2);
			Ms[0, 2] = Math.Pow(n1, 2);
			Ms[0, 3] = 2 * l1 * m1;
			Ms[0, 4] = 2 * n1 * m1;
			Ms[0, 5] = 2 * l1 * n1;
			//2nd row
			Ms[1, 0] = Math.Pow(l2, 2);
			Ms[1, 1] = Math.Pow(m2, 2);
			Ms[1, 2] = Math.Pow(n2, 2);
			Ms[1, 3] = 2 * l2 * m2;
			Ms[1, 4] = 2 * n2 * m2;
			Ms[1, 5] = 2 * l2 * n2;
			//3rd row
			Ms[2, 0] = Math.Pow(l3, 2);
			Ms[2, 1] = Math.Pow(m3, 2);
			Ms[2, 2] = Math.Pow(n3, 2);
			Ms[2, 3] = 2 * l3 * m3;
			Ms[2, 4] = 2 * n3 * m3;
			Ms[2, 5] = 2 * l3 * n3;
			//4th row
			Ms[3, 0] = l1 * l2;
			Ms[3, 1] = m1 * m2;
			Ms[3, 2] = n1 * n2;
			Ms[3, 3] = l1 * m2 + l2 * m1;
			Ms[3, 4] = l1 * n2 + l2 * n1;
			Ms[3, 5] = m1 * n2 + m2 * n1;
			//5th row
			Ms[4, 0] = l2 * l3;
			Ms[4, 1] = m2 * m3;
			Ms[4, 2] = n2 * n3;
			Ms[4, 3] = l2 * m3 + l3 * m2;
			Ms[4, 4] = l2 * n3 + l3 * n2;
			Ms[4, 5] = m2 * n3 + m3 * n2;
			//6th row
			Ms[5, 0] = l1 * l3;
			Ms[5, 1] = m1 * m3;
			Ms[5, 2] = n1 * n3;
			Ms[5, 3] = l1 * m3 + l3 * m1;
			Ms[5, 4] = l1 * n3 + l3 * n1;
			Ms[5, 5] = m1 * n3 + m3 * n1;
			this.Rotation = Ms;
			Stresses = Ms.Multiply(Stresses);
			return Stresses;
		}
		public bool CheckValidityofSmoothPortion(double[] Stresses)
		{
			bool iseq1 = ISEQUAL(Stresses[0], Stresses[1]);
			bool iseq2 = ISEQUAL(Stresses[1], Stresses[2]);
			bool iseq3 = ISEQUAL(Stresses[0], Stresses[2]);
			bool validity1 = (Stresses[0] > Stresses[1]) || (iseq1);
			bool validity2 = (Stresses[1] >= Stresses[2]) || (iseq2);
			bool validity3 = (Stresses[0] >= Stresses[2]) || (iseq3);
			bool validity = validity1 && validity2 && validity3;
			return validity;
		}
		public bool ISEQUAL(double a, double b)
		{
			bool iseq = Math.Abs(a - b) < Math.Pow(10, -9);
			return iseq;
		}
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

		#endregion
	}
}
