using System;
using System.Linq;
using System.IO;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;
namespace MGroup.Constitutive.Structural.Continuum
{


	public class MohrCoulomb3DNonLinearHardening : IIsotropicContinuumMaterial3D
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
		public double tolerance = Math.Pow(10, -8);
		public Vector PrincipalStresses { get; set; }
		public Matrix PrincipalVectors { get; set; }
		private double youngModulus, shearModulus, poissonRatio, cohesion, friction, dilation, cohesionNew, frictionNew, dilationNew;
		public int npoi;
		private double plasticStrain;
		private double plasticStrainNew;
		bool validitytomainplane = true;
		bool validitytoedges = true;
		public Matrix Rotation = Matrix.CreateFromArray(new double[6, 6]);
		public IMatrixView PrincipalStressConstitutiveMatrix { get; set; }
		public MohrCoulomb3DNonLinearHardening(double youngModulus, double poissonRatio, double cohesion, double friction, double dilation)
		{
			this.tempstresses = new double[6];
			this.youngModulus = youngModulus;
			this.cohesion = cohesion;
			this.friction = friction;
			this.dilation = dilation;
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
		public object Clone()
		{
			//ΜΑΥΒΕ Α PROBLEM
			//var constitutiveMatrixCopy = Matrix.CreateZero(6,6);
			//Array.Copy(constitutiveMatrix, constitutiveMatrixCopy, 36);
			//var strainsCopy = new double[incrementalStrains.Length];
			//Array.Copy(incrementalStrains, strainsCopy, incrementalStrains.Length);
			//var stressesCopy = new double[stresses.Length];
			//Array.Copy(stresses, stressesCopy, stresses.Length);

			var m = new MohrCoulomb3DNonLinearHardening(this.youngModulus, this.poissonRatio, this.cohesion, this.friction, this.dilation)
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
		private Matrix BuildConsistentTangentialConstitutiveMatrix()
		{
			Matrix temp1 = PrincipalStressConstitutiveMatrix.MultiplyLeft(Rotation);
			Matrix final = temp1.MultiplyRight(Rotation.Transpose());
			Matrix check = final - final.Transpose();
			return final;
		}
		public void CalculateNextStressStrainPoint(double[] de, double[] Stresses, IMatrixView ConstitutiveMatrix)
		{
			Stressestrial = compds(de, Stresses);
			this.PrincipalStresses = compprincipalstresses(Stressestrial);
			var fvalue = compf(PrincipalStresses, this.plasticStrain);
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
				this.constitutiveMatrix = this.BuildConsistentTangentialConstitutiveMatrix();
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
		public bool checkFYield(double[] stresses)
		{
			Matrix testtensor = Matrix.CreateFromArray(new double[3, 3]);
			testtensor[0, 0] = stresses[0];
			testtensor[1, 1] = stresses[1];
			testtensor[2, 2] = stresses[2];
			testtensor[0, 1] = stresses[3];
			testtensor[1, 0] = stresses[3];
			testtensor[1, 2] = stresses[4];
			testtensor[2, 1] = stresses[4];
			testtensor[0, 2] = stresses[5];
			testtensor[2, 0] = stresses[5];
			Vector PrStr = Vector.CreateFromArray(new double[3]);
			Matrix PrVr = Matrix.CreateFromArray(new double[3, 3]);
			(PrStr, PrVr) = testtensor.CalcEigensystemSymmetric();
			Matrix Rotation = Matrix.CreateFromArray(new double[3, 3]);
			Rotation[0, 2] = 1;
			Rotation[1, 1] = 1;
			Rotation[2, 0] = 1;
			Vector prstr = Rotation.Multiply(PrStr);
			PrStr = prstr;
			Vector temp1 = Vector.CreateFromArray(new double[3]);
			temp1[0] = PrVr[0, 0];
			temp1[1] = PrVr[1, 0];
			temp1[2] = PrVr[2, 0];
			Vector temp2 = Vector.CreateFromArray(new double[3]);
			temp2[0] = PrVr[0, 2];
			temp2[1] = PrVr[1, 2];
			temp2[2] = PrVr[2, 2];
			PrVr[0, 0] = temp2[0];
			PrVr[1, 0] = temp2[1];
			PrVr[2, 0] = temp2[2];
			PrVr[0, 2] = temp1[0];
			PrVr[1, 2] = temp1[1];
			PrVr[2, 2] = temp1[2];
			double fval = compf(PrStr, this.plasticStrainNew);
			bool isok = Math.Abs(fval) < Math.Pow(10, -5);
			return isok;
		}
		public bool checkFYield2(double[] stresses)
		{
			Matrix testtensor = Matrix.CreateFromArray(new double[3, 3]);
			testtensor[0, 0] = stresses[0];
			testtensor[1, 1] = stresses[1];
			testtensor[2, 2] = stresses[2];
			testtensor[0, 1] = stresses[3];
			testtensor[1, 0] = stresses[3];
			testtensor[1, 2] = stresses[4];
			testtensor[2, 1] = stresses[4];
			testtensor[0, 2] = stresses[5];
			testtensor[2, 0] = stresses[5];
			Vector PrStr = Vector.CreateFromArray(new double[3]);
			Matrix PrVr = Matrix.CreateFromArray(new double[3, 3]);
			(PrStr, PrVr) = testtensor.CalcEigensystemSymmetric();
			Matrix Rotation = Matrix.CreateFromArray(new double[3, 3]);
			Rotation[0, 2] = 1;
			Rotation[1, 1] = 1;
			Rotation[2, 0] = 1;
			Vector prstr = Rotation.Multiply(PrStr);
			PrStr = prstr;
			Vector temp1 = Vector.CreateFromArray(new double[3]);
			temp1[0] = PrVr[0, 0];
			temp1[1] = PrVr[1, 0];
			temp1[2] = PrVr[2, 0];
			Vector temp2 = Vector.CreateFromArray(new double[3]);
			temp2[0] = PrVr[0, 2];
			temp2[1] = PrVr[1, 2];
			temp2[2] = PrVr[2, 2];
			PrVr[0, 0] = temp2[0];
			PrVr[1, 0] = temp2[1];
			PrVr[2, 0] = temp2[2];
			PrVr[0, 2] = temp1[0];
			PrVr[1, 2] = temp1[1];
			PrVr[2, 2] = temp1[2];
			double fval1 = compf(PrStr, this.plasticStrainNew);
			bool isok1 = Math.Abs(fval1) < Math.Pow(10, -5);
			double fval2 = compf2(PrStr, this.plasticStrainNew);
			bool isok2 = Math.Abs(fval2) < Math.Pow(10, -5);
			bool isok = isok1 && isok2;
			return isok;
		}
		public bool checkFYield3(double[] stresses)
		{
			Matrix testtensor = Matrix.CreateFromArray(new double[3, 3]);
			testtensor[0, 0] = stresses[0];
			testtensor[1, 1] = stresses[1];
			testtensor[2, 2] = stresses[2];
			testtensor[0, 1] = stresses[3];
			testtensor[1, 0] = stresses[3];
			testtensor[1, 2] = stresses[4];
			testtensor[2, 1] = stresses[4];
			testtensor[0, 2] = stresses[5];
			testtensor[2, 0] = stresses[5];
			Vector PrStr = Vector.CreateFromArray(new double[3]);
			Matrix PrVr = Matrix.CreateFromArray(new double[3, 3]);
			(PrStr, PrVr) = testtensor.CalcEigensystemSymmetric();
			Matrix Rotation = Matrix.CreateFromArray(new double[3, 3]);
			Rotation[0, 2] = 1;
			Rotation[1, 1] = 1;
			Rotation[2, 0] = 1;
			Vector prstr = Rotation.Multiply(PrStr);
			PrStr = prstr;
			Vector temp1 = Vector.CreateFromArray(new double[3]);
			temp1[0] = PrVr[0, 0];
			temp1[1] = PrVr[1, 0];
			temp1[2] = PrVr[2, 0];
			Vector temp2 = Vector.CreateFromArray(new double[3]);
			temp2[0] = PrVr[0, 2];
			temp2[1] = PrVr[1, 2];
			temp2[2] = PrVr[2, 2];
			PrVr[0, 0] = temp2[0];
			PrVr[1, 0] = temp2[1];
			PrVr[2, 0] = temp2[2];
			PrVr[0, 2] = temp1[0];
			PrVr[1, 2] = temp1[1];
			PrVr[2, 2] = temp1[2];
			double fval1 = compf(PrStr, this.plasticStrainNew);
			bool isok1 = Math.Abs(fval1) < Math.Pow(10, -5);
			double fval2 = compf3(PrStr, this.plasticStrainNew);
			bool isok2 = Math.Abs(fval2) < Math.Pow(10, -5);
			bool isok = isok1 && isok2;
			return isok;
		}
		public double compf(Vector principalstress, double plasticstrain)
		{
			cohesion = GetYieldBackStressFromPlasticStrain(plasticstrain, this.IsotropicHardeningCurve);
			double fval = principalstress[0] - principalstress[2] + (principalstress[0] + principalstress[2]) * Math.Sin(friction) - 2 * cohesion * Math.Cos(friction);
			return fval;
		}
		public double compf2(Vector principalstress, double plasticstrain)
		{
			cohesion = GetYieldBackStressFromPlasticStrain(plasticstrain, this.IsotropicHardeningCurve);
			double fval = principalstress[1] - principalstress[2] + (principalstress[1] + principalstress[2]) * Math.Sin(friction) - 2 * cohesion * Math.Cos(friction);
			return fval;
		}
		public double compf3(Vector principalstress, double plasticstrain)
		{
			cohesion = GetYieldBackStressFromPlasticStrain(plasticstrain, this.IsotropicHardeningCurve);
			double fval = principalstress[0] - principalstress[1] + (principalstress[0] + principalstress[1]) * Math.Sin(friction) - 2 * cohesion * Math.Cos(friction);
			return fval;
		}
		public double[] ReturnMapping(double[] Stresses, double[] StressesTrial)
		{
			Stresses = ReturntoMainPlane(Stresses, Stressestrial);
			if (validitytomainplane == false)
			{
				bool rightcorner = (1 - Math.Sin(dilation)) * this.PrincipalStresses[0] + (1 + Math.Sin(dilation)) * this.PrincipalStresses[2] - 2 * this.PrincipalStresses[1] > 0;
				if (rightcorner)
				{
					Stresses = ReturntoRightEdge(Stresses, Stressestrial);
				}
				else
				{
					Stresses = ReturntoLeftEdge(Stresses, Stressestrial);
				}
			}
			if (validitytoedges == false)
			{
				Stresses = ReturntoApex(Stresses, Stressestrial);
			}
			return Stresses;
		}
		public double[] ReturntoApex(double[] Stresses, double[] StressesTrial)
		{
			double dgammavol = 0.0;
			double cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double residual = cohesion / Math.Tan(friction) - (this.PrincipalStresses[0] + this.PrincipalStresses[1] + this.PrincipalStresses[2]) / 3;
			double slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double bulk = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			double d = slope1 * Math.Cos(friction) / (Math.Sin(dilation) * Math.Tan(friction)) + bulk;
			double dgammavolprev = dgammavol;
			dgammavol = dgammavol - residual / d;
			double convergencerate = (dgammavol - dgammavolprev) / dgammavol;
			double p = 0.0;
			bool hasconverged = Math.Abs(convergencerate) - tolerance < 0;
			while (!hasconverged)
			{
				cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain + dgammavol * Math.Cos(friction) / Math.Sin(dilation), this.IsotropicHardeningCurve);
				p = (this.PrincipalStresses[0] + this.PrincipalStresses[1] + this.PrincipalStresses[2]) / 3 - bulk * dgammavol;
				residual = cohesion / Math.Tan(friction) - p;
				slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain + dgammavol * Math.Cos(friction) / Math.Sin(dilation), this.IsotropicHardeningCurve);
				d = slope1 * Math.Cos(friction) / (Math.Sin(dilation) * Math.Tan(friction)) + bulk;
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
			Stresses = RotatefromPrincipaltoCartesianStresses(Stresses, this.PrincipalStresses, this.PrincipalVectors);
			this.plasticStrainNew = this.plasticStrain + dgammavol * Math.Cos(friction) / Math.Sin(dilation);
			this.cohesionNew = GetYieldBackStressFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			PrConstMatr.Clear();
			PrConstMatr[0, 0] = bulk * (1 - bulk / (bulk + slope1 * Math.Cos(friction) / (Math.Sin(dilation) * Math.Tan(friction))));
			PrConstMatr[0, 1] = bulk * (1 - bulk / (bulk + slope1 * Math.Cos(friction) / (Math.Sin(dilation) * Math.Tan(friction))));
			PrConstMatr[0, 2] = bulk * (1 - bulk / (bulk + slope1 * Math.Cos(friction) / (Math.Sin(dilation) * Math.Tan(friction))));
			PrConstMatr[1, 0] = bulk * (1 - bulk / (bulk + slope1 * Math.Cos(friction) / (Math.Sin(dilation) * Math.Tan(friction))));
			PrConstMatr[1, 1] = bulk * (1 - bulk / (bulk + slope1 * Math.Cos(friction) / (Math.Sin(dilation) * Math.Tan(friction))));
			PrConstMatr[1, 2] = bulk * (1 - bulk / (bulk + slope1 * Math.Cos(friction) / (Math.Sin(dilation) * Math.Tan(friction))));
			PrConstMatr[2, 0] = bulk * (1 - bulk / (bulk + slope1 * Math.Cos(friction) / (Math.Sin(dilation) * Math.Tan(friction))));
			PrConstMatr[2, 1] = bulk * (1 - bulk / (bulk + slope1 * Math.Cos(friction) / (Math.Sin(dilation) * Math.Tan(friction))));
			PrConstMatr[2, 2] = bulk * (1 - bulk / (bulk + slope1 * Math.Cos(friction) / (Math.Sin(dilation) * Math.Tan(friction))));
			PrConstMatr[3, 3] = this.shearModulus;
			PrConstMatr[4, 4] = this.shearModulus;
			PrConstMatr[5, 5] = this.shearModulus;
			this.PrincipalStressConstitutiveMatrix = PrConstMatr;
			return Stresses;
		}
		public double[] ReturntoMainPlane(double[] Stresses, double[] StressesTrial)
		{
			double dgamma = 0;
			double cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double bulk = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			double phi = this.PrincipalStresses[0] - this.PrincipalStresses[2] + (this.PrincipalStresses[0] + this.PrincipalStresses[2]) * Math.Sin(friction) - 2 * cohesion * Math.Cos(friction);
			double slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double slope2 = -4 * this.shearModulus * (1 + Math.Sin(dilation) * Math.Sin(friction) / 3);
			double slope3 = -4 * bulk * Math.Sin(dilation) * Math.Sin(friction);
			double slope4 = -4 * slope1 * Math.Pow(Math.Cos(friction), 2);
			double slope = slope2 + slope3 + slope4;
			double dgammaprev = dgamma;
			dgamma = dgamma - phi / slope;
			double convergencerate = (dgamma - dgammaprev) / dgamma;
			bool hasconverged = Math.Abs(convergencerate) - tolerance < 0;
			while (!hasconverged)
			{
				cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain + dgamma * 2 * Math.Cos(friction), this.IsotropicHardeningCurve);
				phi = this.PrincipalStresses[0] - this.PrincipalStresses[2] + (this.PrincipalStresses[0] + this.PrincipalStresses[2]) * Math.Sin(friction);
				phi = phi - dgamma * (4 * this.shearModulus * (1 + Math.Sin(dilation) * Math.Sin(friction) / 3) + 4 * bulk * Math.Sin(dilation) * Math.Sin(friction));
				phi = phi - 2 * cohesion * Math.Cos(friction);
				slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain + dgamma * 2 * Math.Cos(friction), this.IsotropicHardeningCurve);
				slope2 = -4 * this.shearModulus * (1 + Math.Sin(dilation) * Math.Sin(friction) / 3);
				slope3 = -4 * bulk * Math.Sin(dilation) * Math.Sin(friction);
				slope4 = -4 * slope1 * Math.Pow(Math.Cos(friction), 2);
				slope = slope2 + slope3 + slope4;
				dgammaprev = dgamma;
				dgamma = dgamma - phi / slope;
				convergencerate = (dgamma - dgammaprev) / dgamma;
				hasconverged = Math.Abs(convergencerate) - tolerance < 0;
			}
			Stresses[0] = this.PrincipalStresses[0] - dgamma * (2.0 * this.shearModulus * (1.0 + Math.Sin(dilation) / 3.0) + 2.0 * bulk * Math.Sin(dilation));
			Stresses[1] = this.PrincipalStresses[1] + dgamma * Math.Sin(dilation) * (4.0 * this.shearModulus / 3.0 - 2.0 * bulk);
			Stresses[2] = this.PrincipalStresses[2] + dgamma * (2.0 * this.shearModulus * (1.0 - Math.Sin(dilation) / 3.0) - 2.0 * bulk * Math.Sin(dilation));
			Stresses[3] = 0;
			Stresses[4] = 0;
			Stresses[5] = 0;
			PrConstMatr.Clear();
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[1]) < 0.01)
			{
				PrConstMatr[3, 3] = this.shearModulus;
			}
			else
			{
				PrConstMatr[3, 3] = this.shearModulus * (Stresses[0] - Stresses[1]) / (this.PrincipalStresses[0] - this.PrincipalStresses[1]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[2]) < 0.01)
			{
				PrConstMatr[4, 4] = this.shearModulus;
			}
			else
			{
				PrConstMatr[4, 4] = this.shearModulus * (Stresses[0] - Stresses[2]) / (this.PrincipalStresses[0] - this.PrincipalStresses[2]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[1] - this.PrincipalStresses[2]) < 0.01)
			{
				PrConstMatr[5, 5] = this.shearModulus;
			}
			else
			{
				PrConstMatr[5, 5] = this.shearModulus * (Stresses[1] - Stresses[2]) / (this.PrincipalStresses[1] - this.PrincipalStresses[2]); //Modification Matrix
			}
			validitytomainplane = CheckValidityofMainPlane(Stresses);
			if (validitytomainplane)
			{
				Stresses = RotatefromPrincipaltoCartesianStresses(Stresses, this.PrincipalStresses, this.PrincipalVectors);
			}
			this.plasticStrainNew = this.plasticStrain + dgamma * 2 * Math.Cos(friction);
			this.cohesionNew = GetYieldBackStressFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			//bool checkfyield = checkFYield(Stresses);
			double H = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			double SinXSin = Math.Sin(friction) * Math.Sin(dilation);
			double ALPHA = 4.0 * this.shearModulus * (1 + SinXSin / 3.0) + 4.0 * bulk * SinXSin;
			double DENOM = -ALPHA - H * 4.0 * Math.Pow(Math.Cos(friction), 2); //plus or minus
			double B1 = (2.0 * this.shearModulus * (1.0 + Math.Sin(dilation) / 3.0) + 2.0 * bulk * Math.Sin(dilation)) / DENOM;
			double B2 = (4.0 * this.shearModulus / 3.0 - 2.0 * bulk) * Math.Sin(dilation) / DENOM;
			double B3 = (2.0 * this.shearModulus * (1.0 - Math.Sin(dilation) / 3.0) - 2.0 * bulk * Math.Sin(dilation)) / DENOM;
			PrConstMatr[0, 0] = 2 * this.shearModulus * (2.0 / 3.0 + B1 * (1 + Math.Sin(friction) / 3)) + bulk * (1.0 + 2.0 * B1 * Math.Sin(friction));
			PrConstMatr[0, 1] = 1.0 / 3.0 * (3.0 * bulk - 2.0 * this.shearModulus) * (1.0 + 2.0 * B1 * Math.Sin(friction));
			PrConstMatr[0, 2] = 2.0 * this.shearModulus * (-1.0 / 3.0 - B1 * (1.0 - Math.Sin(friction) / 3.0)) + bulk * (1.0 + 2.0 * B1 * Math.Sin(friction));
			PrConstMatr[1, 0] = 2.0 * this.shearModulus * (-1.0 / 3.0 - B2 * (1.0 + Math.Sin(friction) / 3.0)) + bulk * (1.0 - 2.0 * B2 * Math.Sin(friction));
			PrConstMatr[1, 1] = 4.0 * this.shearModulus / 3.0 * (1.0 + B2 * Math.Sin(friction)) + bulk * (1.0 - 2.0 * B2 * Math.Sin(friction));
			PrConstMatr[1, 2] = 2.0 * this.shearModulus * (-1.0 / 3.0 + B2 * (1.0 - Math.Sin(friction) / 3.0)) + bulk * (1.0 - 2.0 * B2 * Math.Sin(friction));
			PrConstMatr[2, 0] = 2.0 * this.shearModulus * (-1.0 / 3.0 - B3 * (1.0 + Math.Sin(friction) / 3.0)) + bulk * (1.0 - 2.0 * B3 * Math.Sin(friction));
			PrConstMatr[2, 1] = 1.0 / 3.0 * (3.0 * bulk - 2.0 * this.shearModulus) * (1.0 - 2.0 * B3 * Math.Sin(friction));
			PrConstMatr[2, 2] = 2.0 * this.shearModulus * (2.0 / 3.0 + B3 * (1.0 - Math.Sin(friction) / 3.0)) + bulk * (1.0 - 2.0 * B3 * Math.Sin(friction));
			Matrix testPr = Matrix.CreateFromArray(new double[6, 6]);
			double K1 = 1 - (2 * this.shearModulus * (1 + Math.Sin(dilation) / 3) + 2 * bulk * Math.Sin(dilation)) * (1 + Math.Sin(friction)) / (ALPHA + H * 4.0 * Math.Pow(Math.Cos(friction), 2));
			double K2 = (2 * this.shearModulus * (1 + Math.Sin(dilation) / 3) + 2 * bulk * Math.Sin(dilation)) * (1 - Math.Sin(friction)) / (ALPHA + H * 4.0 * Math.Pow(Math.Cos(friction), 2));
			double K3 = (4 * this.shearModulus / 3 - 2 * bulk) * Math.Sin(dilation) * (1 + Math.Sin(friction)) / (ALPHA + H * 4.0 * Math.Pow(Math.Cos(friction), 2));
			double K4 = (4 * this.shearModulus / 3 - 2 * bulk) * Math.Sin(dilation) * (-1 + Math.Sin(friction)) / (ALPHA + H * 4.0 * Math.Pow(Math.Cos(friction), 2));
			double K5 = (2 * this.shearModulus * (1 - Math.Sin(dilation) / 3) - 2 * bulk * Math.Sin(dilation)) * (1 + Math.Sin(friction)) / (ALPHA + H * 4.0 * Math.Pow(Math.Cos(friction), 2));
			double K6 = 1 + (2 * this.shearModulus * (1 - Math.Sin(dilation) / 3) - 2 * bulk * Math.Sin(dilation)) * (-1 + Math.Sin(friction)) / (ALPHA + H * 4.0 * Math.Pow(Math.Cos(friction), 2));
			double D1 = elasticConstitutiveMatrix[0, 0];
			double D2 = elasticConstitutiveMatrix[0, 1];
			testPr[0, 0] = D1 * K1 + D2 * K2;
			testPr[0, 1] = D2 * K1 + D2 * K2;
			testPr[0, 2] = D2 * K1 + D1 * K2;
			testPr[1, 0] = D1 * K3 + D2 * K4 + D2;
			testPr[1, 1] = D1 + D2 * K3 + D2 * K4;
			testPr[1, 2] = D2 * K3 + D2 + D1 * K4;
			testPr[2, 0] = D1 * K5 + D2 * K6;
			testPr[2, 1] = D2 * K5 + D2 * K6;
			testPr[2, 2] = D2 * K5 + D1 * K6;
			testPr[3, 3] = PrConstMatr[3, 3];
			testPr[4, 4] = PrConstMatr[4, 4];
			testPr[5, 5] = PrConstMatr[5, 5];
			Matrix check = PrConstMatr - testPr;
			this.PrincipalStressConstitutiveMatrix = testPr;
			return Stresses;
		}
		public double[] ReturntoRightEdge(double[] Stresses, double[] StressesTrial)
		{
			double dgammaa = 0.0;
			double dgammab = 0.0;
			double cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double bulk = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			double sa = this.PrincipalStresses[0] - this.PrincipalStresses[2] + (this.PrincipalStresses[0] + this.PrincipalStresses[2]) * Math.Sin(friction);
			double sb = this.PrincipalStresses[0] - this.PrincipalStresses[1] + (this.PrincipalStresses[0] + this.PrincipalStresses[1]) * Math.Sin(friction);
			double slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double alpha = 4 * this.shearModulus * (1 + Math.Sin(friction) * Math.Sin(dilation) / 3) + 4 * bulk * Math.Sin(friction) * Math.Sin(dilation);
			double beta = 2 * this.shearModulus * (1 + Math.Sin(friction) + Math.Sin(dilation) - Math.Sin(friction) * Math.Sin(dilation) / 3) + 4 * bulk * Math.Sin(friction) * Math.Sin(dilation);
			Matrix d = Matrix.CreateFromArray(new double[2, 2]);
			d[0, 0] = -alpha - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
			d[1, 1] = -alpha - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
			d[1, 0] = -beta - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
			d[0, 1] = -beta - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
			Matrix dinv = Matrix.CreateFromArray(new double[2, 2]);
			double det = 1.0 / (Math.Pow(d[0, 0], 2) - Math.Pow(d[0, 1], 2));
			dinv[0, 0] = det * d[1, 1];
			dinv[1, 1] = det * d[0, 0];
			dinv[1, 0] = -det * d[1, 0];
			dinv[0, 1] = -det * d[1, 0];
			double dgammaaprev = dgammaa;
			double dgammabprev = dgammab;
			double phia = sa - 2 * cohesion * Math.Cos(friction);
			double phib = sb - 2 * cohesion * Math.Cos(friction);
			Vector phi = Vector.CreateFromArray(new double[2]);
			phi[0] = phia;
			phi[1] = phib;
			Vector dinvphi = dinv.Multiply(phi);
			dgammaa = dgammaa - dinvphi[0];
			dgammab = dgammab - dinvphi[1];
			double convergencea = (dgammaa - dgammaaprev) / dgammaa;
			double convergenceb = (dgammab - dgammabprev) / dgammab;
			bool hasconverged = (Math.Abs(convergencea) - tolerance < 0) && (Math.Abs(convergenceb) - tolerance < 0);
			while (!hasconverged)
			{
				cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain + (dgammaa + dgammab) * 2 * Math.Cos(friction), this.IsotropicHardeningCurve);
				phia = sa - alpha * dgammaa - beta * dgammab - 2 * cohesion * Math.Cos(friction);
				phib = sb - beta * dgammaa - alpha * dgammab - 2 * cohesion * Math.Cos(friction);
				slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain + (dgammaa + dgammab) * 2 * Math.Cos(friction), this.IsotropicHardeningCurve);
				d[0, 0] = -alpha - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
				d[1, 1] = -alpha - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
				d[1, 0] = -beta - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
				d[0, 1] = -beta - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
				det = 1.0 / (Math.Pow(d[0, 0], 2) - Math.Pow(d[0, 1], 2));
				dinv[0, 0] = det * d[1, 1];
				dinv[1, 1] = det * d[0, 0];
				dinv[1, 0] = -det * d[1, 0];
				dinv[0, 1] = -det * d[1, 0];
				phi[0] = phia;
				phi[1] = phib;
				dgammaaprev = dgammaa;
				dgammabprev = dgammab;
				dinvphi = dinv.Multiply(phi);
				dgammaa = dgammaa - dinvphi[0];
				dgammab = dgammab - dinvphi[1];
				convergencea = (dgammaa - dgammaaprev) / dgammaa;
				convergenceb = (dgammab - dgammabprev) / dgammab;
				hasconverged = (Math.Abs(convergencea) - tolerance < 0) && (Math.Abs(convergenceb) - tolerance < 0);
			}
			double A1 = 2.0 * this.shearModulus * (1.0 - Math.Sin(dilation) / 3.0) - 2.0 * bulk * Math.Sin(dilation);
			double A2 = (4.0 * this.shearModulus / 3.0 - 2.0 * bulk) * Math.Sin(dilation);
			double B1 = 2.0 * this.shearModulus * (1.0 + Math.Sin(dilation) / 3.0) + 2.0 * bulk * Math.Sin(dilation);
			Stresses[0] = this.PrincipalStresses[0] - (dgammaa + dgammab) * B1;
			Stresses[1] = this.PrincipalStresses[1] + dgammaa * A2 + dgammab * A1;
			Stresses[2] = this.PrincipalStresses[2] + dgammaa * A1 + dgammab * A2;
			Stresses[3] = 0;
			Stresses[4] = 0;
			Stresses[5] = 0;
			PrConstMatr.Clear();
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[1]) < 0.01)
			{
				PrConstMatr[3, 3] = this.shearModulus;
			}
			else
			{
				PrConstMatr[3, 3] = this.shearModulus * (Stresses[0] - Stresses[1]) / (this.PrincipalStresses[0] - this.PrincipalStresses[1]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[2]) < 0.01)
			{
				PrConstMatr[4, 4] = this.shearModulus;
			}
			else
			{
				PrConstMatr[4, 4] = this.shearModulus * (Stresses[0] - Stresses[2]) / (this.PrincipalStresses[0] - this.PrincipalStresses[2]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[1] - this.PrincipalStresses[2]) < 0.01)
			{
				PrConstMatr[5, 5] = this.shearModulus;
			}
			else
			{
				PrConstMatr[5, 5] = this.shearModulus * (Stresses[1] - Stresses[2]) / (this.PrincipalStresses[1] - this.PrincipalStresses[2]); //Modification Matrix
			}
			//this.plasticStrainNew = this.plasticStrain + (dgammaa + dgammab) * 2 * Math.Cos(friction);
			//bool checkfyield = checkFYield3(Stresses);
			validitytoedges = CheckValidityofMainPlane(Stresses);
			if (validitytoedges)
			{
				Stresses = RotatefromPrincipaltoCartesianStresses(Stresses, this.PrincipalStresses, this.PrincipalVectors);
			}
			this.plasticStrainNew = this.plasticStrain + (dgammaa + dgammab) * 2 * Math.Cos(friction);
			this.cohesionNew = GetYieldBackStressFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			//bool checkfyield = checkFYield3(Stresses);
			double H = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			double alpha1 = (4 * this.shearModulus / 3 - 2 * bulk) * Math.Sin(dilation);
			double beta1 = 2 * this.shearModulus * (1 + Math.Sin(dilation) / 3) + 2 * bulk * Math.Sin(dilation);
			double beta2 = 2 * this.shearModulus * (1 - Math.Sin(dilation) / 3) - 2 * bulk * Math.Sin(dilation);
			double ALPHA = -alpha - 4 * H * Math.Pow(Math.Cos(friction), 2);
			double BETA = -beta - 4 * H * Math.Pow(Math.Cos(friction), 2);
			det = Math.Pow(ALPHA, 2) - Math.Pow(BETA, 2);
			PrConstMatr[0, 0] = bulk + 4 * this.shearModulus / 3 + beta1 * (2 * ALPHA - 2 * BETA) * (2 * this.shearModulus + Math.Sin(friction) * (2 * bulk + 2 * this.shearModulus / 3)) / det;
			PrConstMatr[0, 1] = bulk - 2 * this.shearModulus / 3 + beta1 * (2 * this.shearModulus * (BETA - ALPHA) + ((2 * ALPHA - 2 * BETA) * (2 * bulk + 2 * this.shearModulus / 3) + (BETA - ALPHA) * 2 * this.shearModulus) * Math.Sin(friction)) / det;
			PrConstMatr[0, 2] = bulk - 2 * this.shearModulus / 3 + beta1 * (2 * this.shearModulus * (BETA - ALPHA) + ((2 * ALPHA - 2 * BETA) * (2 * bulk + 2 * this.shearModulus / 3) + (BETA - ALPHA) * 2 * this.shearModulus) * Math.Sin(friction)) / det;
			PrConstMatr[1, 0] = bulk - 2 * this.shearModulus / 3 + (alpha1 * (BETA - ALPHA) + beta2 * (BETA - ALPHA)) * (2 * this.shearModulus + (2 * bulk + 2 * this.shearModulus / 3) * Math.Sin(friction)) / det;
			PrConstMatr[1, 1] = bulk + 4 * this.shearModulus / 3 + (alpha1 * ((2 * bulk * (BETA - ALPHA) + (BETA * 2 * this.shearModulus / 3 + ALPHA * 4 * this.shearModulus / 3)) * Math.Sin(friction) - BETA * 2 * this.shearModulus) + beta2 * (ALPHA * 2 * this.shearModulus + (2 * bulk * (BETA - ALPHA) - (ALPHA * 2 * this.shearModulus / 3 + BETA * 4 * this.shearModulus / 3)) * Math.Sin(friction))) / det;
			PrConstMatr[1, 2] = bulk - 2 * this.shearModulus / 3 + (alpha1 * ((2 * bulk * (BETA - ALPHA) - (ALPHA * 2 * this.shearModulus / 3 + BETA * 4 * this.shearModulus / 3)) * Math.Sin(friction) + ALPHA * 2 * this.shearModulus) + beta2 * ((2 * bulk * (BETA - ALPHA) + (ALPHA * 4 * this.shearModulus / 3 + BETA * 2 * this.shearModulus / 3)) * Math.Sin(friction) - BETA * 2 * this.shearModulus)) / det;
			PrConstMatr[2, 0] = bulk - 2 * this.shearModulus / 3 + ((alpha1 + beta2) * (BETA - ALPHA)) * ((2 * bulk + 2 * this.shearModulus / 3) * Math.Sin(friction) + 2 * this.shearModulus) / det;
			PrConstMatr[2, 1] = bulk - 2 * this.shearModulus / 3 + (alpha1 * (((2 * bulk * (BETA - ALPHA) - (ALPHA * 2 * this.shearModulus / 3 + BETA * 4 * this.shearModulus / 3)) * Math.Sin(friction)) + ALPHA * 2 * this.shearModulus) + beta2 * (((2 * bulk * (BETA - ALPHA) + (ALPHA * 4 * this.shearModulus / 3 + BETA * 2 * this.shearModulus / 3)) * Math.Sin(friction)) - BETA * 2 * this.shearModulus)) / det;
			PrConstMatr[2, 2] = bulk + 4 * this.shearModulus / 3 + (alpha1 * (((2 * bulk * (BETA - ALPHA) + (BETA * 2 * this.shearModulus / 3 + ALPHA * 4 * this.shearModulus / 3)) * Math.Sin(friction)) - BETA * 2 * this.shearModulus) + beta2 * (((2 * bulk * (BETA - ALPHA) - (BETA * 4 * this.shearModulus / 3 + ALPHA * 2 * this.shearModulus / 3)) * Math.Sin(friction)) + ALPHA * 2 * this.shearModulus)) / det;
			Matrix testPr = Matrix.CreateFromArray(new double[6, 6]);
			double K1 = 1 - 2 * (1 + Math.Sin(friction)) * beta1 / (BETA + ALPHA);
			double K2 = (-1 + Math.Sin(friction)) * -beta1 / (BETA + ALPHA);
			double K3 = (-1 + Math.Sin(friction)) * -beta1 / (BETA + ALPHA);
			double K4 = (1 + Math.Sin(friction)) * (alpha1 + beta2) / (BETA + ALPHA);
			double K5 = 1 + (-1 + Math.Sin(friction)) * (beta2 * ALPHA - alpha1 * BETA) / det;
			double K6 = (-1 + Math.Sin(friction)) * (alpha1 * ALPHA - beta2 * BETA) / det;
			double K7 = (1 + Math.Sin(friction)) * (beta2 + alpha1) / (BETA + ALPHA);
			double K8 = (-1 + Math.Sin(friction)) * (beta2 * ALPHA - alpha1 * BETA) / det;
			double K9 = 1 + (-1 + Math.Sin(friction)) * (alpha1 * ALPHA - beta2 * BETA) / det;
			double D1 = elasticConstitutiveMatrix[0, 0];
			double D2 = elasticConstitutiveMatrix[0, 1];
			testPr[0, 0] = D1 * K1 + D2 * K2 + D2 * K3;
			testPr[0, 1] = D2 * K1 + D1 * K2 + D2 * K3;
			testPr[0, 2] = D2 * K1 + D2 * K2 + D1 * K3;
			testPr[1, 0] = D1 * K4 + D2 * K5 + D2 * K6;
			testPr[1, 1] = D2 * K4 + D1 * K5 + D2 * K6;
			testPr[1, 2] = D2 * K4 + D2 * K5 + D1 * K6;
			testPr[2, 0] = D1 * K7 + D2 * K8 + D2 * K9;
			testPr[2, 1] = D2 * K7 + D1 * K8 + D2 * K9;
			testPr[2, 2] = D2 * K7 + D2 * K8 + D1 * K9;
			testPr[3, 3] = PrConstMatr[3, 3];
			testPr[4, 4] = PrConstMatr[4, 4];
			testPr[5, 5] = PrConstMatr[5, 5];
			Matrix check = PrConstMatr - testPr;
			this.PrincipalStressConstitutiveMatrix = testPr;
			return Stresses;
		}
		public double[] ReturntoLeftEdge(double[] Stresses, double[] StressesTrial)
		{
			double dgammaa = 0.0;
			double dgammab = 0.0;
			double cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double bulk = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			double sa = this.PrincipalStresses[0] - this.PrincipalStresses[2] + (this.PrincipalStresses[0] + this.PrincipalStresses[2]) * Math.Sin(friction);
			double sb = this.PrincipalStresses[1] - this.PrincipalStresses[2] + (this.PrincipalStresses[1] + this.PrincipalStresses[2]) * Math.Sin(friction);
			double slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double alpha = 4 * this.shearModulus * (1 + Math.Sin(friction) * Math.Sin(dilation) / 3) + 4 * bulk * Math.Sin(friction) * Math.Sin(dilation);
			double beta = 2 * this.shearModulus * (1 - Math.Sin(friction) - Math.Sin(dilation) - Math.Sin(friction) * Math.Sin(dilation) / 3) + 4 * bulk * Math.Sin(friction) * Math.Sin(dilation);
			Matrix d = Matrix.CreateFromArray(new double[2, 2]);
			d[0, 0] = -alpha - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
			d[1, 1] = -alpha - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
			d[1, 0] = -beta - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
			d[0, 1] = -beta - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
			Matrix dinv = Matrix.CreateFromArray(new double[2, 2]);
			double det = 1.0 / (Math.Pow(d[0, 0], 2) - Math.Pow(d[0, 1], 2));
			dinv[0, 0] = det * d[1, 1];
			dinv[1, 1] = det * d[0, 0];
			dinv[1, 0] = -det * d[1, 0];
			dinv[0, 1] = -det * d[1, 0];
			double dgammaaprev = dgammaa;
			double dgammabprev = dgammab;
			double phia = sa - 2 * cohesion * Math.Cos(friction);
			double phib = sb - 2 * cohesion * Math.Cos(friction);
			Vector phi = Vector.CreateFromArray(new double[2]);
			phi[0] = phia;
			phi[1] = phib;
			Vector dinvphi = dinv.Multiply(phi);
			dgammaa = dgammaa - dinvphi[0];
			dgammab = dgammab - dinvphi[1];
			double convergencea = (dgammaa - dgammaaprev) / dgammaa;
			double convergenceb = (dgammab - dgammabprev) / dgammab;
			bool hasconverged = (Math.Abs(convergencea) - tolerance < 0) && (Math.Abs(convergenceb) - tolerance < 0);
			while (!hasconverged)
			{
				cohesion = GetYieldBackStressFromPlasticStrain(this.plasticStrain + (dgammaa + dgammab) * 2 * Math.Cos(friction), this.IsotropicHardeningCurve);
				phia = sa - alpha * dgammaa - beta * dgammab - 2 * cohesion * Math.Cos(friction);
				phib = sb - beta * dgammaa - alpha * dgammab - 2 * cohesion * Math.Cos(friction);
				slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain + (dgammaa + dgammab) * 2 * Math.Cos(friction), this.IsotropicHardeningCurve);
				d[0, 0] = -alpha - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
				d[1, 1] = -alpha - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
				d[1, 0] = -beta - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
				d[0, 1] = -beta - 4 * slope1 * Math.Pow(Math.Cos(friction), 2);
				det = 1.0 / (Math.Pow(d[0, 0], 2) - Math.Pow(d[0, 1], 2));
				dinv[0, 0] = det * d[1, 1];
				dinv[1, 1] = det * d[0, 0];
				dinv[1, 0] = -det * d[1, 0];
				dinv[0, 1] = -det * d[1, 0];
				phi[0] = phia;
				phi[1] = phib;
				dgammaaprev = dgammaa;
				dgammabprev = dgammab;
				dinvphi = dinv.Multiply(phi);
				dgammaa = dgammaa - dinvphi[0];
				dgammab = dgammab - dinvphi[1];
				convergencea = (dgammaa - dgammaaprev) / dgammaa;
				convergenceb = (dgammab - dgammabprev) / dgammab;
				hasconverged = (Math.Abs(convergencea) - tolerance < 0) && (Math.Abs(convergenceb) - tolerance < 0);
			}
			double A1 = 2.0 * this.shearModulus * (1.0 + Math.Sin(dilation) / 3.0) + 2.0 * bulk * Math.Sin(dilation);
			double A2 = (4.0 * this.shearModulus / 3.0 - 2.0 * bulk) * Math.Sin(dilation);
			double B1 = 2.0 * this.shearModulus * (1.0 - Math.Sin(dilation) / 3.0) - 2.0 * bulk * Math.Sin(dilation);
			Stresses[0] = this.PrincipalStresses[0] - dgammaa * A1 + dgammab * A2;
			Stresses[1] = this.PrincipalStresses[1] + dgammaa * A2 - dgammab * A1;
			Stresses[2] = this.PrincipalStresses[2] + (dgammab + dgammaa) * B1;
			Stresses[3] = 0;
			Stresses[4] = 0;
			Stresses[5] = 0;
			PrConstMatr.Clear();
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[1]) < 0.01)
			{
				PrConstMatr[3, 3] = this.shearModulus;
			}
			else
			{
				PrConstMatr[3, 3] = this.shearModulus * (Stresses[0] - Stresses[1]) / (this.PrincipalStresses[0] - this.PrincipalStresses[1]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[2]) < 0.01)
			{
				PrConstMatr[4, 4] = this.shearModulus;
			}
			else
			{
				PrConstMatr[4, 4] = this.shearModulus * (Stresses[0] - Stresses[2]) / (this.PrincipalStresses[0] - this.PrincipalStresses[2]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[1] - this.PrincipalStresses[2]) < 0.01)
			{
				PrConstMatr[5, 5] = this.shearModulus;
			}
			else
			{
				PrConstMatr[5, 5] = this.shearModulus * (Stresses[1] - Stresses[2]) / (this.PrincipalStresses[1] - this.PrincipalStresses[2]); //Modification Matrix
			}
			//this.plasticStrainNew = this.plasticStrain + (dgammaa + dgammab) * 2 * Math.Cos(friction);
			//bool checkfyield = checkFYield2(Stresses);
			validitytoedges = CheckValidityofMainPlane(Stresses);
			if (validitytoedges)
			{
				Stresses = RotatefromPrincipaltoCartesianStresses(Stresses, this.PrincipalStresses, this.PrincipalVectors);
			}
			this.plasticStrainNew = this.plasticStrain + (dgammaa + dgammab) * 2 * Math.Cos(friction);
			this.cohesionNew = GetYieldBackStressFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			double H = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			double alpha1 = (4 * this.shearModulus / 3 - 2 * bulk) * Math.Sin(dilation);
			double beta1 = 2 * this.shearModulus * (1 + Math.Sin(dilation) / 3) + 2 * bulk * Math.Sin(dilation);
			double beta2 = 2 * this.shearModulus * (1 - Math.Sin(dilation) / 3) - 2 * bulk * Math.Sin(dilation);
			double ALPHA = -alpha - 4 * H * Math.Pow(Math.Cos(friction), 2);
			double BETA = -beta - 4 * H * Math.Pow(Math.Cos(friction), 2);
			det = Math.Pow(ALPHA, 2) - Math.Pow(BETA, 2);
			PrConstMatr[0, 0] = bulk + 4 * this.shearModulus / 3 + (beta1 * (((2 * bulk * (ALPHA - BETA) + (BETA * 4 * this.shearModulus / 3 + ALPHA * 2 * this.shearModulus / 3)) * Math.Sin(friction)) + ALPHA * 2 * this.shearModulus) + alpha1 * (((2 * bulk * (BETA - ALPHA) + (ALPHA * 4 * this.shearModulus / 3 + BETA * 2 * this.shearModulus / 3)) * Math.Sin(friction)) + BETA * 2 * this.shearModulus)) / det;
			PrConstMatr[0, 1] = bulk - 2 * this.shearModulus / 3 + (beta1 * (((2 * bulk * (ALPHA - BETA) - (BETA * 2 * this.shearModulus / 3 + ALPHA * 4 * this.shearModulus / 3)) * Math.Sin(friction)) - BETA * 2 * this.shearModulus) + alpha1 * (((2 * bulk * (BETA - ALPHA) - (ALPHA * 2 * this.shearModulus / 3 + BETA * 4 * this.shearModulus / 3)) * Math.Sin(friction)) - ALPHA * 2 * this.shearModulus)) / det;
			PrConstMatr[0, 2] = bulk - 2 * this.shearModulus / 3 + ((beta1 * (ALPHA - BETA) + alpha1 * (BETA - ALPHA)) * (((2 * bulk + 2 * this.shearModulus / 3) * Math.Sin(friction)) - 2 * this.shearModulus)) / det;
			PrConstMatr[1, 0] = bulk - 2 * this.shearModulus / 3 + (beta1 * (((2 * bulk * (ALPHA - BETA) - (ALPHA * 4 * this.shearModulus / 3 + BETA * 2 * this.shearModulus / 3)) * Math.Sin(friction)) - BETA * 2 * this.shearModulus) + alpha1 * (((2 * bulk * (BETA - ALPHA) - (BETA * 4 * this.shearModulus / 3 + ALPHA * 2 * this.shearModulus / 3)) * Math.Sin(friction)) - ALPHA * 2 * this.shearModulus)) / det;
			PrConstMatr[1, 1] = bulk + 4 * this.shearModulus / 3 + (beta1 * (((2 * bulk * (ALPHA - BETA) + (ALPHA * 2 * this.shearModulus / 3 + BETA * 4 * this.shearModulus / 3)) * Math.Sin(friction)) + ALPHA * 2 * this.shearModulus) + alpha1 * (((2 * bulk * (BETA - ALPHA) + (BETA * 2 * this.shearModulus / 3 + ALPHA * 4 * this.shearModulus / 3)) * Math.Sin(friction)) + BETA * 2 * this.shearModulus)) / det;
			PrConstMatr[1, 2] = bulk - 2 * this.shearModulus / 3 + ((beta1 * (ALPHA - BETA) + alpha1 * (BETA - ALPHA)) * (((2 * bulk + 2 * this.shearModulus / 3) * Math.Sin(friction)) - 2 * this.shearModulus)) / det;
			PrConstMatr[2, 0] = bulk - 2 * this.shearModulus / 3 + (beta2 * (((2 * bulk * (2 * BETA - 2 * ALPHA) + (ALPHA - BETA) * 4 * this.shearModulus / 3 + (BETA - ALPHA) * 2 * this.shearModulus / 3) * Math.Sin(friction)) + (BETA - ALPHA) * 2 * this.shearModulus)) / det;
			PrConstMatr[2, 1] = bulk - 2 * this.shearModulus / 3 + (beta2 * (((2 * bulk * (2 * BETA - 2 * ALPHA) + (BETA - ALPHA) * 2 * this.shearModulus / 3 + (ALPHA - BETA) * 4 * this.shearModulus / 3) * Math.Sin(friction)) + (BETA - ALPHA) * 2 * this.shearModulus)) / det;
			PrConstMatr[2, 2] = bulk + 4 * this.shearModulus / 3 + (beta2 * (2 * BETA - 2 * ALPHA) * (((2 * bulk + 2 * this.shearModulus / 3) * Math.Sin(friction)) - 2 * this.shearModulus)) / det;
			Matrix testPr = Matrix.CreateFromArray(new double[6, 6]);
			double K1 = 1 + (1 + Math.Sin(friction)) * (-beta1 * ALPHA - alpha1 * BETA) / det;
			double K2 = (1 + Math.Sin(friction)) * (beta1 * BETA + alpha1 * ALPHA) / det;
			double K3 = (-1 + Math.Sin(friction)) * (alpha1 - beta1) / (BETA + ALPHA);
			double K4 = (1 + Math.Sin(friction)) * (alpha1 * ALPHA + beta1 * BETA) / det;
			double K5 = 1 + (1 + Math.Sin(friction)) * -(alpha1 * BETA + beta1 * ALPHA) / det;
			double K6 = (-1 + Math.Sin(friction)) * (alpha1 - beta1) / (BETA + ALPHA);
			double K7 = (1 + Math.Sin(friction)) * (beta2) / (BETA + ALPHA);
			double K8 = (1 + Math.Sin(friction)) * (beta2) / (BETA + ALPHA);
			double K9 = 1 + 2 * (-1 + Math.Sin(friction)) * (beta2) / (BETA + ALPHA); ;
			double D1 = elasticConstitutiveMatrix[0, 0];
			double D2 = elasticConstitutiveMatrix[0, 1];
			testPr[0, 0] = D1 * K1 + D2 * K2 + D2 * K3;
			testPr[0, 1] = D2 * K1 + D1 * K2 + D2 * K3;
			testPr[0, 2] = D2 * K1 + D2 * K2 + D1 * K3;
			testPr[1, 0] = D1 * K4 + D2 * K5 + D2 * K6;
			testPr[1, 1] = D2 * K4 + D1 * K5 + D2 * K6;
			testPr[1, 2] = D2 * K4 + D2 * K5 + D1 * K6;
			testPr[2, 0] = D1 * K7 + D2 * K8 + D2 * K9;
			testPr[2, 1] = D2 * K7 + D1 * K8 + D2 * K9;
			testPr[2, 2] = D2 * K7 + D2 * K8 + D1 * K9;
			testPr[3, 3] = PrConstMatr[3, 3];
			testPr[4, 4] = PrConstMatr[4, 4];
			testPr[5, 5] = PrConstMatr[5, 5];
			Matrix check = PrConstMatr - testPr;
			this.PrincipalStressConstitutiveMatrix = testPr;
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
		public bool CheckValidityofMainPlane(double[] Stresses)
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
