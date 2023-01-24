using System;
using System.Linq;
using System.IO;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;
namespace MGroup.Constitutive.Structural.Continuum
{


	public class Tresca : IIsotropicContinuumMaterial3D
	{
		#region FieldsConstructor
		// Comments and explanations
		// Ambrosios Savvides. 
		//Fields and properties
		private const string PLASTIC_STRAIN = "Plastic strain";
		private const string STRESS_X = "Stress X";
		private const string STRESS_Y = "Stress Y";
		private const string STRESS_Z = "Stress Z";
		private const string STRESS_XY = "Stress XY";
		private const string STRESS_XZ = "Stress XZ";
		private const string STRESS_YZ = "Stress YZ";
		private GenericConstitutiveLawState currentState;
		public double Zeta { get; set; }
		public double Kmax;
		public double Kmin;
		public double[,] IsotropicHardeningCurve;
		public double[,] KinematicHardeningCurve;
		private Matrix PrConstMatr = Matrix.CreateFromArray(new double[6, 6]);
		public IMatrixView PrincipalStressConstitutiveMatrix { get; set; }
		public Matrix Rotation = Matrix.CreateFromArray(new double[6, 6]);
		public int npoi;
		public double tolerance = Math.Pow(10, -8);
		bool validitytomainplane = true;
		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (this.constitutiveMatrix == null) UpdateMaterial(new double[6]);
				return constitutiveMatrix;
			}
			set { }
		}
		private Matrix ConstMatr = Matrix.CreateFromArray(new double[6, 6]);
		private Matrix ElConstMatr = Matrix.CreateFromArray(new double[6, 6]);
		public double[] Coordinates { get; set; }
		public int ID { get; set; }
		public bool Modified { get; set; }
		public double PoissonRatio { get; set; }
		public double[] Stresses { get; set; }
		private double[] stresses = new double[6];
		public double YoungModulus { get; set; }
		public double[] tempStresses;
		public double[] Stressestrial { get; set; }
		public double[] initialStresses = new double[6];
		public readonly double shearModulus;
		private Matrix elasticConstitutiveMatrix; //the readonly was erased due to the change of the elasticconstitutivematric regarding time
		private Matrix constitutiveMatrix;
		private double[] incrementalStrains = new double[6];
		private double plasticStrain;
		private double plasticStrainNew;
		private double[] stressesNew = new double[6];
		public double yieldstress;
		public bool hasfailed = false;
		private bool modified;
		public Vector PrincipalStresses { get; set; }
		public Matrix PrincipalVectors { get; set; }
		public void UpdateMaterial(double[] strainsIncrement)
		{
			for (int j1 = 0; j1 < 6; j1++)
			{
				incrementalStrains[j1] = strainsIncrement[j1];
			}
			//this.incrementalStrains = strainsIncrement.DeepClone();
			for (int i = 0; i < 6; i++)
				this.Stresses[i] = this.tempStresses[i];
			this.DecideLoad(incrementalStrains, Stresses, ConstitutiveMatrix);
		}
		void ClearState()
		{
			modified = false;
			//this.ConstitutiveMatrix = new Matrix2D<double>(new double[6, 6]);
			Array.Clear(incrementalStrains, 0, 5);
			Array.Clear(Stresses, 0, 5);
			//Array.Clear(stressesNew, 0, stressesNew.Length);
			this.plasticStrain = 0;
			this.plasticStrainNew = 0;
		}
		public void SaveState()
		{
			//this.plasticStrain = this.plasticStrainNew;
			//Array.Copy(this.stressesNew, this.Stresses, 6);
			//for (int i = 0; i < 6; i++)
				//this.tempStresses[i] = this.Stresses[i];
			this.plasticStrain = this.plasticStrainNew;
			//if (this.plasticStrain > 0)
			//{
				//this.yieldstress = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			//}
		}

		public void ResetModified()
		{
			this.modified = false;
		}
		public GenericConstitutiveLawState CreateState()
		{
			this.plasticStrain = this.plasticStrainNew;
	        for (int i = 0; i < 6; i++)
		       this.tempStresses[i] = this.Stresses[i];
		    stresses.CopyFrom(tempStresses);
			if (this.plasticStrain > 0)
			{
				this.yieldstress = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			}
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

			return Stresses;
		}

		public object Clone()
		{
			var strainsCopy = new double[incrementalStrains.Length];
			incrementalStrains.CopyTo(strainsCopy, 0);
			var stressesCopy = new double[Stresses.Length];
			Stresses.CopyTo(stressesCopy, 0);
			this.ConstitutiveMatrix = ConstitutiveMatrix;
			//watch out if you use clone.
			Tresca m = new Tresca(this.YoungModulus,this.PoissonRatio, this.yieldstress)
			{
				modified = this.Modified,
				plasticStrain = this.plasticStrain
			};
			for (int j1 = 0; j1 < 6; j1++)
			{
				m.incrementalStrains[j1] = strainsCopy[j1];
				m.Stresses[j1] = stressesCopy[j1];
			}

			return m;
		}
		public Tresca(double youngModulus, double poissonRatio, double sigmayield)
		{
			this.tempStresses = new double[6];
			this.YoungModulus = youngModulus;
			this.yieldstress = sigmayield;
			if (poissonRatio == 0.5)
			{
				throw new ArgumentException(
					"Poisson ratio cannot be 0.5 (incompressible solid)");
			}

			this.PoissonRatio = poissonRatio;
			readMatrixData("isotropiccurve.txt", out this.IsotropicHardeningCurve);
			readMatrixData("kinematiccurve.txt", out this.KinematicHardeningCurve);
			npoi = this.IsotropicHardeningCurve.GetLength(0);

			this.shearModulus = this.YoungModulus / (2 * (1 + this.PoissonRatio));
			var value1 = this.YoungModulus / ((1 + this.PoissonRatio) * (1 - 2 * this.PoissonRatio));
			var lamda = this.PoissonRatio * value1;
			this.elasticConstitutiveMatrix = Matrix.CreateZero(6, 6);
			this.elasticConstitutiveMatrix[0, 0] = value1 * (1 - this.PoissonRatio);
			this.elasticConstitutiveMatrix[0, 1] = lamda;
			this.elasticConstitutiveMatrix[0, 2] = lamda;
			this.elasticConstitutiveMatrix[1, 0] = lamda;
			this.elasticConstitutiveMatrix[1, 1] = value1 * (1 - this.PoissonRatio);
			this.elasticConstitutiveMatrix[1, 2] = lamda;
			this.elasticConstitutiveMatrix[2, 0] = lamda;
			this.elasticConstitutiveMatrix[2, 1] = lamda;
			this.elasticConstitutiveMatrix[2, 2] = value1 * (1 - this.PoissonRatio);
			this.elasticConstitutiveMatrix[3, 3] = this.shearModulus;
			this.elasticConstitutiveMatrix[4, 4] = this.shearModulus;
			this.elasticConstitutiveMatrix[5, 5] = this.shearModulus;
			this.constitutiveMatrix = this.elasticConstitutiveMatrix;
			this.ConstitutiveMatrix = this.elasticConstitutiveMatrix;
			this.PrincipalStresses = Vector.CreateFromArray(new double[3]);
			this.PrincipalVectors = Matrix.CreateFromArray(new double[3, 3]);
			this.Stresses = new double[6];
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
		#endregion
		#region DecideLoad
		public void DecideLoad(double[] de, double[] Stresses, IMatrixView ConstitutiveMatrix)
		{
			Stressestrial = compds(de, Stresses);
			double[] deviatoricpart = GetStressDeviator(Stressestrial);
			Vector princistr = compprincipalstresses(deviatoricpart);
			var fvalue = compf(princistr, this.plasticStrain);
			if (fvalue <= 0)
			{
				for (int iii = 0; iii < 6; iii++)
				{
					Stresses[iii] = Stressestrial[iii];
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
		private Matrix BuildConsistentTangentialConstitutiveMatrix()
		{
			Matrix temp1 = PrincipalStressConstitutiveMatrix.MultiplyLeft(Rotation);
			Matrix final = temp1.MultiplyRight(Rotation.Transpose());
			Matrix check = final - final.Transpose();
			return final;
		}
		public double[] ReturnMapping(double[] Stresses, double[] StressesTrial)
		{
			Stresses = ReturntoMainPlane(Stresses, Stressestrial);
			if (validitytomainplane == false)
			{
				bool rightcorner = this.PrincipalStresses[0] + this.PrincipalStresses[2] - 2 * this.PrincipalStresses[1] > 0;
				if (rightcorner)
				{
					Stresses = ReturntoRightCorner(Stresses, Stressestrial);
				}
				else
				{
					Stresses = ReturntoLeftCorner(Stresses, Stressestrial);
				}
			}
			return Stresses;
		}
		public double[] ReturntoMainPlane(double[] Stresses, double[] StressesTrial)
		{
			double dgamma = 0;
			double yieldstr = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double bulk = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			double phi = this.PrincipalStresses[0] - this.PrincipalStresses[2] - yieldstr;
			double slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double slope = -4 * this.shearModulus - slope1;
			double dgammaprev = dgamma;
			dgamma = dgamma - phi / slope;
			double convergencerate = (dgamma - dgammaprev) / dgamma;
			bool hasconverged = Math.Abs(convergencerate) - tolerance < 0;
			while (!hasconverged)
			{
				yieldstr = GetYieldBackStressFromPlasticStrain(this.plasticStrain + dgamma, this.IsotropicHardeningCurve);
				phi = this.PrincipalStresses[0] - this.PrincipalStresses[2] - 4 * this.shearModulus * dgamma - yieldstr;
				slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain + dgamma, this.IsotropicHardeningCurve);
				slope = -4 * this.shearModulus - slope1;
				dgammaprev = dgamma;
				dgamma = dgamma - phi / slope;
				convergencerate = (dgamma - dgammaprev) / dgamma;
				hasconverged = Math.Abs(convergencerate) - tolerance < 0;
			}
			Stresses[0] = this.PrincipalStresses[0] - dgamma * (2.0 * this.shearModulus);
			Stresses[1] = this.PrincipalStresses[1];
			Stresses[2] = this.PrincipalStresses[2] + dgamma * (2.0 * this.shearModulus);
			Stresses[3] = 0;
			Stresses[4] = 0;
			Stresses[5] = 0;
			PrConstMatr.Clear();
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[1]) < 0.01 || Math.Abs(Stresses[0] - Stresses[1]) < 0.01)
			{
				PrConstMatr[3, 3] = this.shearModulus;
			}
			else
			{
				PrConstMatr[3, 3] = this.shearModulus * (Stresses[0] - Stresses[1]) / (this.PrincipalStresses[0] - this.PrincipalStresses[1]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[2]) < 0.01 || Math.Abs(Stresses[0] - Stresses[2]) < 0.01)
			{
				PrConstMatr[4, 4] = this.shearModulus;
			}
			else
			{
				PrConstMatr[4, 4] = this.shearModulus * (Stresses[0] - Stresses[2]) / (this.PrincipalStresses[0] - this.PrincipalStresses[2]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[1] - this.PrincipalStresses[2]) < 0.01 || Math.Abs(Stresses[1] - Stresses[2]) < 0.01)
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
				double p = (Stressestrial[0] + StressesTrial[1] + StressesTrial[2]) / 3;
				Stresses[0] += p;
				Stresses[1] += p;
				Stresses[2] += p;
				Stresses = RotatefromPrincipaltoCartesianStresses(Stresses, this.PrincipalStresses, this.PrincipalVectors);
			}
			this.plasticStrainNew = this.plasticStrain + dgamma;
			this.yieldstress = GetYieldBackStressFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			bool checkfyield = checkFYield(Stresses);
			double H = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			double f = 2 * this.shearModulus / (4 * this.shearModulus + H);
			PrConstMatr[0, 0] = bulk + 4 * this.shearModulus / 3 - 2 * this.shearModulus * f;
			PrConstMatr[0, 1] = bulk - 2 * this.shearModulus / 3;
			PrConstMatr[0, 2] = bulk + 2 * this.shearModulus * f - 2 * this.shearModulus / 3;
			PrConstMatr[1, 0] = bulk - 2 * this.shearModulus / 3;
			PrConstMatr[1, 1] = bulk + 4 * this.shearModulus / 3;
			PrConstMatr[1, 2] = bulk - 2 * this.shearModulus / 3;
			PrConstMatr[2, 0] = bulk - 2 * this.shearModulus / 3 + 2 * this.shearModulus * f;
			PrConstMatr[2, 1] = bulk - 2 * this.shearModulus / 3;
			PrConstMatr[2, 2] = bulk + 4 * this.shearModulus / 3 - 2 * this.shearModulus * f;
			Matrix test = PrConstMatr - PrConstMatr.Transpose();
			this.PrincipalStressConstitutiveMatrix = PrConstMatr;
			return Stresses;
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
			double p = (PrStr[0] + PrStr[1] + PrStr[2]) / 3;
			PrStr[0] = PrStr[0] - p;
			PrStr[1] = PrStr[1] - p;
			PrStr[2] = PrStr[2] - p;
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
			double p = (PrStr[0] + PrStr[1] + PrStr[2]) / 3;
			PrStr[0] = PrStr[0] - p;
			PrStr[1] = PrStr[1] - p;
			PrStr[2] = PrStr[2] - p;
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
			double p = (PrStr[0] + PrStr[1] + PrStr[2]) / 3;
			PrStr[0] = PrStr[0] - p;
			PrStr[1] = PrStr[1] - p;
			PrStr[2] = PrStr[2] - p;
			double fval1 = compf(PrStr, this.plasticStrainNew);
			bool isok1 = Math.Abs(fval1) < Math.Pow(10, -5);
			double fval2 = compf3(PrStr, this.plasticStrainNew);
			bool isok2 = Math.Abs(fval2) < Math.Pow(10, -5);
			bool isok = isok1 && isok2;
			return isok;
		}
		public double[] ReturntoRightCorner(double[] Stresses, double[] StressesTrial)
		{
			double dgammaa = 0.0;
			double dgammab = 0.0;
			double yieldstr = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double bulk = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			double sa = this.PrincipalStresses[0] - this.PrincipalStresses[2];
			double sb = this.PrincipalStresses[0] - this.PrincipalStresses[1];
			double slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			Matrix d = Matrix.CreateFromArray(new double[2, 2]);
			d[0, 0] = -4 * this.shearModulus - slope1;
			d[1, 1] = -4 * this.shearModulus - slope1;
			d[1, 0] = -2 * this.shearModulus - slope1;
			d[0, 1] = -2 * this.shearModulus - slope1;
			Matrix dinv = Matrix.CreateFromArray(new double[2, 2]);
			double det = 1.0 / (Math.Pow(d[0, 0], 2) - Math.Pow(d[0, 1], 2));
			dinv[0, 0] = det * d[1, 1];
			dinv[1, 1] = det * d[0, 0];
			dinv[1, 0] = -det * d[1, 0];
			dinv[0, 1] = -det * d[1, 0];
			double dgammaaprev = dgammaa;
			double dgammabprev = dgammab;
			double phia = sa - yieldstr;
			double phib = sb - yieldstr;
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
				yieldstr = GetYieldBackStressFromPlasticStrain(this.plasticStrain + (dgammaa + dgammab), this.IsotropicHardeningCurve);
				phia = sa - 2 * this.shearModulus * (2 * dgammaa + dgammab) - yieldstr;
				phib = sb - 2 * this.shearModulus * (dgammaa + 2 * dgammab) - yieldstr;
				slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain + (dgammaa + dgammab), this.IsotropicHardeningCurve);
				d[0, 0] = -4 * this.shearModulus - slope1;
				d[1, 1] = -4 * this.shearModulus - slope1;
				d[1, 0] = -2 * this.shearModulus - slope1;
				d[0, 1] = -2 * this.shearModulus - slope1;
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
			double p = (Stressestrial[0] + StressesTrial[1] + StressesTrial[2]) / 3;
			Stresses[0] = p + this.PrincipalStresses[0] - (dgammaa + dgammab) * 2 * this.shearModulus;
			Stresses[1] = p + this.PrincipalStresses[1] + dgammab * 2 * this.shearModulus;
			Stresses[2] = p + this.PrincipalStresses[2] + dgammaa * 2 * this.shearModulus;
			Stresses[3] = 0;
			Stresses[4] = 0;
			Stresses[5] = 0;
			PrConstMatr.Clear();
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[1]) < 0.01 || Math.Abs(Stresses[0] - Stresses[1]) < 0.01)
			{
				PrConstMatr[3, 3] = this.shearModulus;
			}
			else
			{
				PrConstMatr[3, 3] = this.shearModulus * (Stresses[0] - Stresses[1]) / (this.PrincipalStresses[0] - this.PrincipalStresses[1]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[2]) < 0.01 || Math.Abs(Stresses[0] - Stresses[2]) < 0.01)
			{
				PrConstMatr[4, 4] = this.shearModulus;
			}
			else
			{
				PrConstMatr[4, 4] = this.shearModulus * (Stresses[0] - Stresses[2]) / (this.PrincipalStresses[0] - this.PrincipalStresses[2]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[1] - this.PrincipalStresses[2]) < 0.01 || Math.Abs(Stresses[1] - Stresses[2]) < 0.01)
			{
				PrConstMatr[5, 5] = this.shearModulus;
			}
			else
			{
				PrConstMatr[5, 5] = this.shearModulus * (Stresses[1] - Stresses[2]) / (this.PrincipalStresses[1] - this.PrincipalStresses[2]); //Modification Matrix
			}
			this.plasticStrainNew = this.plasticStrain;
			bool checkfyield = checkFYield3(Stresses);
			Stresses = RotatefromPrincipaltoCartesianStresses(Stresses, this.PrincipalStresses, this.PrincipalVectors);
			this.plasticStrainNew = this.plasticStrain + (dgammaa + dgammab);
			this.yieldstress = GetYieldBackStressFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			//bool checkfyield = checkFYield3(Stresses);
			double H = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			double help = (this.shearModulus * H) / (3 * this.shearModulus + H);
			PrConstMatr[0, 0] = 4.0 / 3.0 * help + bulk;
			PrConstMatr[0, 1] = -2.0 / 3.0 * help + bulk;
			PrConstMatr[0, 2] = -2.0 / 3.0 * help + bulk;
			PrConstMatr[1, 0] = -2.0 / 3.0 * help + bulk;
			PrConstMatr[1, 1] = 1.0 / 3.0 * help + bulk;
			PrConstMatr[1, 2] = 1.0 / 3.0 * help + bulk;
			PrConstMatr[2, 0] = -2.0 / 3.0 * help + bulk;
			PrConstMatr[2, 1] = 1.0 / 3.0 * help + bulk;
			PrConstMatr[2, 2] = 1.0 / 3.0 * help + bulk;
			Matrix test = PrConstMatr - PrConstMatr.Transpose();
			this.PrincipalStressConstitutiveMatrix = PrConstMatr;
			return Stresses;
		}
		public double[] ReturntoLeftCorner(double[] Stresses, double[] StressesTrial)
		{
			double dgammaa = 0.0;
			double dgammab = 0.0;
			double yieldstr = GetYieldBackStressFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			double bulk = this.YoungModulus / (3 * (1 - 2 * this.PoissonRatio));
			double sa = this.PrincipalStresses[0] - this.PrincipalStresses[2];
			double sb = this.PrincipalStresses[1] - this.PrincipalStresses[2];
			double slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain, this.IsotropicHardeningCurve);
			Matrix d = Matrix.CreateFromArray(new double[2, 2]);
			d[0, 0] = -4 * this.shearModulus - slope1;
			d[1, 1] = -4 * this.shearModulus - slope1;
			d[1, 0] = -2 * this.shearModulus - slope1;
			d[0, 1] = -2 * this.shearModulus - slope1;
			Matrix dinv = Matrix.CreateFromArray(new double[2, 2]);
			double det = 1.0 / (Math.Pow(d[0, 0], 2) - Math.Pow(d[0, 1], 2));
			dinv[0, 0] = det * d[1, 1];
			dinv[1, 1] = det * d[0, 0];
			dinv[1, 0] = -det * d[1, 0];
			dinv[0, 1] = -det * d[1, 0];
			double dgammaaprev = dgammaa;
			double dgammabprev = dgammab;
			double phia = sa - yieldstr;
			double phib = sb - yieldstr;
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
				yieldstr = GetYieldBackStressFromPlasticStrain(this.plasticStrain + (dgammaa + dgammab), this.IsotropicHardeningCurve);
				phia = sa - 2 * this.shearModulus * (2 * dgammaa + dgammab) - yieldstr;
				phib = sb - 2 * this.shearModulus * (dgammaa + 2 * dgammab) - yieldstr;
				slope1 = GetYieldBackSlopeFromPlasticStrain(this.plasticStrain + (dgammaa + dgammab), this.IsotropicHardeningCurve);
				d[0, 0] = -4 * this.shearModulus - slope1;
				d[1, 1] = -4 * this.shearModulus - slope1;
				d[1, 0] = -2 * this.shearModulus - slope1;
				d[0, 1] = -2 * this.shearModulus - slope1;
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
			double p = (Stressestrial[0] + StressesTrial[1] + StressesTrial[2]) / 3;
			Stresses[0] = p + this.PrincipalStresses[0] - dgammaa * 2 * this.shearModulus;
			Stresses[1] = p + this.PrincipalStresses[1] - dgammab * 2 * this.shearModulus;
			Stresses[2] = p + this.PrincipalStresses[2] + (dgammab + dgammaa) * 2 * this.shearModulus;
			Stresses[3] = 0;
			Stresses[4] = 0;
			Stresses[5] = 0;
			PrConstMatr.Clear();
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[1]) < 0.01 || Math.Abs(Stresses[0] - Stresses[1]) < 0.01)
			{
				PrConstMatr[3, 3] = this.shearModulus;
			}
			else
			{
				PrConstMatr[3, 3] = this.shearModulus * (Stresses[0] - Stresses[1]) / (this.PrincipalStresses[0] - this.PrincipalStresses[1]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[0] - this.PrincipalStresses[2]) < 0.01 || Math.Abs(Stresses[0] - Stresses[2]) < 0.01)
			{
				PrConstMatr[4, 4] = this.shearModulus;
			}
			else
			{
				PrConstMatr[4, 4] = this.shearModulus * (Stresses[0] - Stresses[2]) / (this.PrincipalStresses[0] - this.PrincipalStresses[2]); //Modification Matrix
			}
			if (Math.Abs(this.PrincipalStresses[1] - this.PrincipalStresses[2]) < 0.01 || Math.Abs(Stresses[1] - Stresses[2]) < 0.01)
			{
				PrConstMatr[5, 5] = this.shearModulus;
			}
			else
			{
				PrConstMatr[5, 5] = this.shearModulus * (Stresses[1] - Stresses[2]) / (this.PrincipalStresses[1] - this.PrincipalStresses[2]); //Modification Matrix
			}
			this.plasticStrainNew = this.plasticStrain + (dgammaa + dgammab);
			bool checkfyield = checkFYield2(Stresses);
			Stresses = RotatefromPrincipaltoCartesianStresses(Stresses, this.PrincipalStresses, this.PrincipalVectors);
			this.plasticStrainNew = this.plasticStrain + (dgammaa + dgammab);
			this.yieldstress = GetYieldBackStressFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			double H = GetYieldBackSlopeFromPlasticStrain(this.plasticStrainNew, this.IsotropicHardeningCurve);
			double help = (this.shearModulus * H) / (3 * this.shearModulus + H);
			PrConstMatr[0, 0] = 1.0 / 3.0 * help + bulk;
			PrConstMatr[0, 1] = 1.0 / 3.0 * help + bulk;
			PrConstMatr[0, 2] = -2.0 / 3.0 * help + bulk;
			PrConstMatr[1, 0] = 1.0 / 3.0 * help + bulk;
			PrConstMatr[1, 1] = 1.0 / 3.0 * help + bulk;
			PrConstMatr[1, 2] = -2.0 / 3.0 * help + bulk;
			PrConstMatr[2, 0] = -2.0 / 3.0 * help + bulk;
			PrConstMatr[2, 1] = -2.0 / 3.0 * help + bulk;
			PrConstMatr[2, 2] = 4.0 / 3.0 * help + bulk;
			Matrix test = PrConstMatr - PrConstMatr.Transpose();
			this.PrincipalStressConstitutiveMatrix = PrConstMatr;
			return Stresses;
		}
		#endregion
		#region comps
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
		public bool ISEQUAL(double a, double b)
		{
			bool iseq = Math.Abs(a - b) < Math.Pow(10, -9);
			return iseq;
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
		public Vector complargestshearstress(Vector principalstresses)
		{
			Vector largestshearstress = Vector.CreateFromArray(new double[3]);
			largestshearstress[0] = Math.Abs(principalstresses[0] - principalstresses[1]) / 2;
			largestshearstress[1] = Math.Abs(principalstresses[1] - principalstresses[2]) / 2;
			largestshearstress[2] = Math.Abs(principalstresses[0] - principalstresses[2]) / 2;
			return largestshearstress;
		}
		public double compf(Vector deviatoricprincipalstresses, double plasticstrain)
		{
			yieldstress = GetYieldBackStressFromPlasticStrain(plasticstrain, this.IsotropicHardeningCurve);
			double fval = deviatoricprincipalstresses[0] - deviatoricprincipalstresses[2] - yieldstress;
			return fval;
		}
		public double compf2(Vector deviatoricprincipalstresses, double plasticstrain)
		{
			yieldstress = GetYieldBackStressFromPlasticStrain(plasticstrain, this.IsotropicHardeningCurve);
			double fval = deviatoricprincipalstresses[1] - deviatoricprincipalstresses[2] - yieldstress;
			return fval;
		}
		public double compf3(Vector deviatoricprincipalstresses, double plasticstrain)
		{
			yieldstress = GetYieldBackStressFromPlasticStrain(plasticstrain, this.IsotropicHardeningCurve);
			double fval = deviatoricprincipalstresses[0] - deviatoricprincipalstresses[1] - yieldstress;
			return fval;
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
