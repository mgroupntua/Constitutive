using System;
using System.Linq;
using System.IO;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;

namespace MGroup.Constitutive.Structural.Continuum
{
	public class TsaiWu : IIsotropicContinuumMaterial3D
	{
		#region FieldsConstructor
		// Comments and explanations
		// Ambrosios Savvides. For composite materials
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
		public IMatrixView ConstitutiveMatrix { get { return ElConstMatr; } set { } }
		private Matrix ConstMatr = Matrix.CreateFromArray(new double[6, 6]);
		private Matrix ElConstMatr = Matrix.CreateFromArray(new double[6, 6]);
		public double[] Coordinates { get; set; }
		public int ID { get; set; }
		public bool Modified { get; set; }
		public double PoissonRatio { get; set; }
		public double[] Stresses { get; set; }
		public double YoungModulus { get; set; }
		public double[] tempStresses;
		public double[] Stressestrial { get; set; }
		public double[] initialStresses = new double[6];
		public readonly double shearModulus;
		public IMatrixView elasticConstitutiveMatrix { get { return ElConstMatr; } set { } } //the readonly was erased due to the change of the elasticconstitutivematric regarding time
		private double[] incrementalStrains = new double[6];
		private double plasticStrain;
		private double plasticStrainNew;
		private double[] stressesNew = new double[6];
		private double[] stresses = new double[6];
		public double F11;
		public double F22;
		public double F33;
		public double F12;
		public double F23;
		public double F31;
		public double F44;
		public double F55;
		public double F66;
		public double F1;
		public double F2;
		public double F3;
		public bool hasfailed = false;
		private bool modified;
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
			//Array.Clear(plasticStrain, 0, 5);
			//Array.Clear(plasticStrainNew, 0, 5);
		}
		public void SaveState()
		{
			//this.plasticStrain = this.plasticStrainNew;
			//Array.Copy(this.stressesNew, this.Stresses, 6);
			for (int i = 0; i < 6; i++)
				this.tempStresses[i] = this.Stresses[i];
			this.ConstitutiveMatrix = ConstitutiveMatrix;
		}

		public void ResetModified()
		{
			this.modified = false;
		}
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
			for (int ii = 0; ii < 6; ii++)
			{
				stressesNew[ii] = Stresses[ii];
			}
			return stressesNew;
		}

		public void ClearStresses()
		{
			Array.Clear(Stresses, 0, 5);
			Array.Clear(stressesNew, 0, 5);
		}

		public object Clone()
		{
			var strainsCopy = new double[incrementalStrains.Length];
			incrementalStrains.CopyTo(strainsCopy, 0);
			var stressesCopy = new double[Stresses.Length];
			Stresses.CopyTo(stressesCopy, 0);
			this.ConstitutiveMatrix = ConstitutiveMatrix;
			//watch out if you use clone.
			TsaiWu m = new TsaiWu(new double[6], new double[6], new double[6])
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
		public TsaiWu(double[] youngModuli, double[] poissonRatioi, double[] sigmas)
		{
			ElConstMatr = Matrix.CreateFromArray(new double[6, 6]);
			ElConstMatr.MatrixSymmetry = LinearAlgebra.Providers.MatrixSymmetry.Symmetric;
			var dee1 = youngModuli[0];
			var dee2 = youngModuli[1];
			var dee3 = youngModuli[2];
			var g23 = youngModuli[3];
			var g13 = youngModuli[4];
			var g12 = youngModuli[5];
			var ni12 = poissonRatioi[0];
			var ni13 = poissonRatioi[1];
			var ni23 = poissonRatioi[2];
			var ni21 = ni12 * dee2 / dee1;
			var ni31 = ni13 * dee3 / dee1;
			var ni32 = ni23 * dee3 / dee2;
			var nf = ni12 * ni21 + ni23 * ni32 + ni31 * ni13 + 2 * ni21 * ni32 * ni13;
			ElConstMatr[0, 0] = dee1 * (1 - ni23 * ni32) / (1 - nf);
			ElConstMatr[0, 1] = dee2 * (ni12 + ni32 * ni13) / (1 - nf);
			ElConstMatr[0, 2] = dee3 * (ni13 + ni12 * ni23) / (1 - nf);
			ElConstMatr[1, 0] = ElConstMatr[0, 1];
			ElConstMatr[1, 1] = dee2 * (1 - ni13 * ni31) / (1 - nf);
			ElConstMatr[1, 2] = dee3 * (ni23 + ni21 * ni13) / (1 - nf);
			ElConstMatr[2, 0] = ElConstMatr[0, 2];
			ElConstMatr[2, 1] = ElConstMatr[1, 2];
			ElConstMatr[2, 2] = dee3 * (1 - ni12 * ni21) / (1 - nf); ;
			ElConstMatr[3, 3] = g12;
			ElConstMatr[4, 4] = g23;
			ElConstMatr[5, 5] = g13;
			this.elasticConstitutiveMatrix = ElConstMatr;
			this.ConstitutiveMatrix = ElConstMatr;
			var s1t = sigmas[0]; // 1--X 2--Y 3--Z tension
			var s2t = sigmas[1];
			var s3t = sigmas[2];
			var s1c = sigmas[6]; // 1--X 2--Y 3--Z compression
			var s2c = sigmas[7];
			var s3c = sigmas[8];
			var s12 = sigmas[3];
			var s23 = sigmas[4];
			var s31 = sigmas[5];
			var F12norm = sigmas[9];  // Normalized 2D coefficients. 
			var F23norm = sigmas[10];
			var F31norm = sigmas[11];
			F11 = 1 / (s1t * s1c);
			F1 = 1 / s1t - 1 / s1c;
			F22 = 1 / (s2t * s2c);
			F2 = 1 / s2t - 1 / s2c;
			F33 = 1 / (s3t * s3c);
			F3 = 1 / s3t - 1 / s3c;
			F44 = 1 / (Math.Pow(s12, 2));
			F55 = 1 / (Math.Pow(s23, 2));
			F66 = 1 / (Math.Pow(s31, 2));
			F12 = F12norm * Math.Sqrt(F11 * F22);
			F23 = F23norm * Math.Sqrt(F22 * F33);
			F31 = F31norm * Math.Sqrt(F11 * F33);
			this.tempStresses = new double[6];
			this.Stresses = new double[6];
		}

		#endregion

		#region DecideLoad
		public void DecideLoad(double[] de, double[] Stresses, IMatrixView ConstitutiveMatrix)
		{
			if (hasfailed == false)
			{
				Stressestrial = compds(de, Stresses);
				var fvalue = compf(Stressestrial);
				if (fvalue <= 0)
				{
					for (int iii = 0; iii < 6; iii++)
					{
						Stresses[iii] = Stressestrial[iii];
					}
				}
				else
				{
					Stresses = compfail(Stressestrial, Stresses);
					this.elasticConstitutiveMatrix = Matrix.CreateFromArray(new double[6, 6]);
					this.ConstitutiveMatrix = Matrix.CreateFromArray(new double[6, 6]);
					hasfailed = true;
				}
			}
		}
		#endregion

		#region comps
		public double[] compds(double[] de, double[] Stresses)
		{
			for (int iii = 0; iii < 6; iii++)
			{
				var help = 0.0;
				for (int jjj = 0; jjj < 6; jjj++)
				{
					help += this.ConstitutiveMatrix[iii, jjj] * de[jjj];
				}
				Stresses[iii] += help;
			}
			return Stresses;
		}
		public double compf(double[] Stresses)
		{
			double fval = (F11) * Math.Pow(Stresses[0], 2) + (F22) * Math.Pow(Stresses[1], 2) + (F33) * Math.Pow(Stresses[2], 2);
			fval += Stresses[0] * Stresses[1] * (2 * F12) + Stresses[0] * Stresses[2] * (2 * F31) + Stresses[1] * Stresses[2] * (2 * F23);
			fval += (F55) * Math.Pow(Stresses[4], 2) + (F66) * Math.Pow(Stresses[5], 2) + (F44) * Math.Pow(Stresses[3], 2);
			fval += F1 * Stresses[0] + F2 * Stresses[1] + F3 * Stresses[2];
			fval = fval - 1;
			return fval;
		}
		public double[] compfail(double[] Stressestrial, double[] Stresses)
		{
			double lambdafail = 0.0;
			double A1 = Stressestrial[0] - Stresses[0];
			double A2 = Stressestrial[1] - Stresses[1];
			double A3 = Stressestrial[2] - Stresses[2];
			double A4 = Stressestrial[3] - Stresses[3];
			double A5 = Stressestrial[4] - Stresses[4];
			double A6 = Stressestrial[5] - Stresses[5];
			double B1 = Stresses[0];
			double B2 = Stresses[1];
			double B3 = Stresses[2];
			double B4 = Stresses[3];
			double B5 = Stresses[4];
			double B6 = Stresses[5];
			double A = (F11) * Math.Pow(A1, 2) + (F22) * Math.Pow(A2, 2) + (F33) * Math.Pow(A3, 2);
			A += A1 * A2 * (2 * F12) + A1 * A3 * (2 * F31) + A2 * A3 * (2 * F23);
			A += (F55) * Math.Pow(A5, 2) + (F66) * Math.Pow(A6, 2) + (F44) * Math.Pow(A4, 2);
			double B = (F11) * 2 * A1 * B1 + (F22) * 2 * A2 * B2 + (F33) * 2 * A3 * B3;
			B += (A1 * B2 + A2 * B1) * (2 * F12) + (A1 * B3 + B1 * A3) * (2 * F31) + (A2 * B3 + A3 * B2) * (2 * F23);
			B += (F55) * 2 * A5 * B5 + (F66) * 2 * A6 * B6 + (F44) * 2 * A4 * B4;
			B += F1 * A1 + F2 * A2 + F3 * A3;
			double C = (F11) * Math.Pow(B1, 2) + (F22) * Math.Pow(B2, 2) + (F33) * Math.Pow(B3, 2);
			C += B1 * B2 * (2 * F12) + B1 * B3 * (2 * F31) + B2 * B3 * (2 * F23);
			C += (F55) * Math.Pow(B5, 2) + (F66) * Math.Pow(B6, 2) + (F44) * Math.Pow(B4, 2);
			C += F1 * B1 + F2 * B2 + F3 * B3 - 1;
			var lambda1 = (-B + Math.Sqrt(Math.Pow(B, 2) - 4 * A * C)) / (2 * A);
			var lambda2 = (-B - Math.Sqrt(Math.Pow(B, 2) - 4 * A * C)) / (2 * A);
			lambdafail = lambda1;
			Stresses[0] = lambdafail * (A1) + B1;
			Stresses[1] = lambdafail * (A2) + B2;
			Stresses[2] = lambdafail * (A3) + B3;
			Stresses[3] = lambdafail * (A4) + B4;
			Stresses[4] = lambdafail * (A5) + B5;
			Stresses[5] = lambdafail * (A6) + B6;
			return Stresses;
		}
		#endregion
	}
}
