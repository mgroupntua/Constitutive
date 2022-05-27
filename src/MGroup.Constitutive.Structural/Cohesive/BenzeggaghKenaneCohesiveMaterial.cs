using System;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;

namespace MGroup.Constitutive.Structural.Cohesive
{
	/// <summary>
	/// Benzeggagh-Kenane Mixed mode damage propagation criterion. Bilinear traction separation law. 
	/// Authors Gerasimos Sotiropoulos
	/// </summary>
	public class BenzeggaghKenaneCohesiveMaterial : ICohesiveZoneMaterial
	{
		private static readonly string PLASTIC_STRAIN = "Plastic strain";

		private GenericConstitutiveLawState currentState;
		private bool modified;

		public double T_o_3 { get; }
		public double D_o_3 { get; }
		public double D_f_3 { get; }
		public double T_o_1 { get; }
		public double D_o_1 { get; }
		public double D_f_1 { get; }
		public double n_curve { get; }

		private double tol;
		private double d_prev_step;

		public BenzeggaghKenaneCohesiveMaterial(double to3, double do3, double df3, double to1, double do1, double df1, double nCurve)
		{
			D_o_3 = do3;
			D_f_3 = df3;
			T_o_1 = to1;
			D_o_1 = do1;
			D_f_1 = df1;
			T_o_3 = to3;
			n_curve = nCurve;

			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(PLASTIC_STRAIN, 0d),
			});
		}

		/// <summary>
		/// Creates a clone of material object with the same parameters.
		/// </summary>
		/// <returns>The created material clone</returns>
		//ICohesiveZoneMaterial3D ICohesiveZoneMaterial3D.Clone()
		//{
		//	return this.Clone(); 
		//}

		/// <summary>
		/// Creates a clone of material object with the same parameters.
		/// </summary>
		/// <returns>The created material clone</returns>
		object ICloneable.Clone() => Clone();

		/// <summary>
		/// Creates a clone of material object with the same parameters.
		/// </summary>
		/// <returns>The created material clone</returns>
		public BenzeggaghKenaneCohesiveMaterial Clone() => new BenzeggaghKenaneCohesiveMaterial(T_o_3, D_o_3, D_f_3, T_o_1, D_o_1, D_f_1, n_curve);

		private double Delta_s;
		private double lamda;
		private double veta;
		private double B;
		private double D_o;
		private double T_o;
		private double E;
		private double D_f;
		private double D_o_old;
		private double d;

		private Matrix D_tan;
		private double[,] D_tan_prev;
		private double[] T_int;
		private double[,] D_tan_f;

		private double H;
		private double EH;

		private bool matrices_not_initialized = true;
		public void InitializeMatrices()
		{
			D_tan = Matrix.CreateZero(3, 3);
			D_tan_prev = new double[3, 3];
			T_int = new double[3];
			D_tan_f = new double[3, 3];
			matrices_not_initialized = false;
			tol = Math.Pow(10, -19);
		}

		/// <summary>
		/// Updates the material state for a given new strain point.
		/// </summary>
		/// <param name="Delta">The given strain point.</param>
		public double[] UpdateConstitutiveMatrixAndEvaluateResponse(double[] Delta)
		{
			if (matrices_not_initialized)
			{ this.InitializeMatrices(); }

			// Store previous Tangent moduli for ifmaterialsModified check
			for (int k = 0; k < 3; k++)
			{
				for (int j = 0; j < 3; j++)
				{
					D_tan_prev[k, j] = D_tan[k, j];
				}
			}


			Delta_s = Math.Sqrt(Delta[0] * Delta[0] + Delta[1] * Delta[1]);

			if (Delta[2] < tol)
			{
				lamda = Delta_s;
				veta = 1;
			}
			else
			{
				lamda = Math.Sqrt(Delta[0] * Delta[0] + Delta[1] * Delta[1] + Delta[2] * Delta[2]);
				veta = Delta_s / (Delta_s + 0.5 * Delta[2] + Math.Abs(Delta[2]));
			}
			B = Math.Pow(veta, 2) / (1 + 2 * Math.Pow(veta, 2) - 2 * veta);
			D_o = Math.Sqrt(Math.Pow(D_o_3, 2) + Math.Pow(B, n_curve) * (Math.Pow(D_o_1, 2) - Math.Pow(D_o_3, 2)));
			T_o = Math.Sqrt(Math.Pow(T_o_3, 2) + Math.Pow(B, n_curve) * (Math.Pow(T_o_1, 2) - Math.Pow(T_o_3, 2)));
			E = T_o / D_o;
			D_f = (D_f_3 * D_o_3 + (D_f_1 * D_o_1 - D_f_3 * D_o_3) * Math.Pow(B, n_curve)) / D_o;
			D_o_old = (D_o * D_f) / (D_f - d_prev_step * (D_f - D_o));

			if (lamda < D_o_old)
			{
				d = d_prev_step;
				for (int k = 0; k < 3; k++)
				{
					for (int j = 0; j < 3; j++)
					{
						D_tan[k, j] = 0;
					}
				}
				for (int j = 0; j < 3; j++)
				{
					D_tan[j, j] = (1 - d_prev_step) * E;
				}
				if (Delta[2] < 0)
				{
					D_tan[2, 2] += d_prev_step * E;
				}
				for (int j = 0; j < 3; j++)
				{
					T_int[j] = D_tan[j, j] * Delta[j];
				}
			}
			else
			{
				d = (D_f * (lamda - D_o)) / (lamda * (D_f - D_o));
				if (d>1)
				{ d = 1; }
				for (int k = 0; k < 3; k++)
				{
					for (int j = 0; j < 3; j++)
					{
						D_tan_f[k, j] = 0;
					}
				}
				for (int j = 0; j < 3; j++)
				{
					D_tan_f[j, j] = (1 - d) * E;
				}
				if(Delta[2]<0)
				{ D_tan_f[2, 2] += d * E; }
				for (int j = 0; j < 3; j++)
				{
					T_int[j] = D_tan_f[j, j] * Delta[j];
				}


				for (int k = 0; k < 3; k++)
				{
					for (int j = 0; j < 3; j++)
					{
						D_tan[k, j] = 0;
					}
				}
				for (int j = 0; j < 3; j++)
				{
					D_tan[j, j] = (1 - d) * E;
				}

				H = (D_f * D_o) / ((D_f - D_o) * Math.Pow(lamda, 3));
				EH = E * H;
				D_tan[0, 0] += - EH * Math.Pow(Delta[0], 2);
				D_tan[1, 1] += - EH * Math.Pow(Delta[1], 2);
				if (Delta[2]>0)
				{
					D_tan[2,2]+=-EH * Math.Pow(Delta[2], 2);
					D_tan[0, 2] += -EH * Delta[0] * Delta[2];
					D_tan[1, 2] += -EH * Delta[1] * Delta[2];
					D_tan[0, 1] += -EH * Delta[0] * Delta[1];
					D_tan[2, 0] = D_tan[0, 2];
					D_tan[2, 1] = D_tan[1, 2];
					D_tan[1, 0] = D_tan[0, 1];
				}
				else
				{
					D_tan[2, 2] += d * E;
				}
			}
			this.modified = CheckIfConstitutiveMatrixChanged();

			return T_int;
		}

		private bool CheckIfConstitutiveMatrixChanged()
		{
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					if (Math.Abs(D_tan_prev[i, j] - D_tan[i, j]) > 1e-10)
						return true;

			return false;
		}

		/// <summary>
		/// Returns the tractions of this material for the current strain state.
		/// </summary>
		public double[] Tractions 
		{
			get { return T_int; }
		}

		/// <summary>
		/// Returns the constitutive matrix of the material for the current strain state
		/// </summary>
		public IMatrixView ConstitutiveMatrix
		{
			get
			{
				if (D_tan == null) UpdateConstitutiveMatrixAndEvaluateResponse(new double[3]);
				return D_tan;
			}
		}

		/// <summary>
		/// Saves the current stress strain state of the material (after convergence of the iterative solution process
		/// for a given loading step).
		/// </summary>
		public GenericConstitutiveLawState CreateState()
		{
			d_prev_step = d;
			currentState = new GenericConstitutiveLawState(this, new[]
			{
				(PLASTIC_STRAIN, d),
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
				d_prev_step = currentState.StateValues[PLASTIC_STRAIN];
			}
		}

		/// <summary>
		/// Returns a boolean indicating if the constitutive matrix of the material has changed for the current iteratively update 
		/// of the deformation state.
		/// </summary>
		//public bool IsCurrentStateDifferent() => modified;
		//public bool IsCurrentStateDifferent(IHaveState state) => state == null || state.Equals(CurrentState) == false;

		/// <summary>
		/// Resets the boolean that indicates if the constitutive matrix of the material has changed for the current iteratively update 
		/// of the deformation state.
		/// </summary>
		public void ResetModified()
		{
			modified = false;
		}

		/// <summary>
		/// Returns the ID of the material class indicating a specific material law implementation.
		/// </summary>
		public int ID
		{
			get { return 999; }
		}

		/// <summary>
		/// Clears the saved stress strain point connected to the last converged analysis step.
		/// Currently not a valid operation.
		/// </summary>
		public void ClearState() 
		{
			// possibly TODO 
			//if D_tan of initial state is wanted copy calculations from update material for 
			// d_prev_step=0 kai Delta [i] = 0 gia i=0,1 kai 2
			//but don't use it as elastic in other cases
			// maybe in iterative procedure (example provider.Reset ?)
		}

		/// <summary>
		/// Clears tractions - not a valid operation.
		/// </summary>
		public void ClearTractions()
		{
			throw new InvalidOperationException();
		}

		private double youngModulus = 1;
		public double YoungModulus
		{
			get { throw new InvalidOperationException(); } //return youngModulus; }
			set { throw new InvalidOperationException(); }
		}

		private double poissonRatio = 1;
		public double PoissonRatio
		{
			get { throw new InvalidOperationException(); } //return poissonRatio; }
			set { throw new InvalidOperationException(); }
		}

		private double [] coordinates ;
		
		/// <summary>
		/// Returns coordinates. Currently not a valid operation.
		/// </summary>
		public double [] Coordinates
		{

			get { return coordinates; }
			set { throw new InvalidOperationException(); }
		}
	}
}
