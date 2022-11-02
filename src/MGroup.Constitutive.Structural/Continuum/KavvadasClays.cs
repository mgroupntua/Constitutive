#region USING

using System;
using System.Linq;
using System.IO;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Constitutive;
using MGroup.MSolve.DataStructures;
#endregion

namespace MGroup.Constitutive.Structural.Continuum
{

    public class KavvadasClays : IIsotropicContinuumMaterial3D
    {
		#region FieldsConstructor
		// Comments and explanations
		// Based on Kavvadas and Amorosi (2000) Belokas and Kalos investigations. We consider at first the edition of Kavvadas and Amorosi (2000) In the future it will be enriched with next editions. First Edition is on July 2016. Associative plasticity only in the current edition
		//Fields and properties
		private const string PLASTIC_STRAIN = "Plastic strain";
		private const string STRESS_X = "Stress X";
		private const string STRESS_Y = "Stress Y";
		private const string STRESS_Z = "Stress Z";
		private const string STRESS_XY = "Stress XY";
		private const string STRESS_XZ = "Stress XZ";
		private const string STRESS_YZ = "Stress YZ";
		private GenericConstitutiveLawState currentState;
		private double plasticStrainNew;
		public double Zeta { get; set; }
        public double Kmax;
        public double Kmin;
        public bool hasfailed { get; set; }
        public IMatrixView ConstitutiveMatrix { get { return ConstMatr; } set { } }
        private Matrix ConstMatr = Matrix.CreateFromArray(new double[6, 6]);
        private Matrix ElConstMatr = Matrix.CreateFromArray(new double[6, 6]);
        public double[] Coordinates { get; set; }
        public int ID { get; set; }
        public bool Modified { get; set; }
        public double PoissonRatio { get; set; }
        public double criticalstateline { get; set; }
        public double compressibilityfactor { get; set; }
        public double[] Stresses { get; set; }
        public double YoungModulus { get; set; }
        public double ksi { get; set; }
        public double[] Qgrad = new double[6];
        public double[] Pgrad = new double[6];
        public double[] PAR;
        public double[] QH;
        public double[] tempQH;
        public double[] tempStresses;
        public double[] initialStresses = new double[6];
        public readonly double shearModulus;
        public IMatrixView elasticConstitutiveMatrix { get { return ElConstMatr; } set { } } //the readonly was erased due to the change of the elasticconstitutivematric regarding time
        private double[] incrementalStrains = new double[6];
        private double plasticStrain;
        private double[] TotalStrain;
        private double[] stressesNew = new double[6];
        private double Htot;
        private bool modified;
		private double[] stresses;
        public void UpdateMaterial(double[] strainsIncrement)
        {
            incrementalStrains = new double[6];
            for (int j1 = 0; j1 < 6; j1++)
            {
               incrementalStrains[j1]=strainsIncrement[j1];
            }
            //this.incrementalStrains = strainsIncrement.DeepClone();
            for (int i = 0; i < 21; i++)
                this.QH[i] = this.tempQH[i];
            for (int i = 0; i < 6; i++)
                this.Stresses[i] = this.tempStresses[i];
            this.DecideLoad(PAR, incrementalStrains, QH, Stresses, ConstitutiveMatrix);
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

			return stressesNew;
		}
		void ClearState()
        {
            modified = false;
            //this.ConstitutiveMatrix = new Matrix2D<double>(new double[6, 6]);
            //Array.Clear(incrementalStrains,0,5);
            //Array.Clear(Stresses,0,5);
            //Array.Clear(QH, 0, QH.Length);
            //Array.Clear(stressesNew, 0, stressesNew.Length);
            //Array.Clear(plasticStrain,0,5);
        }
        public void SaveState()
        {
            //this.plasticStrain = this.plasticStrainNew;
            //Array.Copy(this.stressesNew, this.Stresses, 6);
            for (int i = 0; i < 21; i++)
                this.tempQH[i] = this.QH[i];
            for (int i = 0; i < 6; i++)
                this.tempStresses[i] = this.Stresses[i];
            this.ConstitutiveMatrix = ConstitutiveMatrix;
            for (int j1 = 0; j1 < 6; j1++)
            {
                TotalStrain[j1] += incrementalStrains[j1];
            }
            if (this.QH[19]+this.QH[20]<0)
            {
                Console.WriteLine("Coordinate X:");
                Console.WriteLine(Coordinates[0]);
                Console.WriteLine("Coordinate Y:");
                Console.WriteLine(Coordinates[1]);
                Console.WriteLine("Coordinate Z:");
                Console.WriteLine(Coordinates[2]);
                Console.WriteLine("Stress Vector:");
                Console.WriteLine(this.Stresses[0]);
                Console.WriteLine(this.Stresses[1]);
                Console.WriteLine(this.Stresses[2]);
                Console.WriteLine(this.Stresses[3]);
                Console.WriteLine(this.Stresses[4]);
                Console.WriteLine(this.Stresses[5]);
                Console.WriteLine("Strains Vector:");
                Console.WriteLine(this.TotalStrain[0]);
                Console.WriteLine(this.TotalStrain[1]);
                Console.WriteLine(this.TotalStrain[2]);
                Console.WriteLine(this.TotalStrain[3]);
                Console.WriteLine(this.TotalStrain[4]);
                Console.WriteLine(this.TotalStrain[5]);
                Console.WriteLine("Plastic volumetric Strain:");
                Console.WriteLine(this.QH[4]);
                Console.WriteLine("Plastic deviatoric Strain:");
                Console.WriteLine(this.QH[5]);
                Console.WriteLine("Plastic Hardening Modulus:");
                Console.WriteLine(this.QH[19] + this.QH[20]);
                this.hasfailed = true;
            }      
        }

        public void ResetModified()
        {
            this.modified = false;
        }


        public void ClearStresses()
        {
            //Array.Clear(Stresses, 0, 5);
            //Array.Clear(stressesNew, 0, 5);
        }

        public object Clone()
        {
            var strainsCopy = new double[incrementalStrains.Length];
            incrementalStrains.CopyTo(strainsCopy, 0);
            var stressesCopy = new double[Stresses.Length];
            Stresses.CopyTo(stressesCopy, 0);
            this.ConstitutiveMatrix = ConstitutiveMatrix;
            //watch out if you use clone.
            KavvadasClays m = new KavvadasClays(this.YoungModulus, this.PoissonRatio, this.shearModulus, this.ksi, this.initialStresses,this.Htot)
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
        public KavvadasClays(double youngModulus, double criticalstateline, double alpha, double ksi)
        {
        }
        public KavvadasClays(double compressibilityfactor, double criticalstateline, double alpha, double ksi, double[] initialStresses,double Htot) : this(compressibilityfactor, criticalstateline, alpha, ksi)
        {
			//this.Coordinates = nodecoordinates;
			    plasticStrain = 0.0;
                TotalStrain = new double[6];
                var gamma = 10; //effective stress
                this.Htot = Htot;
                Zeta = -initialStresses[2] / gamma; // the initialization of the stresses should be in the model class HexaSoil2.cs            
                var Niso = 2.15;//2.08053; //Note that 2*astar is related with the Initial Stresses in X and Y axes (SX=(2*astar-gamma*zeta)/2 for OCR=3)
                this.Stresses = new double[6];
                this.PAR = new double[20];
                this.QH = new double[21];
                this.tempQH = new double[21];
                for (int i = 0; i < 6; i++)
                    Stresses[i] = initialStresses[i];
            if (alpha == 1)
            {
                PAR[3] = compressibilityfactor;
                //Kmax = 1*0.008686;
                //Kmin = youngModulus*0.008686;
            }
            else
            {
                Kmax = 0.5*0.008686;
                Kmin = 0.5*0.008686;
                PAR[3] = (Kmin - Kmax) * Zeta / this.Htot + Kmax;
            }
            PAR[2] = 10 * PAR[3];
                var s0 = (Stresses[0] + Stresses[1] + Stresses[2]) / 3;
                PAR[18] = 4;
                PAR[9] = 4;
            QH[2] = 1.627;
            //QH[6] = 0.5 * Math.Exp((Niso - QH[2] - PAR[3] * Math.Log(Math.Abs(s0))) / (PAR[2] - PAR[3]));
            //QH[6] = Math.Min(0.5 * Math.Exp((Niso - QH[2] - PAR[3] * Math.Log(Math.Abs(s0))) / (PAR[2] - PAR[3])),400);
            QH[6] = 400;
            QH[1] = QH[6] * PAR[9]; //no use of OCR
                PAR[0] = 0;
                PAR[1] = 0.75;
            if (alpha == 1)
            {
                PAR[4] = criticalstateline;
                PAR[5] = criticalstateline;
                PAR[6] = criticalstateline;
                PAR[7] = criticalstateline;
                PAR[8] = criticalstateline;
            }
            else
            {
                PAR[4] = 0.733609251;
                PAR[5] = 0.733609251;
                PAR[6] = 0.733609251;
                PAR[7] = 0.733609251;
                PAR[8] = 0.733609251;
            }
                PAR[10] = 1;
                PAR[11] = 75;
                PAR[12] = 75;
                PAR[13] = 0;
                PAR[14] = 0;
                PAR[15] = ksi;
                PAR[16] = 5;
                PAR[17] = 1;
                PAR[19] = Niso;
                QH[0] = 0;
                QH[3] = 0;
                QH[4] = 0;
                QH[5] = 0;
                QH[7] = 0;
                QH[8] = 0;
                QH[9] = 0;
                QH[10] = 0;
                QH[11] = 0;
                QH[12] = 0;
                initialStresses = TransformToTransformedStresses(initialStresses);
                QH[13] = initialStresses[0];
                QH[14] = initialStresses[1];
                QH[15] = initialStresses[2];
                QH[16] = initialStresses[3];
                QH[17] = initialStresses[4];
                QH[18] = initialStresses[5];
                QH[19] = 0;
                QH[20] = 0;
                initialStresses = TransformToStandardStresses(initialStresses);
                var xk = QH[2] * Math.Abs(s0) / PAR[3];
                var g = PAR[1] * xk / 2;
                var epsilon = (9 * xk * g) / (3 * xk + g);
                var ni = (3 * xk - 2 * g) / (2 * (3 * xk + g));
                var dee = epsilon / ((1 + ni) * (1 - 2 * ni));
                ElConstMatr= Matrix.CreateFromArray(new double[6, 6]);
                ElConstMatr[0, 0] = dee * (1 - ni);
                ElConstMatr[0, 1] = dee * ni;
                ElConstMatr[0, 2] = dee * ni;
                ElConstMatr[1, 0] = dee * ni;
                ElConstMatr[1, 1] = dee * (1 - ni);
                ElConstMatr[1, 2] = dee * ni;
                ElConstMatr[2, 0] = dee * ni;
                ElConstMatr[2, 1] = dee * ni;
                ElConstMatr[2, 2] = dee * (1 - ni);
                ElConstMatr[3, 3] = g;
                ElConstMatr[4, 4] = g;
                ElConstMatr[5, 5] = g;
                this.PoissonRatio = (3 - PAR[1]) / (6 + PAR[1]);
                this.YoungModulus = 3 * (1 - 2 * PoissonRatio) * xk;
                this.elasticConstitutiveMatrix = ElConstMatr;
                this.initialStresses = initialStresses;
                ConstMatr= ElConstMatr;
                stressesNew = new double[6];
                for (int i = 0; i < 21; i++)
                   this.tempQH[i] = this.QH[i];
                this.tempStresses = new double[6] ;
            //Notice that PAR[0] and QH[0] here are not used.
        }

        #endregion
        #region TransformationForStresses
        public double[] TransformToTransformedStresses(double[] Stresses)
        {
            //Transform at first to the Soil Mechanics sign convention
            for (int i=0;i<6;i++)
            {
                Stresses[i] = -Stresses[i];
            }
            var help = new double[6];
            help[0] = (Stresses[0] + Stresses[1] + Stresses[2]) / 3;
            help[1] = (2 * Stresses[1] - Stresses[0] - Stresses[2]) / Math.Sqrt(6);
            help[2] = (Stresses[2] - Stresses[0]) / Math.Sqrt(2);
            help[3] = (Stresses[3]) * Math.Sqrt(2);
            help[4] = (Stresses[4]) * Math.Sqrt(2);
            help[5] = (Stresses[5]) * Math.Sqrt(2);
            Stresses = help;
            return Stresses;
        }
        public double[] TransformToStandardStresses(double[] Stresses)
        {
            var help = new double[6];
            help[0] = Stresses[0] - Stresses[1] / Math.Sqrt(6) - Stresses[2] / Math.Sqrt(2);
            help[1] = Stresses[0] + 2 * Stresses[1] / Math.Sqrt(6);
            help[2] = Stresses[0] - Stresses[1] / Math.Sqrt(6) + Stresses[2] / Math.Sqrt(2);
            help[3] = (Stresses[3]) / Math.Sqrt(2);
            help[4] = (Stresses[4]) / Math.Sqrt(2);
            help[5] = (Stresses[5]) / Math.Sqrt(2);
            Stresses = help;
            //Transform at the end to the engineering convention
            for (int i = 0; i < 6; i++)
            {
                Stresses[i] = -Stresses[i];
            }
            return Stresses;
        }
        #endregion
        #region TransformationForStrains
        public double[] TransformToTransformedStrains(double[] Strains)
        {
            //Transform at first to the Soil Mechanics sign convention
            for (int i = 0; i < 6; i++)
            {
                Strains[i] = -Strains[i];
            }
            var help = new double[6];
            help[0] = (Strains[0] + Strains[1] + Strains[2]);
            help[1] = (2 * Strains[1] - Strains[0] - Strains[2]) / Math.Sqrt(6);
            help[2] = (Strains[2] - Strains[0]) / Math.Sqrt(2);
            help[3] = (Strains[3]) / Math.Sqrt(2);
            help[4] = (Strains[4]) / Math.Sqrt(2);
            help[5] = (Strains[5]) / Math.Sqrt(2);
            Strains = help;
            return Strains;
        }
        public double[] TransformToStandardStrains(double[] Strains)
        {
            var help = new double[6];
            help[0] = Strains[0] / 3 - Strains[1] / Math.Sqrt(6) - Strains[2] / Math.Sqrt(2);
            help[1] = Strains[0] / 3 + 2 * Strains[1] / Math.Sqrt(6);
            help[2] = Strains[0] / 3 - Strains[1] / Math.Sqrt(6) + Strains[2] / Math.Sqrt(2);
            help[3] = (Strains[3]) * Math.Sqrt(2);
            help[4] = (Strains[4]) * Math.Sqrt(2);
            help[5] = (Strains[5]) * Math.Sqrt(2);
            Strains = help;
            //Transform at the end to the engineering convention
            for (int i = 0; i < 6; i++)
            {
                Strains[i] = -Strains[i];
            }
            return Strains;
        }
        #endregion
        #region TransformationForJacobianMatrix 
        //We consider the case of associative plasticity. In other words the plastic potential is the Yield Surface. Therefore the Consistent Constitutive Matrix is symmetric. 
        //The Matrix Mutiplication is done directly from Matlab
        public Matrix TransformationForJacobianMatrix (Matrix ConstitutiveMatrix)
        {
            var m1 = new double[21];
            var s1 = new double[25];
            //Helping assigniments
            m1[0]=ConstitutiveMatrix[0, 0];
            m1[1]=ConstitutiveMatrix[0, 1];
            m1[2]=ConstitutiveMatrix[0, 2];
            m1[3]=ConstitutiveMatrix[0, 3];
            m1[4]=ConstitutiveMatrix[0, 4];
            m1[5]=ConstitutiveMatrix[0, 5];
            m1[6]=ConstitutiveMatrix[1, 1];
            m1[7]=ConstitutiveMatrix[1, 2];
            m1[8]=ConstitutiveMatrix[1, 3];
            m1[9]=ConstitutiveMatrix[1, 4];
            m1[10]=ConstitutiveMatrix[1, 5];
            m1[11]=ConstitutiveMatrix[2, 2];
            m1[12]=ConstitutiveMatrix[2, 3];
            m1[13]=ConstitutiveMatrix[2, 4];
            m1[14]=ConstitutiveMatrix[2, 5];
            m1[15]=ConstitutiveMatrix[3, 3];
            m1[16]=ConstitutiveMatrix[3, 4];
            m1[17]=ConstitutiveMatrix[3, 5];
            m1[18]=ConstitutiveMatrix[4, 4];
            m1[19]=ConstitutiveMatrix[4, 5];
            m1[20]=ConstitutiveMatrix[5, 5];
            s1[0] =m1[1]/Math.Sqrt(6);
            s1[1] =m1[2]/Math.Sqrt(2);
            s1[2] =(m1[2]+Math.Sqrt(6)*m1[7]/3)/Math.Sqrt(2);
            s1[3] = m1[1] * Math.Sqrt(6)/3 ;
            s1[4] = m1[5] / Math.Sqrt(2); 
            s1[5] = m1[4] / Math.Sqrt(2);
            s1[6] = m1[3] / Math.Sqrt(2);
            s1[7] = Math.Sqrt(6) * (m1[1]+m1[6]*Math.Sqrt(6)/3);
            s1[8] = Math.Sqrt(12) * m1[10];
            s1[9] = Math.Sqrt(12) * m1[9];
            s1[10] =Math.Sqrt(12) * m1[8];
            s1[11] =m1[10]/ Math.Sqrt(6);
            s1[12] =m1[14]/ Math.Sqrt(2);
            s1[13] =m1[9]/ Math.Sqrt(6);
            s1[14] =m1[13]/ Math.Sqrt(2);
            s1[15] =m1[8]/ Math.Sqrt(6);
            s1[16] =m1[12]/ Math.Sqrt(2);
            s1[21] = m1[7]/ Math.Sqrt(6);
            s1[22] =m1[11]/ Math.Sqrt(2);
            s1[23] = m1[6] / Math.Sqrt(6);
            s1[24] = m1[7] / Math.Sqrt(2);
            s1[17] = (s1[21]+s1[22]-m1[2]) / Math.Sqrt(2);
            s1[18] = (-s1[21] + s1[22] + m1[2]) / Math.Sqrt(2);
            s1[19] = Math.Sqrt(6) * (s1[24]+s1[23]-m1[1]);
            s1[20] = Math.Sqrt(6) * (s1[24] - s1[23] + m1[1]);
            //Transform
            ConstitutiveMatrix[0, 0] = m1[0]-s1[0]-s1[1]+s1[19]/6+s1[17];
            ConstitutiveMatrix[0, 1] = m1[0] - s1[0] - s1[1] - s1[19] / 3 ;
            ConstitutiveMatrix[0, 2] = m1[0] - s1[0] - s1[1] + s1[19] / 6 - s1[17];
            ConstitutiveMatrix[0, 3] = -(s1[15]+s1[16]-m1[3])/Math.Sqrt(2);
            ConstitutiveMatrix[0, 4] = -(s1[13] + s1[14] - m1[4]) / Math.Sqrt(2);
            ConstitutiveMatrix[0, 5] = -(s1[11] + s1[12] - m1[5]) / Math.Sqrt(2);
            ConstitutiveMatrix[1, 1] = m1[0] + s1[3] + s1[7] / 3;
            ConstitutiveMatrix[1, 2] = m1[0] + s1[3] - s1[7] / 6 + s1[2];
            ConstitutiveMatrix[1, 3] = s1[6] + s1[10] / 6;
            ConstitutiveMatrix[1, 4] = s1[5] + s1[9] / 6;
            ConstitutiveMatrix[1, 5] = s1[4] + s1[8] / 6;
            ConstitutiveMatrix[2, 2] = m1[0] + s1[1] - s1[0] - s1[20] / 6 + s1[18];
            ConstitutiveMatrix[2, 3] = (s1[16] - s1[15] + m1[3]) / Math.Sqrt(2);
            ConstitutiveMatrix[2, 4] = (s1[14] - s1[13] + m1[4]) / Math.Sqrt(2);
            ConstitutiveMatrix[2, 5] = (s1[12] - s1[11] + m1[5]) / Math.Sqrt(2);
            ConstitutiveMatrix[3, 3] = m1[15]/2;
            ConstitutiveMatrix[3, 4] = m1[16]/2;
            ConstitutiveMatrix[3, 5] = m1[17]/2;
            ConstitutiveMatrix[4, 4] = m1[18]/2;
            ConstitutiveMatrix[4, 5] = m1[19]/2;
            ConstitutiveMatrix[5, 5] = m1[20]/2;
            //Apply Symmetry      
            ConstitutiveMatrix[1, 0] = ConstitutiveMatrix[0, 1];
            ConstitutiveMatrix[2, 0] = ConstitutiveMatrix[0, 2];
            ConstitutiveMatrix[3, 0] = ConstitutiveMatrix[0, 3];
            ConstitutiveMatrix[4, 0] = ConstitutiveMatrix[0, 4];
            ConstitutiveMatrix[5, 0] = ConstitutiveMatrix[0, 5];
            ConstitutiveMatrix[2, 1] = ConstitutiveMatrix[1, 2];
            ConstitutiveMatrix[3, 1] = ConstitutiveMatrix[1, 3];
            ConstitutiveMatrix[4, 1] = ConstitutiveMatrix[1, 4];
            ConstitutiveMatrix[5, 1] = ConstitutiveMatrix[1, 5];
            ConstitutiveMatrix[3, 2] = ConstitutiveMatrix[2, 3];
            ConstitutiveMatrix[4, 2] = ConstitutiveMatrix[2, 4];
            ConstitutiveMatrix[5, 2] = ConstitutiveMatrix[2, 5];
            ConstitutiveMatrix[4, 3] = ConstitutiveMatrix[3, 4];
            ConstitutiveMatrix[5, 3] = ConstitutiveMatrix[3, 5];
            ConstitutiveMatrix[5, 4] = ConstitutiveMatrix[4, 5];
            this.ConstitutiveMatrix = ConstitutiveMatrix;
            return ConstitutiveMatrix;
        }
        #endregion
        # region LoadSubroutines
        public void DecideLoad(double[] PAR, double[] de, double[] QH, double[] Stresses, IMatrixView ConstitutiveMatrix)
        {
            for (int i = 0; i < 6; i++)
                Stresses[i] = Stresses[i] + initialStresses[i];
            Stresses = TransformToTransformedStresses(Stresses);
            de = TransformToTransformedStrains(de);
            var d1 = new double[6];
            de.CopyTo(d1, 0);
            var maxde = d1.Max();
            var minde = d1.Min();
            maxde = Math.Abs(maxde);
            minde = Math.Abs(minde);
            maxde = Math.Max(maxde, minde);
            if (maxde > Math.Pow(10, -5))
            {
                this.ConstitutiveMatrix= Matrix.CreateFromArray(new double[6, 6]);
                this.elasticConstitutiveMatrix= Matrix.CreateFromArray(new double[6, 6]);
                this.ConstMatr= Matrix.CreateFromArray(new double[6, 6]);
                this.ElConstMatr= Matrix.CreateFromArray(new double[6, 6]);
                var ndiv = (double)Math.Truncate(maxde / (Math.Pow(10, -5))) + 1;
                for (int i = 0; i < 6; i++)
                {
                    de[i] = de[i] / ndiv;
                }
                for (int i = 0; i < (int)Math.Truncate(ndiv); i++)
                {
                    if (QH[3] == 0)
                    {
                        Stresses = unload(PAR, de, QH, Stresses);
                    }
                    else
                    {
                        var help = 0.0;
                        help = Numer(PAR, de, QH, Stresses);
                        if (help < 0)
                        {
                            Stresses = unload(PAR, de, QH, Stresses);
                        }
                        else if (help > 0)
                        {
                            Stresses = load(PAR, de, QH, Stresses);
                        }
                        else if (help == 0)
                        {
                            QH = loadneutral(PAR, de, QH, Stresses);
                        }
                    }
                    var helpc = cstif(PAR, QH, Stresses);
                    ConstMatr[0, 0] += helpc[0, 0] / ndiv;
                    ConstMatr[0, 1] += helpc[0, 1] / ndiv;
                    ConstMatr[0, 2] += helpc[0, 2] / ndiv;
                    ConstMatr[0, 3] += helpc[0, 3] / ndiv;
                    ConstMatr[0, 4] += helpc[0, 4] / ndiv;
                    ConstMatr[0, 5] += helpc[0, 5] / ndiv;
                    ConstMatr[1, 1] += helpc[1, 1] / ndiv;
                    ConstMatr[1, 2] += helpc[1, 2] / ndiv;
                    ConstMatr[1, 3] += helpc[1, 3] / ndiv;
                    ConstMatr[1, 4] += helpc[1, 4] / ndiv;
                    ConstMatr[1, 5] += helpc[1, 5] / ndiv;
                    ConstMatr[2, 2] += helpc[2, 2] / ndiv;
                    ConstMatr[2, 3] += helpc[2, 3] / ndiv;
                    ConstMatr[2, 4] += helpc[2, 4] / ndiv;
                    ConstMatr[2, 5] += helpc[2, 5] / ndiv;
                    ConstMatr[3, 3] += helpc[3, 3] / ndiv;
                    ConstMatr[3, 4] += helpc[3, 4] / ndiv;
                    ConstMatr[3, 5] += helpc[3, 5] / ndiv;
                    ConstMatr[4, 4] += helpc[4, 4] / ndiv;
                    ConstMatr[4, 5] += helpc[4, 5] / ndiv;
                    ConstMatr[5, 5] += helpc[5, 5] / ndiv;
                    //Apply Symmetry      
                    ConstMatr[1, 0] = ConstMatr[0, 1];
                    ConstMatr[2, 0] = ConstMatr[0, 2];
                    ConstMatr[3, 0] = ConstMatr[0, 3];
                    ConstMatr[4, 0] = ConstMatr[0, 4];
                    ConstMatr[5, 0] = ConstMatr[0, 5];
                    ConstMatr[2, 1] = ConstMatr[1, 2];
                    ConstMatr[3, 1] = ConstMatr[1, 3];
                    ConstMatr[4, 1] = ConstMatr[1, 4];
                    ConstMatr[5, 1] = ConstMatr[1, 5];
                    ConstMatr[3, 2] = ConstMatr[2, 3];
                    ConstMatr[4, 2] = ConstMatr[2, 4];
                    ConstMatr[5, 2] = ConstMatr[2, 5];
                    ConstMatr[4, 3] = ConstMatr[3, 4];
                    ConstMatr[5, 3] = ConstMatr[3, 5];
                    ConstMatr[5, 4] = ConstMatr[4, 5];
                    var xk = QH[2] * Stresses[0] / PAR[3];
                    var g = PAR[1] * xk / 2;
                    var epsilon = (9 * xk * g) / (3 * xk + g);
                    var ni = (3 * xk - 2 * g) / (2 * (3 * xk + g));
                    var dee = epsilon / ((1 + ni) * (1 - 2 * ni));
                    ElConstMatr[0, 0] += dee * (1 - ni) / ndiv;
                    ElConstMatr[0, 1] += dee * ni / ndiv;
                    ElConstMatr[0, 2] += dee * ni / ndiv;
                    ElConstMatr[1, 0] += dee * ni / ndiv;
                    ElConstMatr[1, 1] += dee * (1 - ni) / ndiv;
                    ElConstMatr[1, 2] += dee * ni / ndiv;
                    ElConstMatr[2, 0] += dee * ni / ndiv;
                    ElConstMatr[2, 1] += dee * ni / ndiv;
                    ElConstMatr[2, 2] += dee * (1 - ni) / ndiv;
                    ElConstMatr[3, 3] += g / ndiv;
                    ElConstMatr[4, 4] += g / ndiv;
                    ElConstMatr[5, 5] += g / ndiv;
                    this.PoissonRatio = (3 - PAR[1]) / (6 + PAR[1]);
                    this.YoungModulus = 3 * (1 - 2 * PoissonRatio) * xk;
                }
                ConstMatr = TransformationForJacobianMatrix(ConstMatr);
                this.ConstitutiveMatrix = ConstMatr;
                this.elasticConstitutiveMatrix = ElConstMatr;
            }
            else
            {
                if (QH[3] == 0)
                {
                    Stresses = unload(PAR, de, QH, Stresses);
                }
                else
                {
                    var help = 0.0;
                    help = Numer(PAR, de, QH, Stresses);
                    if (help < 0)
                    {
                        Stresses = unload(PAR, de, QH, Stresses);
                    }
                    else if (help > 0)
                    {
                        Stresses = load(PAR, de, QH, Stresses);
                    }
                    else if (help == 0)
                    {
                        QH = loadneutral(PAR, de, QH, Stresses);
                    }
                }
                this.ConstitutiveMatrix= Matrix.CreateFromArray(new double[6, 6]);
                this.ConstMatr= Matrix.CreateFromArray(new double[6, 6]);
                ConstMatr = cstif(PAR, QH, Stresses);
                ConstMatr = TransformationForJacobianMatrix(ConstMatr);
                this.ConstitutiveMatrix = ConstMatr;
                var xk = QH[2] * Stresses[0] / PAR[3];
                var g = PAR[1] * xk / 2;
                var epsilon = (9 * xk * g) / (3 * xk + g);
                var ni = (3 * xk - 2 * g) / (2 * (3 * xk + g));
                var dee = epsilon / ((1 + ni) * (1 - 2 * ni));
                this.elasticConstitutiveMatrix= Matrix.CreateFromArray(new double[6, 6]);
                ElConstMatr= Matrix.CreateFromArray(new double[6, 6]);
                ElConstMatr[0, 0] = dee * (1 - ni);
                ElConstMatr[0, 1] = dee * ni;
                ElConstMatr[0, 2] = dee * ni;
                ElConstMatr[1, 0] = dee * ni;
                ElConstMatr[1, 1] = dee * (1 - ni);
                ElConstMatr[1, 2] = dee * ni;
                ElConstMatr[2, 0] = dee * ni;
                ElConstMatr[2, 1] = dee * ni;
                ElConstMatr[2, 2] = dee * (1 - ni);
                ElConstMatr[3, 3] = g;
                ElConstMatr[4, 4] = g;
                ElConstMatr[5, 5] = g;
                this.elasticConstitutiveMatrix = ElConstMatr;
                this.PoissonRatio = (3 - PAR[1]) / (6 + PAR[1]);
                this.YoungModulus = 3 * (1 - 2 * PoissonRatio) * xk;
            }
            Stresses = TransformToStandardStresses(Stresses);
            de = TransformToStandardStrains(de);
            for (int i = 0; i < 6; i++)
                Stresses[i] = Stresses[i] - initialStresses[i];
            for (int i = 0; i < 6; i++)
                this.Stresses[i] = Stresses[i];
            for (int i = 0; i < 21; i++)
                this.QH[i] = QH[i];
            this.modified = Modified = true;
            //Console.WriteLine(this.QH[2]);
        }
        public double[] load(double[] PAR,double[] de, double[] QH, double[] Stresses)
        {
            if (QH[3] == -1)
            {
                Stresses = yload(PAR, de, QH, Stresses);
            }
            else
            {
                Stresses = bload(PAR, de, QH, Stresses);
            }
            return Stresses;
        }
        public double[] bload( double[] PAR,double[] de, double[] QH,double[] Stresses)
        {
            var Dpp = 0.0;
            var XLDOT = 0.0;
            var Devp = 0.0;
            var Deqp = 0.0;
            var Devptrue = 0.0;
            var Compevp = 0.0;
            var Compeqp = 0.0;
            var Dexpart = 0.0;
            var Dpart1 = 0.0;
            var DD = new double[6];
            var ds = new double[6];
            var Alpha = 0.0;
            XLDOT = compldot(PAR, de, QH, Stresses);
            Qgrad = compp(PAR, QH, Stresses);
            Dpp = compdotpdot(Qgrad);
            Devp = Math.Abs(XLDOT * Qgrad[0]);
            Deqp = Math.Abs(XLDOT * Dpp);
            Devptrue = Math.Abs(XLDOT) * Qgrad[0];
            Compevp = QH[4] + Devp;
            Compeqp = QH[5] + Deqp;
            Dexpart = Math.Exp(-Math.Abs(PAR[11] * Compevp + PAR[12] * Compeqp));
            Dpart1 = Math.Abs(PAR[9] - PAR[10] + PAR[13] * Compevp + PAR[14] * Compeqp);
            QH[4] = Compevp;
            QH[5] = Compeqp;
            QH[6] = QH[6] * (1 + QH[2] * Devptrue / (PAR[2] - PAR[3]));
            QH[1] = QH[6] * (Dpart1 * Dexpart + PAR[10]);
            for (int i = 0; i < 6; i++)
            {
                DD[i] = de[i] - Math.Abs(XLDOT) * Qgrad[i];
            }
            ds = compds(PAR, DD, QH, Stresses);
            for (int i = 0; i < 6; i++)
            {
                Stresses[i] = Stresses[i] + ds[i];
            }
            Alpha = QH[1];
            QH[2] = QH[2] * (1 - de[0]);
            Stresses=adjbse(PAR, Alpha, Stresses);
            compsl(PAR, QH, Stresses);
            QH[3] = 1;
            this.QH = QH;
            return Stresses;
        }
        public double[] yload( double[] PAR,double[] de,  double[] QH, double[] Stresses)
        {
            var Dpp = 0.0;
            var XLDOT = 0.0;
            var Devp = 0.0;
            var Deqp = 0.0;
            var Devptrue = 0.0;
            var Compevp = 0.0;
            var Compeqp = 0.0;
            var Dexpart = 0.0;
            var Dpart1 = 0.0;
            var QH6 = 0.0;
            var Alpha = 0.0;
            var DA = 0.0;
            var Ratio = 0.0;
            var YPLDOT = 0.0;
            var YLDOT = 0.0;
            var Fcheck = 0.0;
            var RR = 0.0;
            var B = new double[6];
            var ds = new double[6];
            var DD = new double[6];
            var DD1 = new double[6];
            var Stressescheck = new double[6];
            XLDOT = compldot(PAR, de, QH, Stresses);
            Qgrad = compp(PAR, QH, Stresses);
            Dpp = compdotpdot(Qgrad);
            for (int i = 0; i < 6; i++)
            {
                DD[i] = de[i] - Math.Abs(XLDOT) * Qgrad[i];
            }
            ds = compds(PAR, DD, QH, Stresses);
            B[0] = (Stresses[0] - QH[13]) / PAR[15] - (Stresses[0] - QH[1]);
            for (int i = 1; i < 6; i++)
            {
                B[i] = (Stresses[i] - QH[i + 13]) / PAR[15] - Stresses[i];
            }
            Devp = Math.Abs(XLDOT * Qgrad[0]);
            Deqp = Math.Abs(XLDOT * Dpp);
            Devptrue = Math.Abs(XLDOT) * Qgrad[0];
            Compevp = QH[4] + Devp;
            Compeqp = QH[5] + Deqp;
            Dexpart = Math.Exp(-Math.Abs(PAR[11] * Compevp + PAR[12] * Compeqp));
            Dpart1 = Math.Abs(PAR[9] - PAR[10] + PAR[13] * Compevp + PAR[14] * Compeqp);
            QH6 = QH[6] * (1 + QH[2] * Devptrue / (PAR[2] - PAR[3]));
            Alpha = QH6 * (Dpart1 * Dexpart + PAR[10]);
            DA = Alpha - QH[1];
            Ratio = DA / QH[1];
            YPLDOT = 1 + Ratio;
            for (int i = 0; i < 6; i++)
            {
                Stressescheck[i] = Stresses[i] + ds[i];
            }
            Fcheck = FBSE(PAR, Alpha, Stressescheck);
            YLDOT = compyldot(PAR, QH, Stresses, ds, de, DA);
            YPLDOT = 1 + Ratio;
            if (Fcheck<0)
            {
                QH[1] = Alpha;
                for (int i = 0; i < 6; i++)
                {
                    QH[13 + i] = YPLDOT * QH[13 + i] + YLDOT * B[i];          
                }
                for (int i = 0; i < 6; i++)
                {
                    Stresses[i] = Stresses[i] + ds[i];
                }
                QH[4] = Compevp;
                QH[5] = Compeqp;
                QH[6] = QH6;
                QH[2] = QH[2] * (1 - de[0]);
                Stresses=adjpye(PAR,QH,Stresses);
                QH[3] = -1;
            }
            else
            {
                RR = IntBse(PAR, QH, Stresses, ds);
                for (int i = 0; i < 6; i++)
                {
                    DD[i] = (1- RR)*de[i];
                    DD1[i] = (RR) * de[i];
                }
                XLDOT = compldot(PAR, DD1, QH, Stresses);
                Qgrad = compp(PAR, QH, Stresses);
                Dpp = compdotpdot(Qgrad);
                Devp = Math.Abs(XLDOT * Qgrad[0]);
                Deqp = Math.Abs(XLDOT * Dpp);
                Devptrue = Math.Abs(XLDOT) * Qgrad[0];
                Compevp = QH[4] + Devp;
                Compeqp = QH[5] + Deqp;
                Dexpart = Math.Exp(-Math.Abs(PAR[11] * Compevp + PAR[12] * Compeqp));
                Dpart1 = Math.Abs(PAR[9] - PAR[10] + PAR[13] * Compevp + PAR[14] * Compeqp);
                QH[6] = QH[6] * (1 + QH[2] * Devptrue / (PAR[2] - PAR[3]));
                QH[1] = QH[6] * (Dpart1*Dexpart+PAR[10]);
                QH[4] = Compevp;
                QH[5] = Compeqp;
                QH[2] = QH[2] * (1 - DD1[0]);
                for (int i = 0; i < 6; i++)
                {
                    Stresses[i] = Stresses[i] + RR*ds[i];
                }
                Stresses=adjbse(PAR, QH[1], Stresses);
                compsl(PAR, QH, Stresses);
                QH[3] = 1;
                if (RR<1)
                {
                    Stresses=bload(PAR,DD, QH, Stresses);
                }
            }
            this.QH = QH;
            return Stresses;
        }
        public double[] unload( double[] PAR, double[] de,  double[] QH,  double[] Stresses)
        {
            var ds = new double[6];
            var DV = 0.0;
            var Stressescheck = new double[6];
            var Alpha = 0.0;
            var Fcheck = 0.0;
            var Fcheck2 = 0.0;
            var XLamda = 0.0;
            var DD1 = new double[6];
            ds=compds (PAR, de, QH, Stresses);
            DV = -de[0] * QH[2];
            for (int i = 0; i < 6; i++)
            {
                Stressescheck[i] = Stresses[i] + ds[i];
            }
            Alpha = QH[1];
            Fcheck = FPYE(PAR, QH, Stressescheck);
            if (Fcheck<0)
            {
                Stresses = Stressescheck;
                QH[2] = QH[2] + DV;
                QH[3] = 0;
            }
            else
            {
                XLamda=IntPye(PAR, QH, Stresses, ds);
                for (int i = 0; i < 6; i++)
                {
                    Stresses[i] = Stresses[i] + XLamda*ds[i];
                    DD1[i] = (1 - XLamda) * de[i];
                }
                QH[2] = QH[2] + XLamda * DV;
                Stresses=adjpye(PAR, QH, Stresses);
                Fcheck2 = FBSE(PAR,Alpha,Stresses);
                if (Fcheck2<0)
                {
                    QH[3] = -1;
                }
                else
                {
                    Stresses=adjbse(PAR, Alpha, Stresses);
                    compsl(PAR, QH, Stresses);
                    QH[3] = 1;
                }
                for (int i = 0; i < 6; i++)
                {
                    QH[7 + i] = Stresses[i];
                }
                Stresses=load(PAR, DD1, QH, Stresses);
            }
            this.QH = QH;
            return Stresses;
        }
        public double [] loadneutral(double[] PAR,double[] de, double[] QH, double[] Stresses)
        {
            QH[2] = QH[2] * (1 - de[0]);
            return QH;
        }
        #endregion
        #region Numer/FB/FP/INT/ADJSubroutines
        public double Numer (double[] PAR, double[] de,double [] QH,double[] Stresses)
        {
            var Xnum = 0.0;
            Qgrad = compp(PAR, QH, Stresses);
            Xnum = Quad(PAR,de,QH,Stresses,Qgrad);
            return Xnum;
        }
        public double FBSE(double[] PAR,double Alpha, double[] Stresses)
        {
            var F = 0.0;
            F = Math.Pow((Alpha - Stresses[0]), 2) - Math.Pow(Alpha, 2);
            for (int i = 1; i < 6; i++)
            {
                F = F + Math.Pow((Stresses[i] / PAR[i + 3]), 2);
            }
            return F;
        }
        public double FPYE(double[] PAR, double [] QH, double[] Stresses)
        {
            var F = 0.0;
            F = Math.Pow((QH[13] - Stresses[0]), 2) - Math.Pow((PAR[15]*QH[1]), 2);
            for (int i = 1; i < 6; i++)
            {
                F = F + Math.Pow(((Stresses[i]-QH[13+i]) / PAR[i + 3]), 2);
            }
            return F;
        }
        public double[] adjbse(double [] PAR,double Alpha, double[] Stresses)
        {
            var F1 = 0.0;
            var Alpha1 = 0.0;
            var XLamda = 0.0;
            F1 = FBSE(PAR, Alpha, Stresses);
            Alpha1 = Math.Pow(Alpha, 2);
            XLamda = Math.Sqrt(Math.Abs(Alpha1 / (F1 + Alpha1)));
            Stresses[0] = Alpha + XLamda * (Stresses[0] - Alpha);
            for (int i = 1; i < 6; i++)
            {
                Stresses[i] = XLamda * Stresses[i];
            }
            return Stresses;
        }
        public double[] adjpye(double [] PAR,double [] QH, double[] Stresses)
        {
            var F1 = 0.0;
            var Alpha1 = 0.0;
            var XLamda = 0.0;
            F1 = FPYE(PAR, QH, Stresses);
            Alpha1 = Math.Pow((QH[1] * PAR[15]), 2);
            XLamda = Math.Sqrt(Math.Abs(Alpha1 / (F1 + Alpha1)));
            for (int i = 0; i < 6; i++)
            {
                Stresses[i] = XLamda * Stresses[i] + (1 - XLamda) * QH[13 + i];
            }
            return Stresses;
        }
        public double IntBse(double [] PAR,double [] QH, double[] Stresses,double[] Ds)
        {
            var AA = 0.0;
            var BB = 0.0;
            var CC = 0.0;
            var dummy = 0.0;
            var XLamda = 0.0;
            AA = Math.Abs(Math.Pow(Ds[0], 2));
            BB = Ds[0] * (Stresses[0] - QH[1]);
            for (int i = 1; i < 6; i++)
            {
                BB = BB + (Ds[i] * Stresses[i]) / Math.Pow(PAR[3 + i], 2);
                AA = AA + Math.Pow((Ds[i] / PAR[3 + i]), 2);
            }
            if (AA == 0)
            {
                XLamda = 1;
            }
            CC = FBSE(PAR, QH[1], Stresses);
            dummy = Math.Pow(BB, 2) - AA * CC;
            XLamda = (-BB + Math.Sqrt(Math.Abs(dummy))) / (AA);
            if (CC == 0)
            {
                XLamda = 1;
            }
            if (XLamda > 1)
            {
                XLamda = 1;
            }
            if (XLamda < 0)
            {
                XLamda = 0;
            }
            return XLamda;
        }
        public double IntPye(double [] PAR,double[] QH, double[] Stresses,double[] Ds)
        {
            var AA = 0.0;
            var BB = 0.0;
            var CC = 0.0;
            var dummy = 0.0;
            var XLamda = 0.0;
            AA = Math.Abs(Math.Pow(Ds[0], 2));
            BB = Ds[0] * (Stresses[0] - QH[13]);
            for (int i = 1; i < 6; i++)
            {
                BB = BB + (Ds[i] * (Stresses[i] - QH[13 + i])) / Math.Pow(PAR[3 + i], 2);
                AA = AA + Math.Pow((Ds[i] / PAR[3 + i]), 2);
            }
            if (AA == 0)
            {
                XLamda = 1;
            }
            CC = FPYE(PAR, QH, Stresses);
            dummy = Math.Pow(BB, 2) - AA * CC;
            XLamda = (-BB + Math.Sqrt(Math.Abs(dummy))) / (AA);
            if (CC == 0)
            {
                XLamda = 1;
            }
            if (XLamda > 1)
            {
                XLamda = 1;
            }
            if (XLamda < 0)
            {
                XLamda = 0;
            }
            return XLamda;
        }
        #endregion
        #region CompSubroutines
        public void compsl (double [] PAR,double [] QH,double[] Stresses)
        {
            QH[13] = (1 - PAR[15]) * Stresses[0] + PAR[15] * QH[1];
            for (int i=1; i<6; i++)
            {
                QH[13 + i] = (1 - PAR[15]) * Stresses[i];
            }
            this.QH = QH;
        }
        public double compldot (double [] PAR,double[] de,double [] QH,double[] Stresses)
        {
            var Xn = 0.0;
            var Xd = 0.0;
            var H = 0.0;
            var Xldot = 0.0;
            Qgrad = compp(PAR, QH, Stresses);
            Xn = Quad(PAR, de, QH, Stresses, Qgrad);
            var qg = new double[6];
            for (int j1=0;j1<6;j1++)
            {
                qg[j1] = Qgrad[j1];
            }
            Xd = Quad(PAR, qg, QH, Stresses, Qgrad);
            H = comph(PAR,QH,Stresses);
            Xd = Xd + H;
            if (Math.Abs(Xn) < Math.Pow(10, -15))
            {
                Xldot = 0.0;
            }
            else
            {
                if (Xd == 0)
                {
                    Xd = Math.Pow(10, -15);
                }
                Xldot = Xn / Xd;
            }
            return Xldot;
        }
        public double compyldot(double [] PAR,double [] QH,double[] Stresses,double[] ds,double[] de,double DA)
        {
            var yldot = 0.0;
            var Xn = 0.0;
            var Xd = 0.0;
            var Ratio = 0.0;
            Ratio = DA / QH[1];
            for (int i = 1; i < 6; i++)
                {
                  Xn = Xn + (Stresses[i]-QH[13+i]) * (ds[i]-Ratio*Stresses[i]) / (Math.Pow(PAR[3+i],2));
                  Xd = Xd + (Stresses[i] - QH[13 + i]) * (Stresses[i]) / (Math.Pow(PAR[3 + i], 2));
                }
            Xn=Xn+(Stresses[0]-QH[13]) * (ds[0]-Ratio* Stresses[0]);
            Xd=Xd+(Stresses[0] - QH[13]) * (Stresses[0]-QH[1]);
            Xd = PAR[15] * Math.Pow(QH[1], 2) - Xd;
            if (Xd == 0)
            {
                Xd = Math.Pow(10, -15);
            }
            yldot = Xn / Xd;
            if (Math.Abs(Xn) < Math.Pow(10, -15))
            {
                yldot = 0.0;
            }
            return yldot;
        }
        public double[] compp(double[] PAR, double[] QH, double[] Stresses)
        {
            Qgrad[0] = 2*(Stresses[0] - QH[13]);
            for (int i = 1; i < 6; i++)
            {
                Qgrad[i] = 2 * (Stresses[i] - QH[13 + i]) / Math.Pow(PAR[3+i],2);
            }
            return Qgrad;
        }
        public double[] comppbse(double[] PAR, double[] QH,double[] Stresses)
        {
            Pgrad[0] = 2 * (Stresses[0] - QH[1])*PAR[15];
            for (int i = 1; i < 6; i++)
            {
                Pgrad[i] = 2 *PAR[15]* (Stresses[i]) / Math.Pow(PAR[3 + i], 2);
            }
            return Pgrad;
        }
        public double compdotpdot (double [] Pgrad)
        {
            var Dpp = 0.0;
            for (int i=1; i<6; i++)
            {
                Dpp = Dpp + Pgrad[i]*Pgrad[i];
            }
            Dpp = Dpp / 1.5;
            Dpp = Math.Sqrt(Dpp);
            return Dpp;
        }
        public double[] conjugate(double[] PAR,double[] QH,double[] Stresses)
        {
            var SM = new double[6];
            SM[0] =QH[1]+(Stresses[0]-QH[13])/PAR[15];
            for (int i=1;i<6;i++)
            {
                SM[i] =(Stresses[i]-QH[13+i])/PAR[15];
            }
            return SM;
        }
        public double comphfact (double [] PAR,double [] QH,double[] Stresses,double[] SM)
        {
            var AA = 0.0;
            var BB = 0.0;
            var Halfterm1 = 0.0;
            var Halfterm2 = 0.0;
            var Hfact = 0.0;
            var Xd = 0.0;
            for (int i = 1; i<6;i++)
            {
                Halfterm1 =(SM[i]-Stresses[i])/PAR[3+i];
                Halfterm2 = (Stresses[i] - QH[7+i]) / PAR[3 + i];
                AA = AA + Halfterm1 * Halfterm1;
                BB = BB + Halfterm2 * Halfterm2;
            }
            AA = AA + Math.Pow((SM[0]-Stresses[0]), 2);
            BB = BB + Math.Pow((Stresses[0]-QH[7]), 2);
            if (Math.Abs(AA)<Math.Pow(10,-15))
            {
                Hfact = 0.0;
            }
            else
            {
                Qgrad = compp(PAR, QH, Stresses);
                var qg = new double[6];
                for (int j1 = 0; j1 < 6; j1++)
                {
                    qg[j1] = Qgrad[j1];
                }
                Xd = Quad(PAR, qg, QH, Stresses, Qgrad);
                if (Math.Abs(BB)==0)
                {
                    BB = Math.Pow(10, -15);
                }
                Hfact = PAR[16] * Xd * Math.Pow((AA / BB), PAR[17]);
            }
            return Hfact;
        }
        public double comph (double [] PAR,double [] QH,double[] Stresses)
        {
            var H = 0.0;
            var Hf1 = 0.0;
            var Hf2 = 0.0;
            var Dpart1 = 0.0;
            var Dpart2 = 0.0;
            var Dpart3 = 0.0;
            var Dexpart = 0.0;
            var Hpart1 = 0.0;
            var Hpart2 = 0.0;
            var Hpart3 = 0.0;
            var Dpp = 0.0;
            var SM = new double[6];
            var Hfact = 0.0;
            if (QH[3] == 1)
            {
                Qgrad = compp(PAR, QH, Stresses);
                Dpp = compdotpdot(Qgrad);
                Dexpart = Math.Exp(-Math.Abs(PAR[11] * QH[4] + PAR[12] * QH[5]));
                Dpart1 = Math.Abs(PAR[9]-PAR[10]+PAR[13] * QH[4] + PAR[14] * QH[5]);
                Dpart2 = Math.Abs(PAR[11] * Math.Abs(Qgrad[0]) + PAR[12] * Dpp);
                Dpart3 = Math.Abs(PAR[13] * Math.Abs(Qgrad[0]) + PAR[14] * Dpp);
                Hpart1 =(QH[2])*(Qgrad[0])*(Dpart1*Dexpart+PAR[10])/(PAR[2]-PAR[3]);
                Hpart2 = Math.Abs(Dpart3*Dexpart);
                Hpart3 =-Math.Abs(Dpart1*Dpart2*Dexpart);
                H =2*Stresses[0]*QH[6]*PAR[15]*(Hpart1+Hpart2+Hpart3);
                return H;
            }
            else
            {
                SM = conjugate(PAR, QH, Stresses);
                Hfact = comphfact(PAR, QH, Stresses, SM);
                Pgrad = comppbse(PAR, QH, SM);
                Dpp = compdotpdot(Pgrad);
                Dexpart = Math.Exp(-Math.Abs(PAR[11] * QH[4] + PAR[12] * QH[5]));
                Dpart1 = Math.Abs(PAR[9] - PAR[10] + PAR[13] * QH[4] + PAR[14] * QH[5]);
                Dpart2 = Math.Abs(PAR[11] * Math.Abs(Pgrad[0]) + PAR[12] * Dpp);
                Dpart3 = Math.Abs(PAR[13] * Math.Abs(Pgrad[0]) + PAR[14] * Dpp);
                Hpart1 = (QH[2]) * (Pgrad[0]) * (Dpart1 * Dexpart + PAR[10]) / (PAR[2] - PAR[3]);
                Hpart2 = Math.Abs(Dpart3 * Dexpart);
                Hpart3 = -Math.Abs(Dpart1 * Dpart2 * Dexpart);
                Hf1 =2*SM[0]* QH[6] * PAR[15] * (Hpart1 + Hpart2 + Hpart3);
                Hf2 =Hfact;
                QH[19] = Hf1;
                QH[20] = Hf2;
                H = Hf1 + Hf2;
                this.QH = QH;
                return H;
            }
        }
        #endregion
        #region Cstif/Compds/QuadSubroutines
        public Matrix cstif (double [] PAR,double [] QH,double[] Stresses)
        {
            var Xk = 0.0;
            var shear = 0.0;
            var dd = new double[6];
            var ee = new double[6];
            var omega = 0.0;
            var H = 0.0;
            Matrix paramConstMatr = Matrix.CreateFromArray(new double[6, 6]);
            Xk = QH[2] * Stresses[0] / PAR[3];
            shear = Xk * PAR[1];
            paramConstMatr[0, 0] = Xk;
            for (int i = 1; i<6;i++)
            {
                paramConstMatr[i, i] = shear;
            }
            if (QH[3] ==0)
            {
                return paramConstMatr;
            }
            else
            {
                Qgrad=compp(PAR, QH, Stresses);
                var qg = new double[6];
                for (int j1 = 0; j1 < 6; j1++)
                {
                    qg[j1] = Qgrad[j1];
                }
                dd = compds(PAR, qg, QH, Stresses);
                for (int i = 0; i< 6;i++)
                {
                    omega = omega +Qgrad[i]*dd[i];
                }
                H = comph(PAR, QH, Stresses);
                omega = omega + H;
                var qg1 = new double[6];
                for (int j1 = 0; j1 < 6; j1++)
                {
                    qg1[j1] = Qgrad[j1];
                }
                ee = compds(PAR, qg1, QH, Stresses);
                for (int i = 0; i < 6; i++)
                    for (int j = 0; j < 6; j++)
                    {
                        paramConstMatr[i, j] = paramConstMatr[i, j] - dd[i] * ee[j] / omega;
                    }
                return paramConstMatr;
            }
        }
        public double[] compds(double [] PAR,double[] de,double [] QH,double[] Stresses)
        {
            var ds = new double[6];
            var Xk = 0.0;
            var shear = 0.0;
            Xk = QH[2] * Stresses[0] / PAR[3];
            shear = Xk * PAR[1];
            ds[0] = Xk * de[0];
            for (int i = 1; i < 6; i++)
            {
                ds[i] = shear * de[i];
            }
            return ds;
        }
        public double Quad (double [] PAR,double[] de,double [] QH,double[] Stresses, double [] Qgrad)
        {
            var value = 0.0;
            var dd = new double[6];
            dd = compds(PAR, de, QH, Stresses);
            for (int i=0;i<6;i++)
            {
                value = value + Qgrad[i] * dd[i];
            }
            return value;
        }
        #endregion
        #region ReadandWriteMethods
        public static void readData(string DataFileName, double[] array, int identifier)
        {
            //The identifier is for 3 possibilities for using it in the Constructor
            // 1 Deterministic Analysis and not using the readData for reading input file  0
            // 2 Deterministic Analysis and using the readData for reading input file      1
            // 3 Stochastic Analysis                                                       2
            string dataLine;
            string[] dataFields;
            string[] numSeparators1 = { ":" };
            string[] numSeparators2 = { " " };
            StreamReader rStream;
            rStream = File.OpenText(DataFileName);
            int dim = 1;
            dataLine = rStream.ReadLine();
            dataFields = dataLine.Split(numSeparators1, StringSplitOptions.RemoveEmptyEntries);
            dim = int.Parse(dataFields[1]);
            for (int i = 0; i < dim; i++)
            {
                dataLine = rStream.ReadLine();
                dataFields = dataLine.Split(numSeparators1, StringSplitOptions.RemoveEmptyEntries);
                array[i] = double.Parse(dataFields[1]);
            }
            rStream.Close();
        }
        public static void writeData(double[] array, int identifier)
        {
            string filename, dataLine;
            // The identifier is for telling if you want to write the whole array (1) or the last element (0) (for example the whole displacement curve or the last increment)
            // To insert spaces, use the simple space character " ", not tabs (i.e. "\t"). 
            // the editors do not 'interpret' the tabs in the same way, 
            // so if you open the file with different editors can be a mess.
            //string spaces1 = "        ";
            //string spaces2 = "              ";

            // format specifier to write the real numbers
            string fmtSpecifier = "{0: 0.0000E+00;-0.0000E+00}";

            StreamWriter wStream;
            filename = "displacements.txt";
            wStream = File.CreateText(filename);
            if (identifier == 1)
            {
                for (int i = 0; i < array.GetLength(0); i++)
                {
                    dataLine = String.Format(fmtSpecifier, array[i]);
                    wStream.WriteLine(dataLine);
                }
                wStream.Close();
            }
            else
            {
                dataLine = String.Format(fmtSpecifier, array[array.GetLength(0) - 1]);
                wStream.WriteLine(dataLine);
                wStream.Close();
            }
        }
        #endregion 
    }
}
