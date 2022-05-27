using System.Collections.Generic;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.PorousMedia
{
    public class ElementPoreMassProvider : IElementMatrixProvider
    {
        private readonly IElementMatrixProvider solidMassProvider;
        private readonly double massCoefficient;

        public ElementPoreMassProvider(IElementMatrixProvider solidMassProvider, double massCoefficient)
        {
            this.solidMassProvider = solidMassProvider;
            this.massCoefficient = massCoefficient;
        }

        private IMatrix PorousMatrix(IElementType element)
        {
            IPorousElementType elementType = (IPorousElementType)element;
            int dofs = 0;
            foreach (IList<IDofType> dofTypes in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                foreach (IDofType dofType in dofTypes) dofs++;
            var poreMass = SymmetricMatrix.CreateZero(dofs);

            IMatrix mass = solidMassProvider.Matrix(element);

            int matrixRow = 0;
            int solidRow = 0;
            foreach (IList<IDofType> dofTypesRow in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                foreach (IDofType dofTypeRow in dofTypesRow)
                {
                    int matrixCol = 0;
                    int solidCol = 0;
                    foreach (IList<IDofType> dofTypesCol in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                        foreach (IDofType dofTypeCol in dofTypesCol)
                        {
                            if (dofTypeCol == PorousMediaDof.Pressure)
                            {
                            }
                            else
                            {
                                if (dofTypeRow != PorousMediaDof.Pressure)
                                    poreMass[matrixRow, matrixCol] = mass[solidRow, solidCol] * massCoefficient;
                                solidCol++;
                            }
                            matrixCol++;
                        }

                    if (dofTypeRow != PorousMediaDof.Pressure) solidRow++;
                    matrixRow++;
                }

            return poreMass;
        }

        #region IElementMatrixProvider Members

        public IMatrix Matrix(IElementType element)
        {
            if (element is IPorousElementType)
                return PorousMatrix(element);
            else
            {
                IMatrix massMatrix = solidMassProvider.Matrix(element);
                massMatrix.Scale(massCoefficient);
                return massMatrix;
            }
        }

        #endregion
    }
}
