using System.Collections.Generic;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Providers;

namespace MGroup.Constitutive.PorousMedia
{
    public class ElementPoreDampingProvider : IElementMatrixProvider
    {
        private readonly IElementMatrixProvider solidDampingProvider;
        private readonly double dampingCoefficient;

        public ElementPoreDampingProvider(IElementMatrixProvider solidDampingProvider, double dampingCoefficient)
        {
            this.solidDampingProvider = solidDampingProvider;
            this.dampingCoefficient = dampingCoefficient;
        }

        private IMatrix PorousMatrix(IElementType element)
        {
            IPorousElementType elementType = (IPorousElementType)element;
            int dofs = 0;
            foreach (IList<IDofType> dofTypes in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                foreach (IDofType dofType in dofTypes) dofs++;
            var poreDamping = SymmetricMatrix.CreateZero(dofs);

            IMatrix damping = solidDampingProvider.Matrix(element);
            IMatrix saturation = elementType.SaturationMatrix();
            IMatrix coupling = elementType.CouplingMatrix();

            int matrixRow = 0;
            int solidRow = 0;
            int fluidRow = 0;
            foreach (IList<IDofType> dofTypesRow in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                foreach (IDofType dofTypeRow in dofTypesRow)
                {
                    int matrixCol = 0;
                    int solidCol = 0;
                    int fluidCol = 0;
                    foreach (IList<IDofType> dofTypesCol in elementType.DofEnumerator.GetDofTypesForMatrixAssembly(element))
                        foreach (IDofType dofTypeCol in dofTypesCol)
                        {
                            if (dofTypeCol == PorousMediaDof.Pressure)
                            {
                                if (dofTypeRow == PorousMediaDof.Pressure)
                                    poreDamping[matrixRow, matrixCol] = -saturation[fluidRow, fluidCol];
                                else
                                    poreDamping[matrixRow, matrixCol] = coupling[fluidCol, solidRow];
                                fluidCol++;
                            }
                            else
                            {
                                if (dofTypeRow != PorousMediaDof.Pressure)
                                    poreDamping[matrixRow, matrixCol] = damping[solidRow, solidCol] * dampingCoefficient;
                                else
                                    poreDamping[matrixRow, matrixCol] = coupling[fluidRow, solidCol];
                                solidCol++;
                            }
                            matrixCol++;
                        }

                    if (dofTypeRow == PorousMediaDof.Pressure)
                        fluidRow++;
                    else
                        solidRow++;
                    matrixRow++;
                }

            return poreDamping;
        }

        #region IElementMatrixProvider Members

        public IMatrix Matrix(IElementType element)
        {
            if (element is IPorousElementType)
                return PorousMatrix(element);
            else
            {
                IMatrix dampingMatrix = solidDampingProvider.Matrix(element);
                dampingMatrix.Scale(dampingCoefficient);
                return dampingMatrix;

            }
        }

        #endregion
    }
}
