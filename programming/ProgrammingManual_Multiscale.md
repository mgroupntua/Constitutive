# MGroup.Multiscale
## Programming Manual 

#### Usage and programming guide
In the following instructions are given regarding the following analysis schemes available in the 
multiscale project of MSolve:
1. Concurrent (FE2) Multiscale analysis of microscopically heterogeneous structures
2. Stress strain history analysis of heterogenous material models
##### 1. Computational Homogenization for 2D problems-(stess-strain history analysis)
In this problem a finite element model representing the microstructure of a 2D 
structure is built. In a multiscale analysis context, described in over mentioned 
publications, appropriate boundary conditions (i.e. displacements), emanating from
a macroscopic strain state, are applied in the model.
After solving the microstructural boundary value problem, integration is performed
over the boundary of the rve model to calculate the macroscopic stress and constitutive
matrix.

###### Steps:
###### 1.1 Definition of the rve model 
The rve model is defined by use of an appropriate for 2D macrostructure rve-Builder,
which is an object of a class implementing the interface
```c#
public interface IdegenerateRVEbuilder : IRVEbuilder
```
In the present example we will make use of a homogeneous rve builder
```c#
IdegenerateRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderLinearAndDegenerate();
```
###### 1.2 Definition of the macroscopic material
Next, it is necessary to define an microstructure2D object that can be used as a material for
a 2D macroscopic model:
```c#
IContinuumMaterial2D microstructure3 = new Microstructure2DplaneStress(homogeneousRveBuilder1, 
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1);
```
###### 1.3 Imposition of the macroscopic strain state
The deformation state of the microstructure can be updated for an update of
the strain state of a macroscopic material point in the same way as for any other material:
```c#
microstructure3.UpdateMaterial(new double[3] { 0.010, 0, 0 });
```
###### 1.4 Calculation of macroscopic stress and constitutive tensor
Calculation of stress tensor and constitutive tensor is performed automatically after convergence of 
the bvp of the microstructure and then they can be accessed in a standard manner.
the strain state of a macroscopic material point in the same way as for any other material. 
For example stresses can be accessed as follows :
```c#
double[] stressesCheck3 = new double[3] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2] };
```
###### 1.5 Calculation of a stress - strain history 
Instead of applying more strain increments manually an appropriate method can be used
```c#
private static (double[][] stressHistory, double[][,] constitutiveMatrixHistory) StressStrainHistory(double[][] strainHistory, IContinuumMaterial2D testedMaterial)
```

###### 1.6 Entire example
*This example can be found in MSolve.Multiscale.Tests as:*
```c#
public static void Check2DscaleTransitionsAndMicrostructure()
```

##### 2. Computational Homogenization for 3D problems-(stess-strain history analysis)
In this problem a finite element model representing the microstructure of a 3D 
structure is built. Again, in a multiscale analysis context, described in the over mentioned 
publications, appropriate boundary conditions (i.e. displacements), emanating from
a macroscopic strain state, are applied in the model.
After solving the microstructural boundary value problem, integration is performed
over the boundary of the rve model to calculate the macroscopic stress and constitutive
matrix.

###### Steps:
###### 2.1 Definition of the rve model 
In this example a 3D rve builder should be used for the creation of the model of the 
microstructure:
```c#
IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderLinear();
```
which implements the following interface: 
```c#
public interface IRVEbuilder
```
###### 2.2 Definition of the macroscopic material
Next, it is necessary to define an microstructure3D object because in this case the heterogeneous
material of a 3D macrostructure is modeled
```c#
IContinuumMaterial3D microstructure3 = new Microstructure3D(homogeneousRveBuilder1, 
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1);
```
###### 2.3 Imposition of the macroscopic strain state
The deformation state of the microstructure can be updated for an update of
the strain state of a macroscopic material point in the same way as for any other material:
```c#
microstructure3.UpdateMaterial(new double[6] { 0.010, 0, 0, 0, 0, 0 });
```
###### 2.4 Calculation of macroscopic stress and constitutive tensor
Calculation of stress tensor and constitutive tensor is performed automatically after convergence of 
the bvp of the microstructure and then they can be accessed in a standard manner.
the strain state of a macroscopic material point in the same way as for any other material. 
For example stresses can be accessed as follows :
```c#
double[] stressesCheck3 = new double[6] { microstructure3.Stresses[0], microstructure3.Stresses[1], microstructure3.Stresses[2], microstructure3.Stresses[3], microstructure3.Stresses[4], microstructure3.Stresses[5] };
```
###### 2.5 Calculation of a stress - strain history 
Instead of applying more strain increments manually an appropriate method can be used
```c#
private static (double[][] stressHistory, double[][,] constitutiveMatrixHistory) StressStrainHistory(double[][] strainHistory, IContinuumMaterial3D testedMaterial)
```

###### 2.6 Entire example
*This example can be found in MSolve.Multiscale.Tests as:*
```c#
public static void Check3DscaleTransitionsAndMicrostructure()
```
##### 3. Computational Homogenization for thin shell problems-(stess-strain history analysis)
In this problem a finite element model representing the microstructure of a thin shell
structure,  is built. Appropriate boundary conditions (i.e. displacements), emanating from
a macroscopic strain state, and plane stress loading conditions, are applied in the model.
After solving the microstructural boundary value problem, integration is performed
over the boundary of the rve model to calculate the macroscopic stress and constitutive
matrix.

###### Steps:
###### 3.1 Definition of the rve model 
Similarily to the 1st example (2D macrostructure) the  rve model is defined by use of an appropriate for 2D macrostructure rve-Builder,
which is an object of a class implementing the interface
```c#
public interface IdegenerateRVEbuilder : IRVEbuilder
```
In the present example we will make use of a homogeneous rve builder
```c#
IdegenerateRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderLinearAndDegenerate();
```
###### 3.2 Definition of the macroscopic material
Next, it is necessary to define an MicrostructureShell2D object because in this case the heterogeneous
material of a shell macrostructure is modeled
```c#
var material4 = new MicrostructureShell2D(homogeneousRveBuilder1, 
                model => (new SkylineSolver.Builder()).BuildSolver(model), false, 1)
            {
                TangentVectorV1 = new double[3] { Vec1[0], Vec1[1], Vec1[2] },
                TangentVectorV2 = new double[3] { Vec2[0], Vec2[1], Vec2[2] }
            };
```
It is observed that in the case of a shell material it is also nesessary to define the orentation of the section of the shell.
The tangetial plane to the midsurface of the structure is defined by two tangent vectors in this example they are defined as: 
```c#
var Vec1 = Vector.CreateFromArray(new double[3] { 1, 0, 0 });
            var Vec2 = Vector.CreateFromArray(new double[3] { 0.5, 2, 0 });
```


###### 3.3 Imposition of the macroscopic strain state
The deformation state of the microstructure can be updated for an update of
the strain state of a macroscopic material point in the same way as for any other shell material:
```c#
var strain = new double[3] { 0.01, 0, 0 }; 
material4.UpdateMaterial(strain);
```
###### 3.4 Calculation of macroscopic stress and constitutive tensor
Calculation of stress tensor and constitutive tensor is performed automatically after convergence of 
the bvp of the microstructure and then they can be accessed in a standard manner.
the strain state of a macroscopic material point in the same way as for any other material. 
For example stresses can be accessed as follows :
```c#
double[] stressesCheck4 = new double[3] { material4.Stresses[0], material4.Stresses[1], material4.Stresses[2] };
```
###### 3.5 Calculation of a stress - strain history 
Instead of applying more strain increments manually an appropriate method can be used
```c#
private static (double[][] stressHistory, double[][,] constitutiveMatrixHistory) StressStrainHistory(double[][] strainHistory, IShellMaterial testedMaterial)
```

###### 3.6 Entire example
*This example can be found in MSolve.Multiscale.Tests as:*
```c#
public static void CheckShellScaleTransitionsAndMicrostructure()
```
##### 4. Use of the rve templates library
A rve builder can be chosen from a variety of model builders in the namespace:
namespace MGroup.Multiscale.RveTemplates
```c#
namespace MGroup.Multiscale.RveTemplates
```

###### Steps:
###### 4.1 Chosing an appropriate Microstructure
 Depending on the model of the macrostructure to be build (2d, 3d or shell etc.) an appropriate 
microstructure should be chosen
```c#
public class Microstructure2DplaneStress : StructuralProblemsMicrostructureBase, IContinuumMaterial2D
public class Microstructure3D : StructuralProblemsMicrostructureBase, IContinuumMaterial3D
public class MicrostructureShell2D : StructuralProblemsMicrostructureBase, IShellMaterial 
```

###### 4.2 Chosing the Rve builder 
In the constructor of a microstructure class it is revealed what kind of a rveBuilder the microstructure object can
be used with (for example: IdegenerateRVEbuilder, IRVEbuilder etc.)
```c#
public Microstructure2DplaneStress(IdegenerateRVEbuilder rveBuilder, Func<Model, ISolver> createSolver, bool EstimateOnlyLinearResponse, int database_size)
public Microstructure3D(IRVEbuilder rveBuilder, Func<Model, ISolver> createSolver, bool EstimateOnlyLinearResponse, int database_size)
public MicrostructureShell2D(IdegenerateRVEbuilder rveBuilder, Func<Model, ISolver> createSolver, 
            bool EstimateOnlyLinearResponse, int database_size)
```
After identifying the interface the rveBuilder should implement an appropriate rve builder can be chosen from the RveTemplates folder

##### 5. Concurrent multi-scale modeling of structures (FEn analysis)
In this problem a concurrent two scale analysis (FE2) is performed for a small cantilever
structure, consisting of hexa8 elements. Appropriate boundary conditions (i.e. displacements), emanating from
each macroscopic strain state are applied in the RVE models assigned on each integration point of the macrostructure.
The macroscopic deformation is updated consecutively until equilibrium for a certain external load is reached. 
For each trial iterative update of the macroscopic deformation state an microstructural non linear bvp is solved foreach macroscopic
integration point.

###### Steps:
###### 5.1 Definition of the rve model (material object) 
Following the steps described in examples (1-4) we difine the models of the microstructure that in this certain example
will be used in each integration point of the macrostructure during FE2 multiscale analysis.
We start by definining an rve builder
```c#
IRVEbuilder homogeneousRveBuilder1 = new HomogeneousRVEBuilderNonLinear();
```
And then the macrostructure object that will be used as a material 
```c#
IContinuumMaterial3DDefGrad material1 = new MicrostructureDefGrad3D(homogeneousRveBuilder1,
                m => (new SkylineSolver.Builder()).BuildSolver(m), false, 1);
```
###### 5.2 Definition of the macroscopic model
Next the model of the macrostructure is built in a standard manner but instead of using a material 
we assign in each element the above defined micro-model:
```c#
Element e1 = new Element()
                {
                    ID = nElement + 1,
                    ElementType = new Hexa8NonLinearDefGrad(material1, GaussLegendre3D.GetQuadratureWithOrder(2, 2, 2)) 
                };
```
###### 5.3 definition of a nonlinear analysis
it is nesessary to define a nonlinear analyzer that will handle the nonlinear analysis a the macroscopic level:
```c#
var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);
```
To this point the definition of the non linear problem in the macroscopic level is not different from the definition of 
a standard non-linear problem. oOn the microscopic level it is not nesessary to define an analyzer, as the non-llinear analysis will be handled
by the microstructure class and its analyzers.

###### 5.4 Solution of the two-level multiscale problem
The solution process of the macroscopic problem, and hence the two level multiscale problem follows:
```c#
parentAnalyzer.Initialize();
parentAnalyzer.Solve();
```
###### 5.6 Entire example
*This example can be found in MSolve.Multiscale.Tests as:*
```c#
public static void CheckElasticMultiscaleCantilever()
```













