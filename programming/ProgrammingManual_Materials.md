# MGroup.Materials
## Programming Manual 

#### Usage and programming guide
In the following, instructions are given regarding the definition of material models and general dynamic properties 
for each element of the model
##### 1. General instructions 
Upon definition of an element of the model the material model of the element is chosen and it is given as
an argument in the constructor of the element. The material object is usually cloned inside the created element
when more than one gauss points are used for integration and it is not nesessary to assign the material in each 
gauss point separately. For example:

```c#
ElasticMaterial3D material1 = new ElasticMaterial3D()
                {
                    YoungModulus = 3.5,
                    PoissonRatio = 0.4,
                };

Element e1 = new Element()
                {
                    ID = hexaID,
                    ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };
Element e2 = new Element()
                {
                    ID = hexaID,
                    ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };
```

In this code segmet an elastic material is defined and the it is defined in hexa8 finite element only once.
The same material "material1" can then be assigned to another element as it will be cloned for the gauss points of
the second element as well, and its properties will not be overwritenn during the analysis by two elements. Separate
clonde objects of class material will be created and used for separate gauss points automatically.



######  1. Specific material models:
Next we give instructions on how to define and use a material object using the material objects of Msolve.Constitutive library.
###### 1.1 Elastic Material 3D
It will now be demonstrated how to define an elastic isotropic material for 3D problems which is an object of class:
```c#
public class ElasticMaterial3D : IIsotropicContinuumMaterial3D
```

Using the constructor of the class the Young's modulus and the Poisson ratio values are defined
```c#
ElasticMaterial3D material1 = new ElasticMaterial3D()
                {
                    YoungModulus = 3.5,
                    PoissonRatio = 0.4,
                };
```

###### 1.2 Elastic Material 2D
It will now be demonstrated how to define an elastic isotropic material for 2D problems which is an object of class:
```c#
public class ElasticMaterial2D : IIsotropicContinuumMaterial2D
```

Using the constructor of the class the Young's modulus and the Poisson ratio values are defined and it is also defined
that this material will be used in a 2D plane stress analysis.
```c#
var material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };
```
In the same manner, a material for plane strain analysis is defined
```c#
var material = new ElasticMaterial2D(StressState3D.PlaneStrain)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };
```
###### 1.3 Elastic isotropic Material for Kirchhoff Love shell formulations
It will now be demonstrated how to define an elastic isotropic material for 2D problems which is an object of class:
```c#
public class ShellElasticMaterial2Dtransformationb : IShellMaterial
```

Using the constructor of the class the Young's modulus and the Poisson ratio values are defined and it is also defined
that this material will be used in a 2D plane stress analysis.
```c#
var material3 = new ShellElasticMaterial2Dtransformationb() { YoungModulus = E_disp, PoissonRatio = ni_disp, TangentVectorV1 = new double[3] { 1,0,0 }, TangentVectorV2 = new double[3] { 0,1,0 } };
```
The vectors tangent to the midsurface of the shell element are also defined.
The correct values of the shell midsurface tangent vectors will be calculated in the shell element and then assigned to the cloned 
ShellElasticMaterial2Dtransformationb object and the values given above will not be used. Unless the material object is defined 
for purpose other than an analysis of a model, and it is defined separately so that its tangent moduli can be calculated etc.

###### 1.4 Elastic isotropic Material for Reisner Mindlin shell formulations
It will now be demonstrated how to define an elastic isotropic material for 2D problems which is an object of class:
```c#
public class ShellElasticMaterial3D : IShellMaterial 
```

Using the constructor of the class the Young's modulus and the Poisson ratio values are defined and it is also defined
that this material will be used in a 2D plane stress analysis.
```c#
var material3 new ShellElasticMaterial3D()
            {
               YoungModulus=3.4,
               PoissonRatio=0.5,
               ShearCorrectionCoefficientK= (double)(5/6)              
            };
```
The vectors tangent and normal to the midsurface of the shell element can also be defined (public properties).
```c#
public double[] NormalVectorV3 { get; set; }
public double[] TangentVectorV1 { get; set; }
public double[] TangentVectorV2 { get; set; }
```
As it is explained above, the correct values of the shell midsurface tangent vectors will be calculated in the shell element and then assigned to the cloned 
ShellElasticMaterial3D object and the values given above will not be used. Unless the material object is defined 
for purpose other than an analysis of a model, and it is defined separately so that its stifness matrix can be calculated etc.

###### 1.5 Elastic Material (Deformation Gradient based)
It will now be demonstrated how to define an elastic isotropic material for 3D problems, using the deformation gradient based formulation, which is an object of class:
```c#
public class ElasticMaterial3DDefGrad : IContinuumMaterial3DDefGrad
```

Using the constructor of the class the Young's modulus and the Poisson ratio values are defined
```c#
IContinuumMaterial3DDefGrad material1 = new ElasticMaterial3DDefGrad() { PoissonRatio = 0.3, YoungModulus = 1353000 };
```

###### 1.6 Bond Slip Traction Separation Law
It will now be demonstrated how to define an bond slip traction separation law for cohesive elements, which is an object of class:
```c#
public class BondSlipCohMat : ICohesiveZoneMaterial3D
```

Using the constructor of the class the stiffness in the shear separation , the stifness after yielding,
the normal separation stiffness , the yield shear stress and the initial stress and hardening values , and the tolerance of the internal
iterative procedure for the calculation of the next stress-strain point are defined
```c#
BondSlipCohMat material1 = new BondSlipCohMat(100, 10, 100, 100, new double[2], new double[2], 1e-10);
```

###### 2. Useful methods 
Next some useful methods of the material classes are presented. They are called by the element classes but they can be usefull in general.


###### 2.1 Imposition of the macroscopic strain state
The deformation state of the microstructure can be updated for an update of
the strain state of a macroscopic material point using the update state method:
```c#
material1.UpdateMaterial(new double[3] { 0.010, 0, 0 });
```
###### 1.4 Calculation of macroscopic stress and constitutive tensor
Calculation of stress tensor and constitutive tensor follows the strain state update of the material.
They can be accesed as follows:
```c#
double[] stressesCheck3 = new double[3] { material1.Stresses[0], material1.Stresses[1], material1.Stresses[2] };
constitutive[m, n] = material1.ConstitutiveMatrix[m, n];
```














