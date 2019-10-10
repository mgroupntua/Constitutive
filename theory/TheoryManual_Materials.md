# MGroup.Materials
In this chapter, the material laws available in the MSolve material library and the model properties classes are described.
## Material models

### Elastic isotropic 
It is an implementation of a linear elastic isotropic constitutive law
that can be expresses as 

<a href="https://www.codecogs.com/eqnedit.php?latex=_0^t\sigma=\mathbf{C}&space;\&space;_0^t\mathbf{e}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?_0^t\sigma=\mathbf{C}&space;\&space;_0^t\mathbf{e}" title="_0^t\sigma=\mathbf{C} \ _0^t\mathbf{e}" /></a>

in which 

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?_0^t\sigma="  /></a>  engineering stresses and   <a target="_blank"><img src="https://latex.codecogs.com/gif.latex?_0^te="  /></a> engineering strains 
when the elastic material is used with a standard small deformations continuum element formulation: 
```c#
public class ContinuumElement3D : IStructuralFiniteElement, ICell<Node>
```

or it can be expressed as:

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?_0^t\/S=\mathbf{C}&space;\&space;_0^t\mathbf{\epsilon}"  /></a>

in which 

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?_0^t\/S="  /></a>  second Piola-Kirchhoff stresses and   <a target="_blank"><img src="https://latex.codecogs.com/gif.latex?_0^t\epsilon="  /></a> Green Lagrange strains 
when the elastic material is used in a total Lagrangian formulation implemented in element type: 
```c#
public class Hexa8NonLinear : IStructuralFiniteElement, IEmbeddedHostElement
```
The constitutive matrix is given as 

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?C=\frac{E(1-\nu)}{(1+\nu)(1-2\nu)}\begin{bmatrix} 
1 & \frac{\nu}{1-\nu} & \frac{\nu}{1-\nu} & & & \\
\frac{\nu}{1-\nu} & 1 & \frac{\nu}{1-\nu} & & & \\
\frac{\nu}{1-\nu} & \frac{\nu}{1-\nu} & 1 & & & \\
& & & \frac{1- 2\nu}{2(1-\nu)} & & \\
& & & & \frac{1- 2\nu}{2(1-\nu)} & \\
& & & & & \frac{1- 2\nu}{2(1-\nu)} \\
\end{bmatrix}"  /></a>

for 3 dimensional problems and as 

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?C=\frac{E}{1-\nu^2}\begin{bmatrix} 
1 & \nu &  \\
\nu & 1 &  \\
& & \frac{1-\nu}{\nu}\\
\end{bmatrix}"  /></a>

or as

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?C=\frac{E(1-\nu)}{(1+\nu)(1-2\nu)}\begin{bmatrix} 
1 & \frac{\nu}{(1-\nu)} &  \\
\frac{\nu}{(1-\nu)} & 1 &  \\
& & \frac{1-2\nu}{2(1-\nu)}\\
\end{bmatrix}"  /></a>

for two dimensional plane stress and plane strain problems respectively.
The material constants used in the constitutive matrix C are:

  E= Young's modudlus, ν = Poisson's ratio.

References
- "Finite Element Procedures" Klaus-Jurgen Bathe, 1996 Prentice Hall, Inc. 

### Elastic isotropic with orientation trasnformations for shells

#### 1. Kirchhoff Love Shell Formulation 

According to the Kirchhoff Love shell element formulation implmented in Msolve(see ref.1 pp.38) the Green Lagrange strain coefficients
 Eij are computed from the metric coefficients in the actual and reference configuration as:

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?E_{ij}=\frac{1}{2}(g_{ij}-G_{ij})"  /></a>

and they refer to the contravariant basis <a target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{G^{i}}\otimes\boldsymbol{G^{j}}"  /></a> of the undeformed configuration.

For shell analysis the continuum is reduced to the midsurface of the shell and transverse normal stresses are neglected. For a Kirchhoff–
Love shell there is additionally the assumption that the director remains straight and normal to the midsurface during
deformation. This implies that the strain is assumed to be linear through the thickness and the transverse shear strains vanish,
Ea3=0 (see ref. 2 pp. 3906). 

The resulting constitutive matrix is identical to the plane stress 2D elasticity tensor when expressed 
in a local coordinate system, with normalised uints (ref.1 pp.43), hence it needs to be transformed to a 
to the same coordinate system as the green Lagrange strain coefficients by use of the following transformation:

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?\bar{E}_{\gamma\delta} =E_{\alpha \beta }(\boldsymbol{E}_{\gamma } \boldsymbol{G}^{\alpha})(\boldsymbol{G}^{\beta}\boldsymbol{E}_{\delta })"  /></a>

(ref. 3.41)

that results in the following transformation matrix:
<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{bmatrix} 
(\boldsymbol{E}_{1 } \boldsymbol{G}^{1})(\boldsymbol{G}^{1}\boldsymbol{E}_{1}) & (\boldsymbol{E}_{1} \boldsymbol{G}^{2})(\boldsymbol{G}^{2}\boldsymbol{E}_{1}) & (\boldsymbol{E}_{1} \boldsymbol{G}^{1})(\boldsymbol{G}^{2}\boldsymbol{E}_{1}) \\
(\boldsymbol{E}_{2} \boldsymbol{G}^{1})(\boldsymbol{G}^{1}\boldsymbol{E}_{2}) & (\boldsymbol{E}_{2} \boldsymbol{G}^{2})(\boldsymbol{G}^{2}\boldsymbol{E}_{2}) & (\boldsymbol{E}_{2} \boldsymbol{G}^{1})(\boldsymbol{G}^{2}\boldsymbol{E}_{2}) \\
2(\boldsymbol{E}_{2} \boldsymbol{G}^{1})(\boldsymbol{G}^{1}\boldsymbol{E}_{1}) & 2(\boldsymbol{E}_{2} \boldsymbol{G}^{2})(\boldsymbol{G}^{2}\boldsymbol{E}_{1}) & (\boldsymbol{E}_{2} \boldsymbol{G}^{1})(\boldsymbol{G}^{2}\boldsymbol{E}_{1})+(\boldsymbol{E}_{2} \boldsymbol{G}^{2})(\boldsymbol{G}^{1}\boldsymbol{E}_{1}) \\
\end{bmatrix}"  /></a>


References:

 - J. Kiendl  Isogeometric Analysis and Shape Optimal Design of Shell Structures, PhD thesis
 - J. Kiendl et al. / Comput. Methods Appl. Mech. Engrg. 198 (2009) 3902–3914


#### 2. Reissner Mindlin shell formulation

The element behaviour is based on the assumptions that, straight line defined by tge nodal director vectors (which, usually, give lines that in the
original configurationare close to normal to the midsurface of the shell) remain straight during the element deformations (but not nesessarily normal to the midsurfaec of
the shell) and that no transverse normal stress is developed in the dorections of the director vectors. 
The resulting constitutive matrix is given as (See ref. 1 Bathe):

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?C=\frac{E}{(1-\nu^2)}\begin{bmatrix} 
1 & \nu & 0 & 0 & 0 & 0\\
  & 1 & 0 & 0 & 0 & 0 \\
  &  & 0  & 0 & 0 & 0 \\
& & & \frac{1- \nu}{2} & 0 & 0 \\
& & & & k\frac{1- \nu}{2} & 0\\
& & & & & k\frac{1- \nu}{2} \\
\end{bmatrix}"  /></a>

and it refers to a local cartesian coordinate defined by the normal vectors and the tangent of
the midsurface of the shell.

In the formulation used it is transformed to a global cartesian system using the following transformation matrix:

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?Q_{sh}=\begin{bmatrix} 
l_1^2 & m_1^2 & n_1^2 & l_1m_1 & m_1n_1 & n_1l_1\\
l_2^2 & m_2^2 & n_2^2 & l_2m_2 & m_2n_2 & n_2l_2\\
l_3^2 & m_3^2 & n_2^2 & l_2m_2 & m_2n_2 & n_2l_2\\
2l_1l_2 & 2m_1m_2 & 2n_1n_2 & l_1m_2+l_2m_1 & m_1n_2+m_2n_1 & n_1l_2+n_2l_1\\
2l_2l_3 & 2m_2m_3 & 2n_2n_3 & l_2m_3+l_3m_2 & m_2n_3+m_3n_2 & n_2l_3+n_3l_2\\
2l_3l_1 & 2m_3m_1 & 2n_3n_1 & l_3m_1+l_1m_3 & m_3n_1+m_1n_3 & n_3l_1+n_1l_3\\
\end{bmatrix}"  /></a>

where:

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?l_1=cos(e_x,e_r) m_1=cos(e_y,e_r) n_1=cos(e_z,e_r)
l_2=cos(e_x,e_s)  m_2=cos(e_y,e_s) n_2=cos(e_z,e_s)
l_3=cos(e_x,e_t) m_3=cos(e_y,e_t) n_3=cos(e_z,e_t)
"  /></a>

the direction cosines of the local cartesian system basis vectors with respect to the
global cartesian coordinate system.

References
- "Finite Element Procedures" Klaus-Jurgen Bathe, 1996 Prentice Hall, Inc. 

### Bond slip traction separation law for cohesive material
This model describes an uncoupled behaviour in normal <a target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta_{\nu} "  /></a> and tangential direction
 <a target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta_{\tau}  "  /></a>.
In the normal direction a fully elastic behaviour is implemented. For the tangential separation an elastoplastic behaviour with ductile
unloading is exhibited. 
Shear stresses are calculated as:

 <a target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\tau=D:(\delta _{\tau}-\delta_{\tau}^p)"  /></a>.

The yield function is given as: 

<a href="https://www.codecogs.com/eqnedit.php?latex=f=\sqrt{(\tau_x-\alpha_x)^2&plus;(\tau_y-\alpha_y)^2}-\tau_{max}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f=\sqrt{(\tau_x-\alpha_x)^2&plus;(\tau_y-\alpha_y)^2}-\tau_{max}" title="f=\sqrt{(\tau_x-\alpha_x)^2+(\tau_y-\alpha_y)^2}-\tau_{max}" /></a>

Plastic straining occurs for:

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?f=0"  /></a> and <a target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dot{f}=0"  /></a>

It is assummed for plastic straining that 

<a target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dot{\delta}_p=\lambda m "  /></a> where <a target="_blank"><img src="https://latex.codecogs.com/gif.latex?m=\frac{\partial f }{\partial \tau}"  /></a>

### Thermal properties:
1. Specific Heat (Heat Capacity):
The heat capacity of a material is defined as the amount of heat required to raise its temperature by 1°. The heat capacity per unit mass, of material is defined as its specific heat. 

$ c= \frac{1}{m}\cdot \frac{dQ}{dT}

Where, m = Mass,
T = Temperature,
Q = Energy content, and
dQ = Energy (heat) added or subtracted to produce the temperature change dT.

For unit mass per degree change in temperature specific heat c = dQ, the quantity of heat that must be added per unit mass of a solid to raise its temperature by one degree. 

2. Thermal Conductivity:
It is defined as the amount of heat conducted in a unit time through a unit area normal to the direction of heat flow. Heat conduction through isotropic solids is expressed by Fourier’s law:

$ q=- k\frac{dT}{dx}

q = Rate of heat flow/unit area normal to the direction of flow,
T = Temperature,
x = Distance measured in the direction of flow, and
k = Thermal conductivity.
Heat flow through solids is due to elastic vibration of atoms or molecules or due to transfer of energy by the free electrons.

3. Density
The density, or more precisely, the volumetric mass density, of a substance is its mass per unit volume
Mathematically, density is defined as mass divided by volume:
 
$ \rho =\frac{m}{V}

Where, m = Mass,
V=volume

References 
1) Fundamentals of the finite element method for Heat and Fluid flow, R. Lewis, P. Nithiarasu and K. Seetharamu, Wiley, 2004
2) Thermal properties measurement of materials, Y. Jannot and A. Degiovanni, Wiley, 2018













