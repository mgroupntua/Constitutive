{% include mathjax.html %}
# MGroup.Multiscale
## Main capabilities 
1. Concurrent (FE2) Multiscale analysis of microscopically heterogeneous structures
2. Stress strain history analysis of heterogenous material models
### Multiscale analysis of structures
####  Theory

The multiscale analysis scheme implemented in Msolve is a two scale homogenization scheme, based on a scale
separation hypothesis, that results in a nested solution implementation. The material configuration to be considered is 
assummed to be macroscopically heterogeneous (continuum mechanics theory is suitable to describe the macroscopic behaviour), 
but microscopically heterogeneous (the morphology isconsists of distinguishablecomponents as e.g grains, interfaces, cavvities that
are described in finite element model of a representative volume element). Local periodicity is assummed, i.e. the microstructure can 
have differrent morphologies (and hence different RVE realizations) corresponding to different macroscopic point, while it repeats 
itslef in a small vicinity of each individual macroscopic point. 

The macroscopic material points are connected to the models of the microstructure based on two averaging relations:

The first one is imposed as a constraint for the bvp of the microstructure:

$ \bar{F}_{M}:=\frac{1}{V_{0}}\int_{y\in  V_{0}}^{ }F_mdV_{0}

and the other one is exploited for the calculation of the macroscopic stress tensor:

$ \bar{P}_{M}:=\frac{1}{V_{0}}\int_{y\in  V_{0}}^{ }P_mdV_{0}

In view of these relations then, the hill mandel macrohomogeneity condition is indeed satisfied and the local variation of work in 
the macroscale equals the volume average of variation of work performed on the RVE:

$ \bar{P}_{M}:\delta F_{M}=\frac{1}{V_{0}}\int_{V_{0}}^{ }P_m :  \delta F_{m}dV_{0}

In a nested solution scheme, calculation of stresses and the constitutive matrix of each macroscopi integration point for an iterative update
of the macroscopic deformation state requires first he solution of non linear in general boundary value problem in each rve. The iterative solution
scheme conecting both scales, implement in Msolve, can be found in [3],[5] and [6]
  

  - [1] Miehe C, Koch A. Computational micro-to-macro transitions of discretized
microstructures undergoing small strains. Arch Appl Mech 2002;72:300–17.
http://dx.doi.org/10.1007/s00419-002-0212-2.
  - [2] Kouznetsova V, Baaijens F. Approach to micro-macro modeling of heterogeneous materials,
 Comp. Mech. 2001 Jan. DOI: 10.1007/s004660000212 
  - [3] Kouznetsova V, Computational homogenization for the multi-scale analysis of multi-phase materials
 PhD thesis, Eidhoven University of Technology, DOI: 10.6100/IR560009
  - [4] Feyel F. A multilevel finite element method (FE2) to describe the response of
highly non-linear structures using generalized continua. Comput Methods
Appl Mech Eng 2003;192:3233–44. http://dx.doi.org/10.1016/S0045-7825(03)
00348-7.
  - [5]Schröder J. (2014) A numerical two-scale homogenization scheme: the FE2-method.
  In: Schröder J., Hackl K. (eds) Plasticity and Beyond. CISM International Centre for Mechanical Sciences, vol 550. Springer, Vienna
  - [6] Papadopoulos et al., The impact of interfacial properties on the macroscopic performance of
   carbon nanotube composites. A FE2-based multiscale study. Comp. Structures. 136 (2016) 582–592













