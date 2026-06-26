# Knowledge Graph — mini-continuum

## L1: Definitions
- StressTensor (Cauchy stress tensor components)
- StrainTensor (small-strain tensor components)
- ElasticMaterial (isotropic elastic constants E, nu, rho)
- DeformationGradient (F = dx/dX, finite deformation)
- RotationTensor (axis-angle representation)
- MaterialSymmetry (isotropic → triclinic)
- StressInvariants (I1, J2, J3)
- SymTensor2 (generic symmetric 2nd-order tensor)

## L2: Core Concepts
- Generalized Hooke's law (isotropic, orthotropic, anisotropic)
- Elastic constant conversions (E,nu ↔ lambda,mu ↔ K,G)
- Plane stress / plane strain
- Elastic waves (P, S, Rayleigh, Love, Lamb, Stoneley)
- Beam theory (Euler-Bernoulli, Timoshenko)
- Plate/shell theory (Kirchhoff, Mindlin)
- Yield criteria (von Mises, Tresca, Mohr-Coulomb, Drucker-Prager)
- Finite deformation (polar decomposition, strain measures)
- Conservation laws (mass, momentum, energy)
- Clausius-Duhem inequality

## L3: Mathematical Structures
- Symmetric 2nd-order tensor algebra
- 3x3 matrix operations (inverse, determinant, eigenvalues)
- Voigt notation (6-component ↔ 3x3 tensor)
- Polar decomposition (Newton iteration)
- Tensor invariants
- Bond transformation (stiffness rotation)
- 6x6 stiffness/compliance matrices

## L4: Fundamental Laws
- Hooke's law (generalized)
- Cauchy's equations of motion
- Clausius-Duhem inequality (2nd law of thermodynamics)
- Fourier's law of heat conduction
- Euler buckling formula
- Griffith fracture criterion
- Paris fatigue crack growth law
- Basquin/Coffin-Manson fatigue laws

## L5: Computational Methods
- 1D FDTD for wave equation
- FDM for Euler-Bernoulli beam
- FDM for Terzaghi consolidation
- FDM for heat equation
- Radial return algorithm (J2 plasticity)
- Boltzmann superposition integral
- Eigenvalue computation (Cardano's method)

## L6: Canonical Systems
- Simply supported / cantilever beams
- Timoshenko beam (shear deformation)
- Kirchhoff circular plate (SS/clamped)
- Rectangular plate vibration
- Euler column buckling
- Hertz contact (sphere-sphere)
- Elastic wave reflection/transmission
- Membrane vibration
- Cylindrical shell

## L7: Applications
- Seismic wave analysis (P/S/Rayleigh speeds, Snell)
- Structural beam design (FDM, uniform/point loads)
- Column stability (multiple BCs)
- Plate vibration design (frequency, boundary effects)
- Fatigue life prediction (Basquin, Miner)
- Creep life assessment (Larson-Miller, Norton)
- Fracture safety assessment (K_I vs K_IC)

## L8: Advanced Topics
- Finite deformation theory (Green-Lagrange, Almansi, Hencky)
- Anisotropic elasticity (orthotropic, transversely isotropic)
- Composite mechanics (Rule of Mixtures, Halpin-Tsai, laminate ABD)
- Plasticity (J2, Drucker-Prager, Cam-Clay, hardening)
- Viscoelasticity (Maxwell, Kelvin-Voigt, SLS, Prony series)
- Contact mechanics (Hertz, JKR, DMT, Mindlin)
- Fracture mechanics (K-factor, J-integral, CTOD, Weibull)
- Micromechanics (Eshelby, Mori-Tanaka, dislocations, Hall-Petch)
- Poroelasticity (Biot theory, Terzaghi consolidation, Gassmann)
- Thermoelasticity (thermal stress, TED, Biot number)

## L9: Research Frontiers
- Fractional derivative viscoelasticity
- Kachanov continuum damage mechanics
- Eshelby configurational forces
- Crystal plasticity (Schmid law, Taylor hardening)
- Nano-mechanics (Peierls stress)
- Multi-scale homogenization

## Cross-module Links
```
mini-newtonian → mini-continuum → mini-fluid-dynamics
                                → mini-solid-state
                                → mini-computational-physics
```
