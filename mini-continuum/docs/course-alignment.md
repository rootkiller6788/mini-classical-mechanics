# Course Alignment — mini-continuum

## MIT 2.071/2.080 Mechanics of Materials / Structural Mechanics
- Stress tensor, strain tensor, Hooke's law: types.jl, elasticity.jl
- Elastic constants, plane stress/strain: elasticity.jl
- Euler-Bernoulli beams, Timoshenko beams: beams.jl
- Plate theory (Kirchhoff): plates.jl
- Buckling, stability: beams.jl, stability.jl
- Yield criteria, fatigue: failure.jl

## Stanford ME 340 Theory of Elasticity
- Finite deformation kinematics: types.jl (DeformationGradient, polar decomposition)
- Green-Lagrange, Almansi strain measures: types.jl, elasticity.jl
- Cauchy/Piola-Kirchhoff stress tensors: conservation.jl
- Anisotropic elasticity: elasticity.jl
- Conservation laws: conservation.jl

## Berkeley CE 230 Structural Mechanics
- Beam bending FDM solver: beams.jl
- Plate vibration, natural frequencies: plates.jl
- Shell theory (cylindrical): plates.jl
- Contact mechanics (Hertz): contact.jl
- Structural stability: stability.jl

## Caltech Ae/AM/CE/ME 102 Mechanics of Solids
- Continuum mechanics foundations: types.jl
- Elastic wave propagation: waves.jl
- Plasticity, J2 flow theory: plasticity.jl
- Fracture mechanics: fracture.jl

## Princeton MAE 542 Continuum Mechanics
- Balance laws (mass, momentum, energy): conservation.jl
- Constitutive theory (elastic, plastic, viscoelastic): elasticity.jl, plasticity.jl, viscoelasticity.jl
- Clausius-Duhem inequality: conservation.jl

## Cambridge Part II/III Theoretical Physics
- Landau-Lifshitz Theory of Elasticity: all modules
- Wave propagation in elastic solids: waves.jl

## Oxford C6 Solid Mechanics
- Poroelasticity (Biot theory): poroelasticity.jl
- Thermoelasticity: thermoelastic.jl
- Composites (micromechanics): composites.jl, micromechanics.jl

## ETH Zurich 151-0515 Continuum Mechanics
- Tensor algebra and calculus: types.jl
- Finite strain theory: types.jl
- Viscoelasticity: viscoelasticity.jl
- Creep mechanics: creep.jl

## Tokyo University Applied Solid Mechanics
- Failure criteria (Mohr-Coulomb, Drucker-Prager): failure.jl, plasticity.jl
- Fatigue analysis: failure.jl
- Micromechanics, dislocations: micromechanics.jl

## References
1. Landau & Lifshitz, Theory of Elasticity (Vol. 7)
2. Timoshenko & Goodier, Theory of Elasticity
3. Timoshenko & Gere, Theory of Elastic Stability
4. Achenbach, Wave Propagation in Elastic Solids
5. Johnson, Contact Mechanics
6. Anderson, Fracture Mechanics
7. Mura, Micromechanics of Defects in Solids
8. Christensen, Theory of Viscoelasticity
9. Wang, Theory of Linear Poroelasticity
10. Hill, The Mathematical Theory of Plasticity
