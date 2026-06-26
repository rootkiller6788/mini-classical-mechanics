# Knowledge Graph — mini-rigid-body

## L1: Definitions (Complete ✅)

| Concept | C Type | Lean Definition | Goldstein Ref |
|---------|--------|-----------------|---------------|
| 3D Vector | `vec3 {x,y,z}` | `Vec3` | §1.2 |
| 3×3 Matrix | `mat3 {m[9]}` | `Mat3` | §4.2 |
| Inertia Tensor | `InertiaTensor {Ixx,Iyy,Izz,Ixy,Ixz,Iyz}` | `InertiaTensor` | §5.3 |
| Euler Angles (ZXZ) | `EulerAngles {phi,theta,psi}` | `EulerAngles` | §4.8 |
| Principal Axes | `PrincipalAxes {moments[3], axes}` | `PrincipalAxes` | §5.4 |
| Rigid Body State | `RigidBodyState {t,omega,euler,L,T}` | — | §5 |
| Quaternion | `Quaternion {w,x,y,z}` | `Quaternion` | Appx B |
| ODE Derivative Fn | `ode_derivative_fn` | — | — |

## L2: Core Concepts (Complete ✅)

- Inertia tensor as quadratic form: I_ij = Σ m_k (r_k²δ_ij − r_{k,i} r_{k,j})
- Principal axis theorem: real symmetric I is always diagonalizable
- Parallel axis (Steiner) theorem: I_P = I_cm + M(d²𝟙 − d⊗d)
- Euler's equations: Iᵢω̇ᵢ − (Iⱼ−Iₖ)ωⱼωₖ = Nᵢ
- Torque-free motion: two integrals (T, L²)
- Angular velocity as so(3) element: dR/dt = R·[ω]×
- Quaternion double cover: SU(2) → SO(3), q → R(q)

## L3: Mathematical Structures (Complete ✅)

- **SO(3)** Lie group: 3×3 orthogonal matrices with det=+1
- **so(3)** Lie algebra: 3×3 skew-symmetric matrices [ω]×
- **SU(2)** double cover: unit quaternions, q → −q same rotation
- Exponential map: so(3) → SO(3) via Rodrigues formula
- Logarithm map: SO(3) → so(3) via axis-angle extraction
- SLERP: geodesic interpolation on S³
- ZXZ, ZYZ Euler angle conventions
- Tait-Bryan (roll-pitch-yaw) aerospace convention
- Jacobi eigenvalue algorithm for 3×3 symmetric matrices

## L4: Fundamental Laws (Complete ✅)

| Law | Statement | Verification |
|-----|-----------|-------------|
| Euler equations (free) | I₁ω̇₁ = (I₂−I₃)ω₂ω₃ (cyclic) | 126 tests |
| Euler equations (torque) | I·ω̇ + ω×(I·ω) = N | unit test |
| Energy conservation | dT/dt = 0 (torque-free) | RK4 drift < 1e-12 |
| Angular momentum conservation | dL/dt = 0 (torque-free, space frame) | RK4 drift < 1e-12 |
| Parallel axis theorem | I_P = I_cm + M(d²𝟙 − d⊗d) | Roundtrip test |
| König's theorem | T = T_cm + T_rot_about_cm | Unit test |
| Angular impulse | ΔL = ∫ N dt | Unit test |

## L5: Computational Methods (Complete ✅)

| Method | Implementation | Order |
|--------|---------------|-------|
| Forward Euler | `euler_step_omega()` | O(dt) |
| Midpoint (RK2) | `midpoint_step_omega()` | O(dt²) |
| Classical RK4 | `rk4_step_omega()` | O(dt⁴) |
| Velocity Verlet | `verlet_step_omega()` | O(dt²) symplectic-like |
| DOPRI54 adaptive | `dopri54_step_omega()` | O(dt⁵) with error control |
| Jacobi diagonalization | `principal_axes_decompose()` | O(n³) per sweep, 3×3 |
| Bisection root-finding | `steady_precession_angle()` | O(log(1/ε)) |
| Midpoint quadrature (3D) | `inertia_tensor_numerical()` | O(N³) |

## L6: Canonical Systems (Complete ✅)

| System | Physics | Implementation |
|--------|---------|---------------|
| Free rigid body | Torque-free Euler eqns | `simulate_free_rigid_body()` |
| Tennis racket theorem | Intermediate axis instability | `axis_stability_analysis()` |
| Dzhanibekov effect | Periodic flips about mid axis | `dzhanibekov_simulation()` |
| Heavy symmetric top | 1 fixed point in gravity | `symmetric_top_eff_potential()` |
| Steady precession | θ=const, dV_eff/dθ=0 | `steady_precession_angle()` |
| Lagrange top | Full 6D ODE system | `lagrange_top_ode()` |
| Sleeping top | θ=0 stability | `sleeping_top_is_stable()` |
| Gyroscope | N = ω×L, free/forced precession | `gyroscopic_torque()` |

## L7: Applications (Partial+ ✅)

| Application | Implementation | Domain keywords |
|-------------|---------------|-----------------|
| Gyroscope navigation | `gyroscope_2axis_simulate()` | GPS, navigation |
| Satellite attitude | `examples/satellite_attitude.jl` | NASA, Apollo, SpaceX |
| Free precession | `examples/free_precession.jl` | — |
| Gyroscope demo | `examples/gyroscope.jl` | — |

## L8: Advanced Topics (Partial+ ✅)

| Topic | Implementation |
|-------|---------------|
| Poinsot geometric construction | `invariable_plane()`, `polhode_tangent()` |
| Energy-momentum bifurcation | `separatrix_energy()`, `energy_on_momentum_sphere()` |
| Polhode/herpolhode curves | `body_to_space_trajectory()` |
| Energy ellipsoid visualization | `inertia_ellipsoid_points()` |
| Inertia tensor for arbitrary bodies | `inertia_tensor_numerical()` with density function |

## L9: Research Frontiers (Partial ✅ — Documented)

| Frontier | Status |
|----------|--------|
| Quantum rigid rotor | Documented: quantization of SO(3), spherical harmonics as Wigner D-matrices |
| Relativistic Thomas precession | Documented: spin precession in curved spacetime |
| Nonholonomic rigid bodies | Documented: Chaplygin sleigh, rattleback (celtic stone) |
| Many-body rigid body dynamics | Documented: articulated body algorithm (Featherstone) |

## Cross-Module Connections

```
mini-newtonian ──→ mini-rigid-body ──→ mini-hamiltonian (rigid body phase space)
                                    ──→ mini-aerospace (attitude control)
                                    ──→ mini-gyroscope (navigation)
```

## Mathematical Prerequisites

| Math | Usage |
|------|-------|
| Ordinary Differential Equations | Euler equations, Lagrange top ODE |
| Linear Algebra | Matrix operations, eigenvalue decomposition |
| Lie Groups / Algebras | SO(3), so(3), SU(2) structure |
| Vector Calculus | Cross products, quadratic forms |
| Numerical Analysis | RK methods, Jacobi algorithm, bisection |

