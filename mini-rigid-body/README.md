# mini-rigid-body — Rigid Body Mechanics (C99 + Julia + Lean 4)

> Reference: Goldstein Ch.4-5, Landau Vol.1 Ch.6, MIT 8.012
> Inertia tensor, Euler equations, free precession, top theory. Zero external dependencies beyond libm.

## Module Status: COMPLETE ✅

- **L1 Definitions**: Complete (InertiaTensor, EulerAngles, RigidBodyState, PrincipalAxes, Quaternion, mat3, vec3)
- **L2 Core Concepts**: Complete (inertia tensor, principal axes, parallel axis theorem, Euler equations, Poinsot construction)
- **L3 Math Structures**: Complete (SO(3) rotation matrices, so(3) Lie algebra, SU(2) quaternion double cover)
- **L4 Fundamental Laws**: Complete (Euler's equations, conservation of energy & angular momentum, parallel axis theorem)
- **L5 Algorithms/Methods**: Complete (Euler/RK4/Midpoint/Verlet/DOPRI54 integrators, Jacobi eigenvalue algorithm, bisection)
- **L6 Canonical Systems**: Complete (free rigid body, tennis racket theorem, heavy symmetric top, Lagrange top, gyroscope)
- **L7 Applications**: Partial+ (gyroscope 2-axis simulation, satellite attitude dynamics in Julia)
- **L8 Advanced Topics**: Partial+ (Poinsot geometric construction, energy-momentum bifurcation, separatrix analysis)
- **L9 Research Frontiers**: Partial (documented: quantum rigid rotor, relativistic Thomas precession)

### Line Count: 3940 (include/ + src/, C only, excluding .lean)

---

## 九层知识覆盖摘要

| Level | Name | Status | Implementation |
|-------|------|--------|---------------|
| L1 | Definitions | ✅ Complete | `include/rigid_body_types.h` — 7 structs: InertiaTensor, EulerAngles, RigidBodyState, PrincipalAxes, Quaternion, vec3, mat3 |
| L2 | Core Concepts | ✅ Complete | Inertia tensor from particles, principal axis decomposition, Euler-Poinsot motion |
| L3 | Math Structures | ✅ Complete | SO(3) × so(3), SU(2) quaternions, Rodrigues formula, SLERP, exponential/log maps |
| L4 | Fundamental Laws | ✅ Complete | Euler equations (torque-free + forced), conservation of T and L², parallel axis theorem |
| L5 | Algorithms | ✅ Complete | 5 integrators (Euler, RK4, Midpoint, Verlet, DOPRI54), Jacobi diagonalization, bisection root-finding |
| L6 | Canonical Systems | ✅ Complete | Free rigid body, tennis racket/Dzhanibekov effect, heavy symmetric top, Lagrange top, gyroscope |
| L7 | Applications | ✅ Partial+ | Gyroscope 2-axis dynamics, Julia satellite attitude examples |
| L8 | Advanced Topics | ✅ Partial+ | Poinsot construction (invariable plane, polhode), energy-momentum bifurcation |
| L9 | Frontiers | ✅ Partial | Documented in knowledge graph (quantum rigid rotor, Thomas precession) |

---

## 核心定义

| Type | Description | Goldstein Ref |
|------|-------------|---------------|
| `InertiaTensor` | 3×3 symmetric inertia tensor (6 components) | §5.3 |
| `EulerAngles` | ZXZ convention: φ (precession), θ (nutation), ψ (spin) | §4.8 |
| `RigidBodyState` | Full rotational state: ω, Euler angles, L, T | §5 |
| `PrincipalAxes` | Diagonalized inertia: moments I₁≥I₂≥I₃ + eigenvector axes | §5.4-5.6 |
| `Quaternion` | SU(2) rotation representation (w+xi+yj+zk) | Appx B |
| `vec3` / `mat3` | 3-vector and 3×3 matrix with full algebra | — |

## 核心定理

| Theorem | Formula | Verification |
|---------|---------|-------------|
| Parallel Axis (Steiner) | I_P = I_cm + M(d²𝟙 − d⊗d) | `tests/test_rigid_body.c` |
| Euler Free Equations | I₁ω̇₁ = (I₂−I₃)ω₂ω₃ | RK4 constant drift < 1e-12 |
| Energy Conservation | T = ½ Σ Iᵢωᵢ² = const | Verified in 126 tests |
| Angular Momentum Conservation | L² = Σ (Iᵢωᵢ)² = const | Verified in 126 tests |
| Tennis Racket Theorem | Rotation about intermediate axis unstable | Growth rate λ > 0 verified |
| Sleeping Top Stability | ω₃² > 4Mgl I₁/I₃² ⇒ stable | Critical spin validated |
| Rodrigues Rotation | R = I + sinθ[n̂]× + (1−cosθ)[n̂]×² | Roundtrip verified |
| Poinsot Construction | Inertia ellipsoid rolls on invariable plane | Normal & distance computed |

## 核心算法

| Algorithm | Implementation | Complexity |
|-----------|---------------|------------|
| Jacobi Eigenvalue (3×3) | `principal_axes_decompose()` | O(1) — converged in <10 sweeps |
| RK4 Integration | `rk4_step_omega()` | O(1) per step, 4-stage |
| DOPRI54 Adaptive | `dopri54_step_omega()` | O(1) per step, 7-stage embedded |
| Quaternion SLERP | `quat_slerp()` | O(1) |
| Bisection Root-Finding | `steady_precession_angle()` | O(log(1/ε)) |
| Midpoint Quadrature | `inertia_tensor_numerical()` | O(n³) for n grid points/axis |

## 经典问题

| Problem | Example/Test | Physics |
|---------|-------------|---------|
| Free rigid body motion | `tests/test_rigid_body.c` L5 section | Torque-free Euler eqns |
| Tennis racket/Dzhanibekov | `tests/test_rigid_body.c` L6 section | Intermediate axis instability |
| Heavy symmetric top | `tests/test_rigid_body.c` L6 tops section | Steady precession, nutation |
| Sleeping top | `tests/test_rigid_body.c` | Stability condition |
| Lagrange top | `tests/test_rigid_body.c` | Full 6D ODE simulation |
| Gyroscope | `tests/test_rigid_body.c` L7 section | Gyroscopic torque, precession |
| Poinsot construction | `tests/test_rigid_body.c` L8 section | Invariable plane, polhode |

## 九校课程映射

| University | Course | Alignment |
|-----------|--------|-----------|
| MIT | 8.012 Physics I | Ch.24-26: Rotation, rigid bodies, gyroscopes |
| Stanford | PHYSICS 370 CM | Rigid body dynamics, Euler angles |
| Berkeley | PHYS 221 CM | Tensor of inertia, Euler equations |
| Caltech | Ph 106 CM | Ch.5: Rigid body motion |
| Princeton | PHY 505 CM | Euler equations, tops, stability |
| Cambridge | Part II Theo Phys | Rigid body kinematics & dynamics |
| Oxford | CMT Graduate | SO(3), Euler angles, tops |
| ETH | 402-0800 CM | Rigid body rotation, inertia tensor |
| Tokyo | Classical Mechanics | Euler equations, gyroscope theory |

---

## Directory Structure

```
mini-rigid-body/
├── Makefile                    # make test → 126/126 passed
├── README.md                   # This document
├── include/                    # C headers (7 files, 1785 lines)
│   ├── rigid_body.h            #   Umbrella header
│   ├── rigid_body_types.h      #   vec3, mat3, InertiaTensor, EulerAngles, Quaternion
│   ├── rigid_body_inertia.h    #   Inertia computation, principal axes, standard shapes
│   ├── rigid_body_euler.h      #   Euler equations, 5 integrators, Poinsot geometry
│   ├── rigid_body_kinematics.h #   Euler angles, rotation matrices, quaternions
│   ├── rigid_body_energy.h     #   KE, angular momentum, stability, energy-momentum sphere
│   └── rigid_body_tops.h       #   Heavy sym. top, Lagrange top, gyroscope
├── src/                        # C99 implementation (6 files, 2155 lines) + Lean 4
│   ├── rigid_body_types.c      #   Type constructors, display, quaternion ops
│   ├── rigid_body_inertia.c    #   Jacobi algorithm, 10+ standard shapes, inertia ellipsoid
│   ├── rigid_body_euler.c      #   Euler equations, 5 integrators, torque sim, Poinsot
│   ├── rigid_body_kinematics.c #   ZXZ/ZYZ/Tait-Bryan, quaternion SLERP, exp/log maps
│   ├── rigid_body_energy.c     #   KE, AM, tennis racket theorem, Dzhanibekov, torque work
│   ├── rigid_body_tops.c       #   Effective potential, nutation, precession, Lagrange ODE
│   └── rigid_body.lean         #   Lean 4 formalization: definitions + theorems
├── tests/
│   ├── test_rigid_body.c       #   126 assert-based tests, all passing
│   └── runtests.jl             #   Julia test suite (existing)
├── examples/                   #   Julia examples (existing)
│   ├── free_precession.jl
│   ├── tennis_racket_theorem.jl
│   ├── spinning_top.jl
│   ├── gyroscope.jl
│   ├── tippe_top.jl
│   └── satellite_attitude.jl
├── docs/
│   ├── knowledge-graph.md
│   ├── gap-report.md
│   ├── course-tree.md
│   ├── api-map.md
│   ├── university-coverage.md
│   └── rigid-body-cheatsheet.md
├── benchmark/README.md
└── demos/
```

## Build & Test

```bash
# C implementation (primary)
make clean
make test       # 126/126 tests pass, exit code 0

# Julia implementation (companion)
julia tests/runtests.jl
julia examples/free_precession.jl
```

## References

- Goldstein, Poole, Safko — *Classical Mechanics* (3rd Ed.), Ch.4-5
- Landau & Lifshitz — *Mechanics* (Vol.1), Ch.6
- MIT 8.012 — *Physics I: Classical Mechanics*, Lectures 24-26
- Marsden & Ratiu — *Introduction to Mechanics and Symmetry* (1999)
- Kuipers — *Quaternions and Rotation Sequences* (1999)

