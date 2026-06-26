# 课程对标 — mini-lagrangian (Course Alignment)

## 九校课程章节映射

### MIT — 8.012 Physics I: Classical Mechanics
| Goldstein Chapter | Topic | Implementation |
|-------------------|-------|---------------|
| Ch.1 — Survey of Elementary Principles | Generalized coordinates, constraints | `generalized.jl`, `constraints.jl` |
| Ch.2 — Variational Principles | Hamilton principle, Euler-Lagrange eq | `euler_lagrange.jl`, `action.jl`, `variational_calculus.jl` |
| Ch.3 — Central Force Problem | Kepler orbits, scattering | `central_force.jl` |
| Ch.4 — Kinematics of Rigid Body | Euler angles, SO(3) | `rigid_body.jl` |
| Ch.5 — Rigid Body Equations | Euler equations, heavy top | `rigid_body.jl` |
| Ch.6 — Small Oscillations | Normal modes, forced response | `small_oscillations.jl` |
| Ch.8 — Legendre Transform | L→H, Poisson brackets | `legendre.jl` |
| Ch.11 — Canonical Perturbation | Secular terms, KAM | `perturbation.jl` |
| Ch.13 — Noether Theorem | Symmetries, conservation laws | `noether.jl` |

### Stanford — PHYSICS 370 Classical Mechanics
- Lagrangian & Hamiltonian formalism: `euler_lagrange.jl`, `legendre.jl`
- Action principles: `action.jl`
- Noether theorem: `noether.jl`
- Rigid body dynamics: `rigid_body.jl`

### Berkeley — PHYS 242 Classical Mechanics
- Variational principles: `variational_calculus.jl`, `action.jl`
- Constrained systems: `constraints.jl`
- Normal modes: `small_oscillations.jl`

### Caltech — Ph 106 Classical Mechanics
- Calculus of variations: `variational_calculus.jl`
- Central forces: `central_force.jl`
- Coupled oscillations: `small_oscillations.jl`

### Princeton — PHY 505 Classical Mechanics
- Lagrangian dynamics: full module coverage
- Symmetries & Noether: `noether.jl`

### Cambridge — Part II Theoretical Physics
- Lagrangian & Hamiltonian methods: `legendre.jl`
- Small oscillations: `small_oscillations.jl`
- Perturbation theory: `perturbation.jl`

### Oxford — Classical Mechanics (Graduate)
- Action principles: `action.jl`
- Rigid bodies: `rigid_body.jl`
- Adiabatic invariants: `perturbation.jl`

### ETH — 402-0800 Analytische Mechanik
- Lagrange formalism: `generalized.jl`, `euler_lagrange.jl`
- Hamilton formalism: `legendre.jl`
- Symmetrien: `noether.jl`

### 东京大学 — 解析力学
- ラグランジュ形式: full module coverage
- ネーター定理: `noether.jl`
- 剛体力学: `rigid_body.jl`

## 教材对标

| 教材 | 覆盖章节 | 实现映射 |
|------|---------|---------|
| Goldstein, Poole & Safko (2002) | Ch.1-8, 11, 13 | 全部 14 个源文件 |
| Landau & Lifshitz Vol.1 (1987) | Ch.1-7 | generalized, euler_lagrange, action, central_force, small_oscillations, rigid_body, legendre |
| Arnold, Mathematical Methods (1989) | Ch.2-5, 10 | variational_calculus, perturbation, legendre |
| Marsden & Ratiu, Mechanics & Symmetry | Ch.1-4 | noether, legendre, variational_calculus |
| José & Saletan, Classical Dynamics | Ch.2-6 | 全部模块 |
