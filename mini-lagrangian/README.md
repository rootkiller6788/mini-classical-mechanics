# mini-lagrangian — 拉格朗日力学 (Julia)

> 参考 Goldstein, Landau & Lifshitz Vol.1, MIT 8.012, Stanford PHYSICS 370
> 从最小作用量原理推导运动方程，实现完整九层知识覆盖。

## Module Status: COMPLETE ✅

| Level | Status | Description |
|-------|--------|-------------|
| L1 Definitions | Complete | 10+ struct types covering all core entities |
| L2 Core Concepts | Complete | Least action, EL equations, symmetry, gauge invariance |
| L3 Math Structures | Complete | SO(3), symplectic geometry, Poisson brackets, Christoffel |
| L4 Fundamental Laws | Complete | EL eq, Hamilton eq, Noether theorem, Euler eq, Kepler laws |
| L5 Computational Methods | Complete | RK4/RK45, Verlet, symplectic integrators, Floquet, averaging |
| L6 Canonical Systems | Complete | 11 canonical systems solved with examples/ |
| L7 Applications | Partial+ | EM particle tracking, coupled oscillators, GPS relativity |
| L8 Advanced Topics | Complete | Relativistic L, gauge theory, adiabatic invariants, KAM |
| L9 Research Frontiers | Partial | KAM nondegeneracy, Berry phase, AB effect, Lyapunov chaos |

**Score: 16/18** (L1-L6 Complete + L7 Partial+ + L8 Complete + L9 Partial)
**准入: include/ + src/ = 3346 lines (>= 3000)**

## 模块对照

| 子模块 | 内容 | Goldstein | 实现 |
|--------|------|-----------|------|
| Generalized Coords | 广义坐标、坐标变换、质量矩阵、Christoffel 符号 | Ch.1 | `generalized.jl` |
| Euler-Lagrange Eq | EL方程数值求解器、耗散系统、显含时系统 | Ch.2 | `euler_lagrange.jl` |
| Action Principle | 作用量泛函、离散变分积分器、高阶辛积分器 | Ch.2 | `action.jl` |
| Constraints | 拉格朗日乘子法、完整/非完整约束 | Ch.2 | `constraints.jl` |
| Legendre Transform | L→H 变换、Poisson 括号、Routhian、Fenchel 共轭 | Ch.8 | `legendre.jl` |
| Small Oscillations | 简正模、受迫振动、频率响应、Rayleigh 阻尼 | Ch.6 | `small_oscillations.jl` |
| Noether's Theorem | 对称性与守恒量、场论 Noether、EM 守恒流 | Ch.13 | `noether.jl` |
| ODE Integrators | RK4/RK45、Verlet、Leapfrog、辛Euler、DOPRI5 | — | `integrators.jl` |
| Central Force | Kepler轨道、Binet公式、Rutherford散射、Bertrand定理 | Ch.3 | `central_force.jl` |
| Rigid Body | Euler角、SO(3)、转动惯量、Euler方程、重力陀螺 | Ch.4-5 | `rigid_body.jl` |
| Variational Calculus | Gateaux导数、Maupertuis、Jacobi度规、Beltrami恒等式 | Ch.2 | `variational_calculus.jl` |
| Relativistic | 相对论粒子L、四动量、Lorentz变换、双生子佯谬 | Ch.7 | `relativistic.jl` |
| Gauge Theory | 规范变换、场强张量、Berry相位、Dirac磁单极子 | — | `gauge_theory.jl` |
| Perturbation Theory | 绝热不变量、平均法、Poincaré截面、KAM理论 | Ch.11 | `perturbation.jl` |

## API

### 广义坐标 (L1, L3)

| 函数 | 描述 |
|------|------|
| `mass_matrix(masses, jacobian)` | 构建质量矩阵 M(q): T=½q̇ᵀMq̇ |
| `generalized_momenta(M, qdot)` | 广义动量 p = M·q̇ |
| `polar_to_cartesian(r,θ)` | 极坐标→笛卡尔 |
| `spherical_to_cartesian(r,θ,φ)` | 球坐标→笛卡尔 |
| `cylindrical_to_cartesian(rho,φ,z)` | 柱坐标→笛卡尔 |
| `elliptic_to_cartesian(c,μ,ν)` | 椭圆坐标→笛卡尔 |
| `christoffel_symbols(M_func, q)` | Christoffel 符号 Γⁱ_{jk} |
| `geodesic_acceleration(M_func, q, qdot)` | 测地线加速度 |

### EL 求解器 (L4)

| 函数 | 描述 |
|------|------|
| `EulerLagrangeSystem(n, L, ∂L/∂q, ∂L/∂q̇, M)` | 定义 EL 系统 |
| `el_to_first_order(sys)` | EL → 一阶 ODE |
| `standard_el_to_ode(sys)` | 标准 L=T-U 的 ODE |
| `rayleigh_dissipation_lagrange(R, ...)` | 含耗散函数的 EL |
| `time_dependent_el_to_ode(...)` | 显含时拉格朗日量 |
| `energy_change_rate(L, q, qdot, t)` | dE/dt = -∂L/∂t |

### 变分原理 (L2)

| 函数 | 描述 |
|------|------|
| `action(L, q, qdot, t_span)` | 计算作用量 S[q] |
| `verify_least_action(L, q_true, q_var, ...)` | 验证最小作用量原理 |
| `first_variation(L, ...)` | 一阶变分 δS |
| `second_variation_sign(L, ...)` | 二阶变分符号 (:minimum/:saddle) |
| `variational_midpoint_step(...)` | 离散变分积分器单步 |
| `gauss_legendre_4th_order_step(...)` | 四阶 Gauss-Legendre 辛积分器 |

### 约束系统 (L4)

| 函数 | 描述 |
|------|------|
| `holonomic_constraint(name, f, ∇f)` | 完整约束 f(q,t)=0 |
| `nonholonomic_constraint(name, f, ∇_q, ∇_qdot)` | 非完整约束 |
| `constrained_el_to_ode(sys)` | 含约束的扩展 ODE |
| `spherical_constraint(R)` | 球面约束 |
| `rigid_rod_constraint(L)` | 刚性杆约束 |
| `rolling_disk_constraints(R)` | 纯滚动约束 |

### Legendre 变换 (L3)

| 函数 | 描述 |
|------|------|
| `qdot_to_momentum(lt, q, qdot)` | p = M·q̇ |
| `momentum_to_qdot(lt, q, p)` | q̇ = M⁻¹·p |
| `hamiltonian(lt, q, p)` | H = T + U |
| `hamiltons_equations(lt)` | Hamilton 正则方程 |
| `poisson_bracket(f, g, q, p)` | {f, g} 数值计算 |
| `routhian(L, ..., cyclic)` | Routhian 混合表述 |

### Noether 定理 (L4)

| 函数 | 描述 |
|------|------|
| `energy_from_lagrangian(L, q, qdot, grad)` | E = p·q̇ - L |
| `noether_charge(grad_L_qdot, Q)` | Noether 荷 = p·Q |
| `translation_symmetry_Q(N, dir)` | 空间平移生成元 |
| `rotation_symmetry_Q_z(positions)` | 旋转生成元 |
| `field_noether_current(...)` | 场论守恒流 (T^μ_ν) |
| `poynting_vector(E, B)` | Poynting 矢量 |
| `maxwell_stress_tensor(E, B)` | Maxwell 应力张量 |

### 小振动 (L6)

| 函数 | 描述 |
|------|------|
| `find_equilibrium(grad_U, hess_U, q0)` | Newton-Raphson 找平衡点 |
| `solve_normal_modes(M, K, q_eq)` | 广义特征值问题 |
| `numerical_hessian(U, q)` | 数值 Hessian |
| `forced_response(M, K, C, F0, ω)` | 受迫振动频率响应 |
| `frequency_sweep(M, K, C, F0, ω_range)` | 扫频分析 |

### 中心力问题 (L6)

| 函数 | 描述 |
|------|------|
| `effective_potential(sys, r, L)` | U_eff = U + L²/(2mr²) |
| `kepler_orbit(m, k, E, L)` | 从 (E,L) 构造 Kepler 轨道 |
| `solve_kepler_equation(e, M)` | Kepler 方程 Newton 求解 |
| `laplace_runge_lenz(m, k, r, v)` | LRL 矢量 |
| `scattering_angle(b, E, k, m)` | Rutherford 散射角 |
| `bertrand_check(...)` | Bertrand 定理验证 |

### 刚体动力学 (L6)

| 函数 | 描述 |
|------|------|
| `rotation_matrix_euler(angles)` | Euler 角 → R ∈ SO(3) |
| `inertia_tensor_point(masses, pos)` | 离散转动惯量张量 |
| `euler_equations_ode(I_body)` | Euler 方程 ODE |
| `torque_free_precession(I, ω0, ...)` | 无力矩进动模拟 |
| `heavy_top_lagrangian(I1, I3, ...)` | 重力陀螺拉格朗日量 |

### ODE 积分器 (L5)

| 函数 | 描述 |
|------|------|
| `rk4_step(f, t, y, h)` | 经典 RK4 (O(h⁴)) |
| `rk4_integrate(f, y0, span, h)` | RK4 全区间积分 |
| `verlet_integrate(a, q0, v0, ...)` | Störmer-Verlet (辛) |
| `leapfrog_integrate(a, q0, v0, ...)` | Leapfrog (二阶辛) |
| `rk45_adaptive_step(f, t, y, h)` | DOPRI5 自适应 |

### 变分法 (L2, L3)

| 函数 | 描述 |
|------|------|
| `gateaux_derivative(F, y, η)` | Gateaux 方向导数 |
| `maupertuis_action(E, U, path)` | Maupertuis 约化作用量 |
| `jacobi_metric(E, U, M, q)` | Jacobi 度规张量 |
| `beltrami_identity(L, dL_dydot, y, ydot)` | Beltrami 恒等式 |
| `isoperimetric_el(...)` | 等周问题 EL 方程 |

### 相对论性拉格朗日量 (L8)

| 函数 | 描述 |
|------|------|
| `relativistic_particle_L(m, v²)` | L = -mc²√(1-v²/c²) |
| `relativistic_energy(m, v²)` | E = γmc² |
| `relativistic_momentum(m, v)` | p = γmv |
| `four_velocity(v)` | 四速度 u^μ |
| `lorentz_boost_matrix(v)` | Lorentz boost 矩阵 |
| `twin_paradox(D, v)` | 双生子佯谬计算 |

### 规范理论 (L8)

| 函数 | 描述 |
|------|------|
| `gauge_transform_A(A, χ, r, t)` | A → A + ∇χ |
| `field_strength_tensor(A, r, t)` | F_{μν} = ∂_μ A_ν - ∂_ν A_μ |
| `magnetic_field_from_A(A, r, t)` | B = ∇ × A |
| `dirac_monopole_A(g, r, θ, φ)` | Dirac 磁单极子矢势 |
| `berry_connection(state, R)` | Berry 联络 A_n(R) |
| `aharonov_bohm_phase(q, A, path)` | AB 相位 Δφ |

### 摄动理论 (L8, L9)

| 函数 | 描述 |
|------|------|
| `adiabatic_invariant_1d(m, E, U, ...)` | 作用量 I = (1/π)∫p dq |
| `period_from_action(m, E, U, ...)` | T = 2π dI/dE |
| `averaging_lagrangian(ω, f, A, φ)` | 平均法缓变方程 |
| `poincare_section(f, y0, ...)` | Poincaré 截面 |
| `maximal_lyapunov_exponent(f, J, y0, ...)` | 最大 Lyapunov 指数 |
| `floquet_multipliers(A, T, n)` | Floquet 乘子 |
| `mathieu_stability(δ, ε)` | Mathieu 方程稳定性 |
| `kam_nondegeneracy_check(H0, I)` | KAM 非退化条件 |

## Examples

| Example | 内容 | 关键概念 |
|---------|------|---------|
| `pendulum_lagrangian.jl` | 单摆：Lagrange vs Newton | 能量函数 = Hamiltonian |
| `double_pendulum_lagrangian.jl` | 双摆 | 质量矩阵的非对角耦合 |
| `bead_rotating_hoop.jl` | 旋转环上珠子 | 有效势、分岔、平衡点稳定性 |
| `charged_particle_em.jl` | 电磁场中带电粒子 | 速度依赖势、正则动量 ≠ 机械动量 |
| `atwood_lagrangian.jl` | Atwood 机 | Lagrange vs Newton 张力分析 |
| `brachistochrone.jl` | 最速降线 | 变分法、Beltrami 恒等式 |
| `coupled_oscillators.jl` | 耦合谐振子 | 简正模、对称性分析 |
| `sliding_block_wedge.jl` | 滑块-楔子系统 | 循环坐标、动量守恒 |

## 九校课程映射

| 学校 | 关键课程 | 覆盖章节 |
|------|---------|---------|
| **MIT** | 8.012 Mechanics | Lectures 20-24, 26-27 (Lagrangian, oscillations, rigid body) |
| **Stanford** | PHYSICS 370 CM | Goldstein Ch.1-8, 11, 13 |
| **Berkeley** | PHYS 242 CM | Lagrangian + Hamiltonian formulations |
| **Caltech** | Ph 106 CM | Classical mechanics complete |
| **Princeton** | PHY 505 CM | Advanced classical mechanics |
| **Cambridge** | Part II Theo Phys | Lagrangian & Hamiltonian mechanics |
| **Oxford** | CMT Graduate | Theoretical mechanics core |
| **ETH** | 402-0800 CM | Analytische Mechanik |
| **东京大学** | 解析力学 | Lagrangian formalism |

## 知识覆盖摘要

| 知识点数 | L1 | L2 | L3 | L4 | L5 | L6 | L7 | L8 | L9 |
|---------|----|----|----|----|----|----|----|----|-----|
| 条目数 | 12 | 8 | 10 | 10 | 9 | 11 | 3 | 8 | 4 |

## 目录结构

```
mini-lagrangian/
├── Project.toml
├── Makefile
├── README.md
├── include/                    # API 接口声明 (4 files)
│   ├── lagrangian_types.jl
│   ├── lagrangian_api.jl
│   ├── lagrangian_api_2.jl
│   └── lagrangian_api_3.jl
├── src/                        # 核心实现 (14 files)
│   ├── Lagrangian.jl
│   ├── generalized.jl
│   ├── euler_lagrange.jl
│   ├── action.jl
│   ├── constraints.jl
│   ├── legendre.jl
│   ├── small_oscillations.jl
│   ├── noether.jl
│   ├── integrators.jl
│   ├── central_force.jl
│   ├── rigid_body.jl
│   ├── variational_calculus.jl
│   ├── relativistic.jl
│   ├── gauge_theory.jl
│   └── perturbation.jl
├── examples/                   # 8 个端到端可运行示例
├── tests/
│   └── runtests.jl
└── docs/
    ├── knowledge-graph.md
    ├── coverage-report.md
    ├── gap-report.md
    ├── course-tree.md
    ├── api-map.md
    ├── university-coverage.md
    └── lagrangian-cheatsheet.md
```

## 运行

```bash
# 运行所有测试
make test

# 运行所有示例
make examples

# 安全检查
make check

# 验证 COMPLETE 状态
make verify
```
