# mini-variational-principle — 变分原理 (C)

> 参考 Gelfand & Fomin "Calculus of Variations", Goldstein Ch.2, Courant & Hilbert
> 泛函分析、Euler-Lagrange方程、约束变分、二阶变分理论、Hamilton原理、
> 最优控制、场论变分、数值变分方法。

## Module Status: COMPLETE ✅

- **include/ + src/ 总行数**: 3264 (≥ 3000 ✅)
- **L1 Definitions**: Complete — 8 struct typedefs, FunctionalContext, FunctionalDiscretization
- **L2 Core Concepts**: Complete — 泛函求值(3种积分方法), 范数, 增量, 能量分解
- **L3 Math Structures**: Complete — Gateaux/Frechet导数, Jacobi形式, 凸性分析
- **L4 Fundamental Laws**: Complete — Euler-Lagrange方程, Beltrami恒等式, Legendre/Jacobi/Weierstrass条件
- **L5 Algorithms**: Complete — Shooting BVP, FEM, Ritz, Galerkin, 谱方法, Newton-Kantorovich
- **L6 Canonical Systems**: Complete — 最速降线, 悬链线, 谐振子, Dido问题, Ginzburg-Landau, Allen-Cahn
- **L7 Applications**: Partial+ — LQR控制, Bellman动态规划, 障碍问题, 相场模型
- **L8 Advanced Topics**: Partial — 多辛结构, HJB方程, 离散Hamilton原理, 拓扑荷密度
- **L9 Research Frontiers**: Partial — 文档中标注, 不做强制实现

## 核心定义 (L1)

| 定义 | 对应类型 | 文件 |
|------|---------|------|
| Lagrangian密度 | `LagrangianDensity` | `include/functional.h` |
| 泛函离散表示 | `FunctionalDiscretization` | `include/functional.h` |
| 泛函上下文 | `FunctionalContext` | `include/functional.h` |
| Legendre/Jacobi条件 | `SufficientConditions` | `include/second_variation.h` |

## 核心定理 (L4)

| 定理 | 公式 | 实现 |
|------|------|------|
| **Euler-Lagrange方程** | ∂F/∂y - d/dx(∂F/∂y') = 0 | `solve_el_bvp()`, `el_residual_norm()` |
| **Beltrami恒等式** | F - y'∂F/∂y' = const (∂F/∂x=0) | `beltrami_identity()`, `verify_beltrami()` |
| **Legendre必要条件** | ∂²F/∂y'² ≥ 0 | `legendre_condition()`, `strong_legendre()` |
| **Jacobi必要条件** | (a,b)内无共轭点 | `find_conjugate_points()`, `jacobi_condition()` |
| **Weierstrass条件** | E(x,y,pₒ,p) ≥ 0 ∀p | `weierstrass_excess()`, `weierstrass_condition()` |
| **Pontryagin极大值原理** | u* = argminᵤ H(x,u,p,t) | `pontryagin_hamiltonian()`, `switching_function()` |
| **LQR Riccati方程** | AᵀP + PA - PBR⁻¹BᵀP + Q = 0 | `lqr_controller()` |

## 核心算法 (L5)

| 算法 | 复杂度 | 文件 |
|------|--------|------|
| 打靶法(Bisection)求解EL BVP | O(n_grid × log(1/ε)) | `euler_lagrange.c` |
| RK4打靶法 | O(n_grid × log(1/ε)) | `euler_lagrange.c` |
| FEM求解线性变分问题 | O(n³) | `euler_lagrange.c` |
| Ritz直接法 | O(n_basis² × n_quad) | `second_variation.c` |
| Galerkin加权残差法 | O(n_basis² × n_quad) | `second_variation.c` |
| 谱Galerkin (Fourier基) | O(n_modes² × n_quad) | `numerical_variational.c` |
| 障碍问题投影算法 | O(n_iter × n) | `numerical_variational.c` |
| Bellman值迭代 | O(n_steps × nx × nu) | `optimal_control.c` |
| Newton-Kantorovich线性化 | O(max_iter × n³) | `numerical_variational.c` |
| Gauss-Legendre 4点求积 | O(n) | `functional.c` |

## 经典问题 (L6)

| 问题 | Lagrangian | 解 | API |
|------|-----------|-----|-----|
| 最速降线 | √(1+y'²)/√(2gy) | 摆线 (cycloid) | `brachistochrone_lagrangian()` |
| 悬链线 | y√(1+y'²) | y = a·cosh(x/a) | `catenary_lagrangian()` |
| 谐振子 | ½my'² - ½ky² | y = A cos(ωt+φ) | `harmonic_oscillator_lagrangian()` |
| Dido问题 | 等周约束 | 圆 R = L/(2π) | `dido_solution()` |
| 测地线 | √(Σẋ_i²) | 直线/大圆 | `geodesic_integrand()` |
| 极小曲面 | √(1+u_x²+u_y²) | 悬链面/螺旋面 | `minimal_surface_integrand()` |
| Ginzburg-Landau | a|ψ|²+½b|ψ|⁴+c|∇ψ|² | 超导序参量 | `solve_ginzburg_landau()` |
| Allen-Cahn | ε²|∇φ|²+W(φ) | 相场界面 | `allen_cahn_rhs()` |

## 模块对照

| 子模块 | 内容 | C实现 |
|--------|------|-------|
| Functional | 泛函求值、Gateaux/Frechet导数、Beltrami恒等式 | `functional.h/c` |
| Euler-Lagrange | BVP求解、多变量EL、高阶导数、场论EL | `euler_lagrange.h/c` |
| Constrained | 等周问题、Lagrange乘子、可变端点、Bolza | `constrained.h/c` |
| Second Variation | Legendre/Jacobi/Weierstrass、Ritz/Galerkin | `second_variation.h/c` |
| Hamilton Principle | 相空间作用量、辛结构、离散变分 | `hamilton_principle.h/c` |
| Optimal Control | Pontryagin、LQR、Bellman DP、MPC | `optimal_control.h/c` |
| Field Theory | Ginzburg-Landau、相场、弹性、电磁 | `field_theory.h/c` |
| Numerical | FEM、谱方法、障碍问题、自适应步长 | `numerical_variational.h/c` |

## 九校课程映射

| 学校 | 课程 | 对应模块 |
|------|------|---------|
| MIT | 8.012 Classical Mechanics | Hamilton Principle, Euler-Lagrange |
| Stanford | PHYSICS 230 Graduate Mechanics | Phase space, Symplectic |
| Berkeley | PHYS 242 Theoretical Mechanics | Gateaux/Frechet, Convexity |
| Caltech | Ph 106 Classical Mechanics | Canonical systems, Variational |
| Princeton | PHY 505 Classical Mechanics | Constrained variations |
| Cambridge | Part II Theoretical Physics | Field theory EL, Beltrami |
| Oxford | CMT Graduate | Legendre/Jacobi theory |
| ETH | 402-0800 Classical Mechanics | Numerical variational methods |
| MIT | 16.323 Optimal Control | Pontryagin, LQR, HJB |

## 编译与测试

```bash
# 编译全部
make

# 运行测试 (24 tests)
make test

# 查看行数
make lines

# 清理
make clean
```

## 即将到来的主题 (L7-L9)

- **L7 应用**: 蛋白质折叠自由能、气候模型变分同化
- **L8 进阶**: Γ-收敛、多尺度均匀化、随机变分法
- **L9 前沿**: 量子变分算法、机器学习中的信息几何

