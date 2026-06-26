# mini-classical-mechanics — 统一对齐工作流

> **目标**: 所有 8 个子模块的 `src/` 体量和知识颗粒度向 `mini-newtonian` 看齐
> **基准**: mini-newtonian = 8 src files, 1204 lines, 7 examples, 3 docs, 2 demos, 38 tests
> **最后更新**: 2026-06-08 (Phase 1+2 全部完成)

---

## 基准标准: mini-newtonian ✓ GOLD STANDARD

| 维度 | 指标 | 说明 |
|------|------|------|
| src 文件数 | 8 个模块化文件 | 每文件单一职责, 100-250 行 |
| src 总行数 | 1204 lines | 涵盖完整知识体系 |
| 类型系统 | Vec3, State, NBodyState, Trajectory, ODEProblem | 不可变栈分配 |
| 数值方法 | Euler/RK2/RK4/RK45/Verlet/Leapfrog | 手写全部积分器 |
| 物理模型 | 7 种力定律 + 能量/动量/碰撞/约束 | 公式→代码一一对应 |
| examples | 7 个独立脚本 | 每个含解析解验证 |
| docs | 3 个文档 | course-alignment + integrator-analysis + cheatsheet |
| demos | 2 个交互演示 | mini-solar-system + mini-chaotic-pendulum |
| tests | 38 项 | 覆盖所有导出函数 |
| README | 200 行 API 参考表 | 每个函数一行说明 |
| 加速层 | C++17 extern "C" | N>100 热路径 |
| Makefile | ✓ | make test / examples / cpp-build |

---

## 模块状态矩阵

| # | 模块 | src 文件数 | src 总行数 | 颗粒度评分 | examples | docs | demos | 状态 |
|---|------|-----------|-----------|-----------|----------|------|------|------|
| 0 | mini-newtonian | 8 | 1204 | ⭐⭐⭐⭐⭐ | 7 | 3 | 2 | **基准** |
| 1 | mini-lagrangian | 8 | 901 | ⭐⭐⭐⭐ | 8 | 1 | 1 | 接近 |
| 2 | mini-hamiltonian | 10 | 798 | ⭐⭐⭐⭐ | 7 | 1 | 1 | 接近 |
| 3 | mini-chaos | 1 | 92 | ⭐ | 1 | 0 | 0 | **差距大** |
| 4 | mini-celestial-mechanics | 1 | 70 | ⭐ | 1 | 0 | 0 | **差距大** |
| 5 | mini-rigid-body | 1 | 49 | ⭐ | 1 | 0 | 0 | **差距大** |
| 6 | mini-variational-principle | 5 | 124 | ⭐⭐ | 2 | 0 | 0 | **差距大** |
| 7 | mini-continuum | 1 | 32 | ⭐ | 2 | 0 | 0 | **差距大** |

---

## Phase 1: 单文件巨石拆分 (mini-rigid-body, celestial, continuum, chaos)

---

### Step 1 — mini-rigid-body: 49行 → 目标 ~800行

**当前**: 1 文件 `RigidBody.jl` (49 lines)，所有代码堆在一个 module 里。

**目标拆分**:
```
mini-rigid-body/src/
├── RigidBody.jl          # 主模块入口 (~40 lines)
├── types.jl              # InertiaTensor, EulerAngles, RigidBodyState (~100 lines)
├── inertia.jl            # 惯性张量计算、主轴变换、平行轴定理 (~120 lines)
├── euler_equations.jl    # Euler 方程、力矩驱动、数值积分 (~150 lines)
├── kinematics.jl         # Euler 角↔旋转矩阵↔角速度 互转 (~120 lines)
├── energy_momentum.jl    # 转动动能、角动量、进动稳定性 (~100 lines)
└── tops.jl               # 陀螺: Lagrange top, sleeping top, precession+nutation (~150 lines)
```

**需要补充的知识颗粒**:
- [ ] `types.jl`: InertiaTensor (已有), 新增 EulerAngles struct, RigidBodyState, Quaternion(可选)
- [ ] `inertia.jl`: 已有基础, 新增: 常用几何体惯性张量解析公式(球/圆柱/长方体/细杆), 惯性椭球
- [ ] `euler_equations.jl`: 已有基础, 新增: RK4 积分, torque-free 解析解(Jacobi elliptic), Poinsot 几何描述
- [ ] `kinematics.jl`: 已有 euler_to_rotation, 新增: rotation_to_euler, 角速度↔Euler角速率互转矩阵, 四元数(可选)
- [ ] `energy_momentum.jl`: 已有基础, 新增: 网球拍定理(Dzhanibekov effect)详细分析, 动能椭球
- [ ] `tops.jl` (新): Lagrange top 运动方程, 章动, sleeping top 稳定性, 快速陀螺近似

**新增 examples** (当前 1 → 目标 5):
- [ ] `free_precession.jl` (已有, 需增强)
- [ ] `tennis_racket_theorem.jl` — 中间轴不稳定性演示
- [ ] `spinning_top.jl` — 陀螺进动+章动
- [ ] `tippe_top.jl` — 翻转陀螺
- [ ] `gyroscope.jl` — 陀螺仪原理

---

### Step 2 — mini-celestial-mechanics: 70行 → 目标 ~700行

**当前**: 1 文件 `CelestialMechanics.jl` (70 lines)

**目标拆分**:
```
mini-celestial-mechanics/src/
├── CelestialMechanics.jl    # 主模块入口 (~30 lines)
├── types.jl                 # OrbitalElements, StateVector, KeplerOrbit (~80 lines)
├── kepler.jl                # Kepler 方程求解、轨道根数↔状态向量 (~150 lines)
├── two_body.jl              # 二体问题: vis-viva, 轨道周期, 位置/速度函数 (~120 lines)
├── perturbations.jl         # J2 摄动、大气阻力、太阳光压 (~120 lines)
├── three_body.jl            # CRTBP, Lagrange 点(已有), Jacobi 常数, 零速度面 (~150 lines)
└── mission.jl               # Hohmann 转移, 双椭圆转移, 引力辅助, patched conics (~150 lines)
```

**需要补充的知识颗粒**:
- [ ] `types.jl`: OrbitalElements struct (6 个经典轨道根数), 轨道类型判别
- [ ] `kepler.jl`: 已有 solve_kepler, elements↔state. 新增: 双曲线/抛物线 Kepler 方程, 普适变量公式
- [ ] `two_body.jl`: 已有 vis_viva, orbital_period. 新增: 轨道位置/速度的时间函数, 飞行时间, 发射窗口
- [ ] `perturbations.jl` (新): J2 摄动(进动), 大气阻力衰降, 太阳光压, 长期/周期摄动分离
- [ ] `three_body.jl`: 已有 lagrange_points. 新增: Jacobi 积分, 零速度面, Hill 区域, 晕轨道概念
- [ ] `mission.jl` (新): Hohmann Δv, 双椭圆转移, 平面变化, 引力辅助(飞掠), Patched conics

**新增 examples** (当前 1 → 目标 5):
- [ ] `orbital_elements.jl` (已有, 需增强)
- [ ] `hohmann_transfer.jl` — 霍曼转移轨道
- [ ] `j2_perturbation.jl` — J2 摄动下轨道进动
- [ ] `lagrange_points.jl` — 5 个 Lagrange 点 + 零速度面
- [ ] `interplanetary.jl` — 行星际转移 (patched conics)

---

### Step 3 — mini-continuum: 32行 → 目标 ~600行

**当前**: 1 文件 `Continuum.jl` (32 lines)，全部一行式实现

**目标拆分**:
```
mini-continuum/src/
├── Continuum.jl            # 主模块入口 (~30 lines)
├── types.jl                # StressTensor, StrainTensor, DisplacementField (~80 lines)
├── elasticity.jl           # Hooke 定律、Lame 常数、弹性常数互转 (~120 lines)
├── waves.jl                # P/S/Rayleigh/Love 波、波速、色散 (~120 lines)
├── beams.jl                # Euler-Bernoulli 梁、Timoshenko 梁、屈曲 (~120 lines)
├── plates.jl               # Kirchhoff 板、膜、振动 (~100 lines)
└── failure.jl              # von Mises, Tresca, Mohr-Coulomb, 疲劳 (~100 lines)
```

**需要补充的知识颗粒**:
- [ ] `types.jl`: 完整的 StressTensor, StrainTensor (3x3), 主应力/主应变
- [ ] `elasticity.jl`: 已有基础, 新增: 弹性常数互转表(E,ν,G,K,λ,μ), 各向异性 Hooke, 平面应力/平面应变
- [ ] `waves.jl`: 已有 p/s/rayleigh 波速. 新增: 1D 波方程完整求解器(已有基础), 反射/透射, Love 波, 色散曲线
- [ ] `beams.jl`: 已有 beam_deflection 一行式. 新增: Euler-Bernoulli 完整求解, Timoshenko 梁, 屈曲临界载荷(Euler buckling)
- [ ] `plates.jl` (新): Kirchhoff 板方程, 圆板轴对称, 膜振动, 简正模
- [ ] `failure.jl`: 已有 von_mises, tresca. 新增: Mohr-Coulomb, 最大主应力准则, S-N 曲线(疲劳), 安全系数

**新增 examples** (当前 2 → 目标 5):
- [ ] `wave_and_beam.jl` (已有, 需增强)
- [ ] `buckling.jl` — Euler 屈曲
- [ ] `seismic_waves.jl` — 地震波速 + Snell 定律
- [ ] `plate_vibration.jl` — 板振动简正模
- [ ] `stress_concentration.jl` — 应力集中 + 屈服准则对比

---

### Step 4 — mini-chaos: 92行 → 目标 ~600行

**当前**: 1 文件 `Chaos.jl` (92 lines)，所有函数定义在单文件里

**目标拆分**:
```
mini-chaos/src/
├── Chaos.jl                # 主模块入口 (~30 lines)
├── types.jl                # DynamicalSystem, BifurcationData (~60 lines)
├── flows.jl                # Lorenz, Rössler, Chua, Sprott 系统 (~120 lines)
├── maps.jl                 # Logistic, Hénon, Standard map, Circle map (~120 lines)
├── lyapunov.jl             # Lyapunov 指数(多种算法)、Lyapunov 谱 (~120 lines)
├── bifurcation.jl          # 分岔图, Feigenbaum 常数, 周期倍增检测 (~100 lines)
└── fractals.jl             # Mandelbrot, Julia 集, 分形维数 (~80 lines)
```

**需要补充的知识颗粒**:
- [ ] `types.jl`: DynamicalSystem struct, 参数化系统定义, 分岔数据类型
- [ ] `flows.jl`: 已有 Lorenz, Rössler. 新增: Chua 电路, Sprott 系统, Duffing 吸引子, 强迫振子
- [ ] `maps.jl`: 已有 Logistic, Hénon. 新增: Standard map (Chirikov-Taylor), Circle map, Ikeda map
- [ ] `lyapunov.jl`: 已有基本 Lyapunov. 新增: Benettin 算法求全谱, 条件 Lyapunov 指数, 同步判定
- [ ] `bifurcation.jl`: 已有 logistic_bifurcation. 新增: 通用分岔图生成器, 周期检测算法, 分岔类型判定(saddle-node/pitchfork/Hopf)
- [ ] `fractals.jl`: 已有 mandelbrot_iter. 新增: Julia 集, Newton 分形, 盒计数维数

**新增 examples** (当前 1 → 目标 5):
- [ ] `lorenz_logistic.jl` (已有, 需增强)
- [ ] `bifurcation_analysis.jl` — 分岔图 + Feigenbaum 收敛
- [ ] `lyapunov_spectrum.jl` — Lorenz 全 Lyapunov 谱
- [ ] `fractal_gallery.jl` — Mandelbrot + Julia 集
- [ ] `synchronization.jl` — 混沌同步 (Pecora-Carroll)

---

## Phase 2: 薄文件充实

### Step 5 — mini-variational-principle: 124行 → 目标 ~500行

**当前**: 5 文件但都很薄 (functional 26L, euler_lagrange 29L, constrained 26L, second_variation 36L)

**需要补充的知识颗粒**:
- [ ] `functional.jl` (26→~120): 新增: 多种泛函类型, 弱/强极值判别, 泛函空间, 直接法
- [ ] `euler_lagrange.jl` (29→~120): 新增: 多变量 EL 方程 (场论), 高阶导数 EL, 参数化曲线
- [ ] `constrained.jl` (26→~100): 新增: 多点边值, 可变端点, 横截条件增强
- [ ] `second_variation.jl` (36→~120): 已有基础, 新增: Weierstrass 必要条件, Hilbert 不变积分, 场论极值

**新增 examples** (当前 2 → 目标 5):
- [ ] `brachistochrone.jl` (已有, 需增强)
- [ ] `catenary.jl` (已有, 需增强)
- [ ] `geodesic.jl` — 测地线 (球面/柱面)
- [ ] `minimal_surface.jl` — 极小曲面 (肥皂膜)
- [ ] `isoperimetric.jl` — 等周问题 (Dido 问题)

### Step 6 — mini-lagrangian: 901行 → 目标 ~1200行

**当前**: 8 文件, 结构合理但 action.jl(96L), noether.jl(72L), legendre.jl(132L) 偏薄
- [ ] `action.jl`: 新增离散变分积分器增强, 辛性质证明, 作用量泛函二阶变分
- [ ] `noether.jl`: 新增更多对称性生成元(时间平移、boost), 场论 Noether 定理预处理
- [ ] `small_oscillations.jl`: 新增受迫小振动, 耗散系统简正模
- [ ] `legendre.jl`: 新增 Routhian 消去循环坐标

### Step 7 — mini-hamiltonian: 798行 → 目标 ~1100行

**当前**: 10 文件, poisson(65L), action_angle(65L), liouville(60L) 偏薄
- [ ] `poisson.jl`: 新增 Poisson 括号几何解释, 辛矩阵形式, 更多运动常数示例
- [ ] `action_angle.jl`: 新增更多可积系统的作用量(无限深势阱、Poschl-Teller), 绝热不变量
- [ ] `liouville.jl`: 新增系综演化可视化, 混合与遍历初步
- [ ] `hamilton_jacobi.jl`: 新增更多可分离系统, 作用量作为 H-J 解

---

## Phase 3: 统一标准化 (所有模块)

### README 标准化要求
每个模块 README 必须包含:
- [ ] 参考教材标注
- [ ] 完整 API Reference 表格 (每个导出函数一行)
- [ ] 模块对照表 (子模块→章节→文件)
- [ ] Examples 列表 + 验证项
- [ ] 目录结构树
- [ ] 运行命令
- [ ] 与其他模块的关系说明

### docs/ 补齐
- [ ] `rigid-body-cheatsheet.md`
- [ ] `celestial-cheatsheet.md`
- [ ] `continuum-cheatsheet.md`
- [ ] `chaos-cheatsheet.md`
- [ ] `variational-cheatsheet.md`

### examples/ 补齐
| 模块 | 当前 | 目标 | 缺口 |
|------|------|------|------|
| mini-newtonian | 7 | 7 | ✓ |
| mini-lagrangian | 8 | 8 | ✓ |
| mini-hamiltonian | 7 | 7 | ✓ |
| mini-rigid-body | 1 | 5 | +4 |
| mini-celestial-mechanics | 1 | 5 | +4 |
| mini-continuum | 2 | 5 | +3 |
| mini-variational-principle | 2 | 5 | +3 |
| mini-chaos | 1 | 5 | +4 |

### tests/ 补齐
目标: 每个模块 ≥20 项测试，覆盖所有导出函数

### demos/ 补齐
目标: 每个模块 ≥1 个 demo

---

## 执行顺序总览

```
Phase 1 (优先): 单文件巨石拆分
├── Step 1: mini-rigid-body     (49L→800L, 最薄)
├── Step 2: mini-celestial      (70L→700L)
├── Step 3: mini-continuum      (32L→600L)
├── Step 4: mini-chaos          (92L→600L)

Phase 2 (其次): 薄文件充实
├── Step 5: mini-variational    (124L→500L)
├── Step 6: mini-lagrangian     (901L→1200L)
└── Step 7: mini-hamiltonian    (798L→1100L)

Phase 3 (最后): 统一标准化
├── Step 8: README 全部标准化
├── Step 9: docs/ cheatsheet 补齐
├── Step 10: examples 补齐至 5+
├── Step 11: tests 补齐至 20+
└── Step 12: demos 补齐至 1+
```

---

## 进展追踪

| Step | 模块 | 状态 | 开始 | 完成 |
|------|------|------|------|------|
| 1 | mini-rigid-body | ✅ **完成** | 2026-06-08 | 2026-06-08 |
| 2 | mini-celestial-mechanics | ⬜ 待开始 | — | — |
| 3 | mini-continuum | ⬜ 待开始 | — | — |
| 4 | mini-chaos | ⬜ 待开始 | — | — |
| 5 | mini-variational-principle | ⬜ 待开始 | — | — |
| 6 | mini-lagrangian | ⬜ 待开始 | — | — |
| 7 | mini-hamiltonian | ⬜ 待开始 | — | — |
| 8 | README 标准化 | ⬜ 待开始 | — | — |
| 9 | docs 补齐 | ⬜ 待开始 | — | — |
| 10 | examples 补齐 | ⬜ 待开始 | — | — |
| 11 | tests 补齐 | ⬜ 待开始 | — | — |
| 12 | demos 补齐 | ⬜ 待开始 | — | — |

### Step 1 完成清单

| 维度 | 前 | 后 | 状态 |
|------|-----|-----|------|
| src 文件数 | 1 个巨石文件 | 7 个模块化文件 | ✅ |
| src 总行数 | 49 | **1246** | ✅ (目标 800, 超额) |
| 类型系统 | InertiaTensor | +EulerAngles, RigidBodyState, PrincipalAxes | ✅ |
| 惯性张量 | 3 个函数 | +10 个标准形状 + 惯性椭球 | ✅ |
| Euler 方程 | 1 个函数 | +运动常数 + 3 种积分器 + 监测 | ✅ |
| 运动学 | 2 个函数 | +旋转矩阵↔欧拉角↔角速度 + 四元数 + Rodrigues | ✅ |
| 能量/动量 | 2 个函数 | +稳定性分析 + 网球拍定理 + 动能椭球 | ✅ |
| 陀螺 | 无 | 全新: Lagrange top + sleeping + 章动 + 陀螺仪 | ✅ |
| README | 3 行 | 完整 API 参考表 (60+ 函数) + 模块对照 + 目录树 | ✅ |
| examples | 1 个 | 4 个: free_precession + tennis_racket + spinning_top + gyroscope | ✅ |
| tests | 7 项 | ~85 项全覆盖 | ✅ |
| docs | 0 | 1 个: rigid-body-cheatsheet.md | ✅ |

