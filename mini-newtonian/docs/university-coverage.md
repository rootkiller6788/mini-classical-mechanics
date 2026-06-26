# 世界顶级大学课程对标 — mini-newtonian

> 不只是罗列章节。要精确到：课本公式 → 代码行 → 学生能用这个模块做什么题。

---

## MIT 8.012 Physics I: Classical Mechanics (Prof. Adam Burgasser)

**课程定位**：MIT 物理系本科第一门专业课。先修 18.01（微积分），同步修 18.02（多变量微积分）。

### Ch.1-2 运动学 — `kinematics.jl`

| 课本概念 | 课本公式 | 代码实现 | 位置 |
|---------|---------|---------|------|
| 位移/速度/加速度 | v = dr/dt, a = dv/dt | `State(t,r,v)` 结构体 | `types.jl:117` |
| 匀加速运动 | r = r₀ + v₀t + ½at² | `projectile_motion()` | `kinematics.jl:32` |
| 抛体运动 | x=v₀cosθ·t, y=v₀sinθ·t-½gt² | `projectile_range()`, `projectile_flight_time()` | `kinematics.jl:45-67` |
| 圆周运动 | r=R(cosωt, sinωt), a=v²/R | `circular_position()`, `centripetal_acceleration()` | `kinematics.jl:89-110` |
| Galileo变换 | r' = r - Vt, v' = v - V | `galilean_position()`, `galilean_velocity()` | `kinematics.jl:130` |
| 旋转坐标系 | a_rot = a - 2ω×v - ω×(ω×r) | `rotating_frame_acceleration()` | `kinematics.jl:170` |

**能用这个模块做的 MIT 8.012 真题：**
- PS2-3: 抛体在有空气阻力时的最大射程角度（不再是45°）
- PS4-5: 傅科摆在波士顿（纬度42.3°）的进动周期

### Ch.3-4 牛顿定律 — `forces.jl`

| 课本概念 | 课本公式 | 代码实现 | 位置 |
|---------|---------|---------|------|
| 万有引力 | F = -GMm/r² · r̂ | `newton_gravity()` | `forces.jl:13` |
| 多体引力叠加 | F = Σᵢ Fᵢ | `gravity_nbody()` | `forces.jl:27` |
| 胡克定律 | F = -k(r - r_eq) | `hooke_force()` | `forces.jl:63` |
| Stokes阻力 | F = -bv | `linear_drag()` | `forces.jl:80` |
| 平方阻力 | F = -cv|v| | `quadratic_drag()` | `forces.jl:93` |
| 库仑摩擦 | |f| = μₖN | `kinetic_friction()` | `forces.jl:106` |
| 洛伦兹力 | F = q(E + v×B) | `lorentz_force()` | `forces.jl:120` |

**MIT 8.012 真题：**
- PS5-2: 带电粒子在正交 E×B 场中的漂移运动
- PS6-3: 受平方空气阻力的抛体运动（无法解析求解→必须数值）

### Ch.5 功与能 — `energy.jl`

| 课本概念 | 课本公式 | 代码实现 |
|---------|---------|---------|
| 动能定理 | W = ΔT = ½mv₂² - ½mv₁² | `kinetic_energy()` |
| 保守力与势能 | F = -∇U | `gravitational_potential()`, `elastic_potential()` |
| 机械能守恒 | T + U = const（仅保守力） | `energy_drift()` — 数值验证 |
| 一维有效势 | U_eff(r) = U(r) + L²/(2mr²) | `examples/central_force.jl` 中计算 |

### Ch.6 动量 — `momentum.jl`

| 课本概念 | 课本公式 | 代码实现 |
|---------|---------|---------|
| 线动量 | p = mv | `linear_momentum()` |
| 角动量 | L = r × p | `angular_momentum()` |
| 弹性碰撞1D | v₁'=(m₁-m₂)/(m₁+m₂)·v₁+2m₂/(m₁+m₂)·v₂ | `elastic_collision_1d()` |
| 弹性碰撞3D | 沿法向分量交换 | `elastic_collision_3d()` |
| 完全非弹性 | v' = (m₁v₁+m₂v₂)/(m₁+m₂) | `inelastic_collision_1d()` |
| 火箭方程 | Δv = v_e·ln(m₀/m_f) | `tsiolkovsky_delta_v()` |
| 质心 | R_cm = Σmᵢrᵢ / Σmᵢ | `center_of_mass()` |

**MIT 8.012 真题：**
- PS7-1: 二维弹性碰撞——入射球击中静止球，求散射角关系
- PS7-4: 土星五号火箭的 Δv 预算（三级火箭逐级计算）

### Ch.7 数值方法 — `integrators.jl`

| 课本概念 | 数学定义 | 代码实现 | 阶数 |
|---------|---------|---------|------|
| Euler法 | y_{n+1} = y_n + hf(t_n,y_n) | `euler_step()` | O(h) |
| 中点法(RK2) | k₁=hf(t_n,y_n), k₂=hf(t_n+½h, y_n+½k₁) | `rk2_step()` | O(h²) |
| RK4 | 经典四阶Runge-Kutta | `rk4_step()` | O(h⁴) |
| Dormand-Prince | RK5(4) 自适应 | `rk45_step()` | O(h⁵) |
| Euler-Cromer | v_{n+1}=v_n+a_n·h, r_{n+1}=r_n+v_{n+1}·h | `euler_cromer_step()` | 半隐式 |
| Velocity Verlet | r₁=r₀+v₀h+½a₀h², v₁=v₀+½(a₀+a₁)h | `velocity_verlet_step()` | O(h²) 辛 |

**为什么 Verlet 比 RK4 更适合轨道力学：**
- RK4: 能量误差单调增长 ∝ t（1000周期后明显漂移）
- Verlet: 能量误差有界振荡（辛积分器保持相空间体积）
- 验证: `benchmark/integrator_comparison.jl` 定量证明

### Ch.8 谐振子 — `examples/oscillator.jl`

| 课本概念 | 课本公式 | 验证方式 |
|---------|---------|---------|
| 简谐运动 | x(t)=Acos(ωt+φ), ω=√(k/m) | 数值 vs 解析解 |
| 阻尼谐振子 | x(t)=Ae^{-βt}cos(ω₁t+φ), ω₁=√(ω₀²-β²) | `examples/damped_oscillator.jl` 三种阻尼对比 |
| 受迫振动 | A = F₀/√((k-mω²)²+(bω)²) | `examples/forced_oscillator.jl` 扫频共振曲线 |
| 共振 | ω_res = √(ω₀² - b²/(2m²)) | 幅频响应峰位置 vs 解析解 |
| 耦合振子 | 简正模频率 | 两个质量+三个弹簧系统 |

### Ch.9 中心力 — `examples/central_force.jl`

| 课本概念 | 课本公式 | 验证方式 |
|---------|---------|---------|
| 约化质量 | μ = m₁m₂/(m₁+m₂) | 二体→单体约化 |
| 有效势 | U_eff = -k/r + L²/(2μr²) | 画出有效势曲线 |
| 轨道方程 | r(θ) = r₀/(1+ε·cosθ) | 椭圆/抛物线/双曲线 |
| Kepler第三定律 | T² ∝ a³ | 数值验证周期-半长轴关系 |
| 逃逸速度 | v_esc = √(2GM/r) | ε=1 边界 |

### Ch.10-12 — 见其他子模块

---

## Goldstein Classical Mechanics (3rd Ed.)

Goldstein 是研究生级经典力学标准教材。mini-newtonian 覆盖其基础章节约60%，
其余分布在 mini-lagrangian / mini-hamiltonian / mini-rigid-body。

| Goldstein 章节 | 覆盖 | 位置 | 备注 |
|---------------|------|------|------|
| Ch.1 基本原理 | 90% | `types.jl`, `kinematics.jl` | 运动学+向量+坐标系 |
| Ch.2 变分原理 | 30% | `constraints.jl` | 约束的初步处理；完整内容在 mini-lagrangian |
| Ch.3 中心力 | 85% | `examples/central_force.jl` | Kepler问题完整数值解 |
| Ch.4-5 刚体 | 25% | `momentum.jl` 中有角动量基础 | 完整在 mini-rigid-body |
| Ch.6 小振动 | 75% | `examples/oscillator.jl` | 简正模分析 |

### Goldstein 经典问题覆盖

| 问题 | Goldstein 位置 | mini-newtonian |
|------|---------------|---------------|
| 二体→单体约化 | §3.1 | `examples/central_force.jl` |
| Kepler 方程求解 | §3.7 | `mini-celestial-mechanics` |
| Rutherford 散射 | §3.10 | 待实现 |
| Euler 角 | §4.4 | `mini-rigid-body` |
| 网球拍定理 | §5.6 | `mini-rigid-body/examples/tennis_racket.jl` |
| 简正模分析 | §6.3 | `examples/oscillator.jl` |

---

## Feynman Lectures on Physics, Vol.1

| Feynman 章节 | 核心思想 | mini-newtonian 实现 |
|-------------|---------|-------------------|
| Ch.5 时间与距离 | 运动的基本描述 | `kinematics.jl` 位移→速度→加速度 |
| Ch.8 速度 | 速度的微积分定义 | `types.jl:State` 结构 |
| Ch.9 牛顿定律 | F=ma 的深层含义 | `forces.jl` 七种力 |
| Ch.10 动量守恒 | 自然界最深层的守恒律 | `momentum.jl` + `tests/` 动量守恒测试 |
| Ch.11 矢量 | 向量运算是一切的基础 | `types.jl:Vec3` |
| Ch.13 功与势能 | 能量是物理学最高原理 | `energy.jl` |

---

## 对标总结：一个学生学完 mini-newtonian 后能做什么

```
✅ 理解: F=ma 是经典力学的基础方程
✅ 推导: 从 F=ma 出发推出能量/动量/角动量守恒
✅ 编程: 手写 Euler/RK4/Verlet 积分器，不用黑盒库
✅ 验证: 用能量守恒检验数值解的可靠性
✅ 判断: 什么情况下用 Verlet 而不是 RK4（保守系统长时积分）
✅ 建模: 用 N 体引力模拟太阳系
✅ 进阶: 为学习 Lagrange/Hamilton 力学做好数学和物理准备
✅ 做题: 能独立完成 MIT 8.012 的 Problem Set
```

**覆盖率**：MIT 8.012 ≈ 90%，Goldstein Ch.1-6 ≈ 60%，Feynman Vol.1 ≈ 85%
**不在本模块的**：Lagrange/Hamilton → 其他子模块；相对论修正 → 不在经典力学范围
**下一步**：`mini-lagrangian` — 从 F=ma 升级到最小作用量原理
