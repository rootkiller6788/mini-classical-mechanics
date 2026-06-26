# mini-chaotic-pendulum — 双摆混沌与非线性动力学

从简单单摆到双摆混沌，完整演示非线性动力学的核心现象：相空间、庞加莱截面、Lyapunov 指数、分岔。

## Overview

单摆是大一物理的经典线性近似案例。但真实物理世界是非线性的——大角度单摆周期偏离 2π√(L/g)，阻尼+驱动产生混沌，双摆则是经典的确定性混沌系统。本 demo 使用 `mini-newtonian` 的 ODE 积分器，逐步展示从线性到混沌的完整路径。

## Physical Models

### 1. 大角度单摆 — 超越小角近似

运动方程：θ'' = -(g/L) sin(θ)

小角度近似 sin(θ)≈θ 给出周期 T₀=2π√(L/g)。精确周期（椭圆积分）：

```
T(θ₀) = T₀ * (2/π) * K(sin²(θ₀/2))
      ≈ T₀ * (1 + θ₀²/16 + 11θ₀⁴/3072 + ...)
```

其中 K 是第一类完全椭圆积分。本 demo 数值验证此级数展开。

### 2. 阻尼+驱动单摆 — 通往混沌之路

运动方程：θ'' + b*θ' + (g/L) sin(θ) = F₀ cos(ω_d t)

参数空间 (b, F₀, ω_d) 中，特定的参数组合产生混沌运动：
- 周期倍增分岔级联 → 混沌
- 奇怪吸引子 → 相空间中的分形结构

### 3. 双摆 — 确定性混沌的典范

运动方程（拉格朗日力学导出）：

```
(m₁+m₂)L₁θ₁'' + m₂L₂θ₂''cos(θ₁-θ₂) + m₂L₂(θ₂')²sin(θ₁-θ₂) + (m₁+m₂)g sin(θ₁) = 0
m₂L₂θ₂'' + m₂L₁θ₁''cos(θ₁-θ₂) - m₂L₁(θ₁')²sin(θ₁-θ₂) + m₂g sin(θ₂) = 0
```

**关键性质**：
- 4 维相空间 (θ₁, θ₂, ω₁, ω₂)
- 无耗散 → 相空间体积守恒（Liouville 定理）
- 两个初始条件相差 10⁻⁶ → 几秒后轨迹完全不相关
- 对初始条件的敏感依赖性 → 蝴蝶效应

### 4. Lyapunov 指数

量化混沌的指标：相邻轨迹的指数分离率。

```
λ = lim_{t→∞} lim_{δ₀→0} (1/t) ln(|δ(t)|/|δ₀|)
```

λ > 0 → 混沌。双摆的最大 Lyapunov 指数 ≈ 0.2-0.5 s⁻¹（取决于参数）。

## Verification Methods

### 相图可视化

- 单摆：相空间 (θ, ω) → 椭圆（小角度）→ 扭曲的闭合曲线（大角度）
- 阻尼单摆：相空间 → 螺旋收敛到不动点 (0,0)
- 驱动单摆：相空间 → 极限环 → 周期倍增 → 奇怪吸引子
- 双摆：4D 相空间投影到 (θ₁, ω₁) 或 (θ₁, θ₂)

### Poincaré 截面

对驱动单摆，以驱动周期 T_d = 2π/ω_d 采样相空间：
- 周期运动 → 截面上有限个点
- 准周期运动 → 截面上闭合曲线
- 混沌运动 → 截面上弥散的点云

### 能量守恒（无阻尼单摆/双摆）

对无阻尼系统，总能量 E = T + U 应为常数。数值漂移反映积分器精度：
- Euler：指数发散
- RK4：线性漂移 ~10⁻¹²/步
- **Velocity Verlet**：有界振荡，无长期漂移 → **本 demo 使用 Verlet**

## Demo Program

### 参数空间探索

```julia
# 从规则到混沌的参数扫描
params = [
    (F0=0.5,  b=0.1,  label="Limit Cycle"),
    (F0=0.9,  b=0.1,  label="Period Doubling"),
    (F0=1.07, b=0.1,  label="Period-4"),
    (F0=1.15, b=0.1,  label="Chaos"),
    (F0=1.50, b=0.1,  label="Strong Chaos"),
]
```

### 双摆初始条件敏感性

```julia
θ1₀, θ2₀ = 90°, 90°          # 基准初始条件
θ1₀', θ2₀' = 90.0001°, 90°   # 偏差 10⁻⁴ 度

# 两个轨迹在 t≈5s 时开始显著分离
# 在 t≈15s 时完全不相关
# → Lyapunov 时间 ≈ 5s
```

### Expected Output

```
=== Chaotic Pendulum Suite ===

--- 1. Simple Pendulum Period vs Amplitude ---
θ₀ (deg) | T (sim) | T (analytic) | Error
---------|---------|-------------|-------
   5     |  2.007  |   2.007     | 0.01%
  30     |  2.044  |   2.044     | 0.02%
  60     |  2.155  |   2.153     | 0.09%
  90     |  2.373  |   2.369     | 0.17%
 120     |  2.785  |   2.778     | 0.25%

--- 2. Damped Driven Pendulum ---
b=0.1, F₀=0.5:  Limit cycle (period=1)
b=0.1, F₀=1.07: Period-4 orbit
b=0.1, F₀=1.15: Chaotic attractor (Lyapunov λ≈0.12)

--- 3. Double Pendulum ---
Initial: θ₁=90°, θ₂=90°
Integration: Velocity Verlet, dt=0.0005s, t=20s
Energy drift: 3.2e-8 (relative max)
Lyapunov exponent: 0.31 s⁻¹
Initial condition sensitivity: diverges at t≈4.8s

--- 4. Poincaré Section (Driven Pendulum) ---
Stroboscopic sampling at ω_d period:
  F₀=0.5:  1 point  (period-1 limit cycle)
  F₀=1.07: 4 points (period-4)
  F₀=1.15: scattered cloud (chaos)
```

## Source Files

| File | Purpose |
|------|---------|
| `demos/mini-chaotic-pendulum/chaotic_pendulum.jl` | 主程序：参数扫描 + 相图生成 |
| `demos/mini-chaotic-pendulum/poincare.jl` | Poincaré 截面分析 |
| `demos/mini-chaotic-pendulum/lyapunov.jl` | Lyapunov 指数计算 |
| `src/integrators.jl` | Velocity Verlet（辛积分器） |
| `src/energy.jl` | 能量守恒验证 |

## Key Concepts

- **线性 vs 非线性**：sin(θ)≈θ 是近似，真实摆在大角度偏离。小角度近似是绝大多数本科物理习题的隐藏假设。
- **确定性 ≠ 可预测**：双摆方程完全确定，但初始条件的不确定性指数放大。这是 Lorenz 发现"蝴蝶效应"的数学根源。
- **Verlet 对混沌系统的优势**：混沌系统指数放大数值误差。辛积分器保持相空间结构不变形，即使具体轨迹不可预测，统计性质（能量分布、Lyapunov 指数）仍可靠。
- **混沌的工程意义**：太阳系长期稳定性、湍流、心律不齐、激光动力学——都涉及混沌。

## References

- Strogatz, S.H. — *Nonlinear Dynamics and Chaos* (2015), Ch.2-3, 8-9
- Baker & Gollub — *Chaotic Dynamics: An Introduction* (1996), Ch.3
- Shinbrot et al. — "Using the Sensitive Dependence of Chaos to Direct Trajectories" (1992)
- MIT 8.012 — Lecture 25: Chaos
- MIT 18.385J / 2.036J — Nonlinear Dynamics and Chaos
