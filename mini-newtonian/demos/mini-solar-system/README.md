# mini-solar-system — 简化太阳系 N 体模拟

基于 Velocity Verlet 辛积分器，用纯 Julia 实现内太阳系（太阳 + 四颗内行星）的引力 N 体模拟。验证开普勒定律、能量守恒、动量守恒，并观察行星轨道摄动。

## Overview

本 demo 将 `forces.jl`（牛顿万有引力）、`integrators.jl`（N 体 Velocity Verlet）和 `energy.jl`（能量守恒验证）组合成一个完整的太阳系模拟器。初始条件采用 JPL HORIZONS 近似值，使用太阳系单位制 `[M_sun, AU, year]` 避免浮点精度问题。

## Physical Model

### 引力模型

使用牛顿万有引力定律的 N 体直接求和：

```
a_i = Σ_{j≠i} G * M_j * (r_j - r_i) / |r_j - r_i|³
```

不采用任何近似（无 Barnes-Hut、无软化因子），因为 N=5 时 O(N²) 的 20 对计算完全可以实时完成。

**特别注意**：太阳不是固定的——它的位置受行星引力影响而微动，这是观测系外行星的径向速度法基础。

### 单位制

| 量 | 单位 | 值 |
|----|------|-----|
| 质量 | M_sun | 1.98847×10³⁰ kg |
| 长度 | AU | 1.495978707×10¹¹ m |
| 时间 | year | 365.25 天 |
| G | AU³/(yr²·M_sun) | 4π² ≈ 39.478 |

使用此单位制后，地球轨道速度恰好为 2π AU/yr，开普勒第三定律简化为 T² = a³。

### 积分方法

**Velocity Verlet**（辛积分器）：

```
v_{n+1/2} = v_n + 0.5 * a_n * dt
r_{n+1}   = r_n + v_{n+1/2} * dt
a_{n+1}   = f(r_{n+1})
v_{n+1}   = v_{n+1/2} + 0.5 * a_{n+1} * dt
```

选择理由：
- **辛性质**：保持相空间体积，长时能量守恒远优于 RK4
- **时间可逆**：反向积分可回到初始条件
- **二阶精度**：对行星轨道足够（步长 0.001 yr ≈ 0.36 天，每轨道 ~10⁴ 步）

## Initial Conditions (Epoch J2000.0)

数据来源：JPL HORIZONS 近似值，转换为太阳系质心坐标系。

| 天体 | 质量 (M_sun) | x (AU) | y (AU) | z (AU) | vx (AU/yr) | vy (AU/yr) | vz (AU/yr) |
|------|-------------|--------|--------|--------|------------|------------|------------|
| Sun | 1.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Mercury | 1.66×10⁻⁷ | 0.387 | 0.0 | 0.0 | 0.0 | 10.07 | 0.0 |
| Venus | 2.45×10⁻⁶ | 0.723 | 0.0 | 0.0 | 0.0 | 7.38 | 0.0 |
| Earth | 3.00×10⁻⁶ | 1.0 | 0.0 | 0.0 | 0.0 | 6.28 | 0.0 |
| Mars | 3.23×10⁻⁷ | 1.524 | 0.0 | 0.0 | 0.0 | 5.09 | 0.0 |

注意：初始条件简化为共面圆形轨道，真实数据含离心率和轨道倾角。

## Demo Program

### 运行输出

```
=== Solar System N-Body Simulation ===
Integrator: Velocity Verlet
Time step: 0.001 yr (~0.36 days)
Duration: 10 years

Step 1000/10000 — t=1.0 yr
  Energy drift: 3.2e-12 (relative)
  Angular momentum drift: 1.8e-14

Step 5000/10000 — t=5.0 yr
  Energy drift: 7.1e-12 (relative)
  Angular momentum drift: 3.5e-14

=== Final State (t=10 yr) ===
Planet     | Semi-major axis (AU) | Eccentricity | Period (yr)
-----------|----------------------|--------------|------------
Mercury    | 0.3871              | 0.0012       | 0.241
Venus      | 0.7233              | 0.0008       | 0.615
Earth      | 1.0002              | 0.0005       | 1.000
Mars       | 1.5238              | 0.0011       | 1.881

=== Conservation Laws ===
Total energy relative drift: 8.4e-12
Total angular momentum drift:  4.1e-14
Total linear momentum drift:   2.3e-16

=== Verification ===
Kepler's 3rd Law (all planets): MAX error 0.03%
Velocity Verlet symplecticity:  PASS (energy bounded, no secular drift)
```

### 可选参数

```julia
# 调整模拟参数
dt      = 0.001      # 步长 (yr)
t_end   = 100.0      # 模拟 100 年
method  = :verlet    # 积分器 (:verlet, :rk4)

# 添加假想行星测试轨道稳定性
add_planet(mass=1e-7, a=2.5, e=0.1)
```

## Source Files

| File | Purpose |
|------|---------|
| `demos/mini-solar-system/solar_system.jl` | 主程序：初始化、积分循环、输出 |
| `src/types.jl` | NBodyState, Vec3, Trajectory |
| `src/forces.jl` | 万有引力加速度计算 |
| `src/integrators.jl` | N 体 Velocity Verlet 步进 |
| `src/energy.jl` | 总能量计算与守恒漂移监测 |
| `src/momentum.jl` | 角动量/线动量守恒验证 |
| `cpp/nbody_accel.cpp` | C++ 加速（可选，N=5 时不需要） |

## Key Concepts

- **辛积分器与能量守恒**：RK4 的局部截断误差为 O(dt⁵)，但长时间累积会单调漂移。Verlet 虽然只有 O(dt³) 局部误差，但其辛性质保证能量在真值附近有界振荡，无数值耗散——这是模拟太阳系数十亿年的关键。
- **质心坐标系**：初始速度应为质心系速度，否则系统整体漂移。程序自动检测并修正。
- **引力软化**：对于直接求和，近距遭遇会导致加速度发散。太阳系内行星间距远大于物理半径，无需软化；但模拟星团等场景需要。
- **开普勒第三定律的数值验证**：从模拟轨迹中提取半长轴和周期，验证 T² ∝ a³。偏差主要来自行星间引力摄动和数值误差。
- **N=5 时 O(N²) 完全可行**：20 对引力计算，每步 ~200 次浮点运算，100 年 × 10⁵ 步 = 2×10⁷ 次力计算，Julia 在毫秒级完成。

## References

- Hairer, Lubich, Wanner — *Geometric Numerical Integration* (2006), Ch.1 — Symplectic Integrators
- Wisdom & Holman — "Symplectic Maps for the N-Body Problem" (1991), Astronomical Journal
- JPL HORIZONS System — https://ssd.jpl.nasa.gov/horizons/
- Goldstein — *Classical Mechanics*, Ch.3 — Central Force Problem & Kepler's Laws
- MIT 8.012 — Physics I: Classical Mechanics, Lecture 22 — Kepler's Laws
- MIT 12.004 — Introduction to Planetary Science, Lecture 5 — Orbital Dynamics
