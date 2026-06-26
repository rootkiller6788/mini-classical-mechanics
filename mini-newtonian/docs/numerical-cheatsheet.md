# Numerical Methods & Physics Cheatsheet

快速参考：所有 ODE 积分器对比、力定律汇总、守恒量公式。对标 complexity-cheatsheet.md。

---

## Integrator Comparison

| Integrator | Order | Stages/Step | Stability | Symplectic | Energy Drift | Best For |
|------------|-------|-------------|-----------|------------|-------------|----------|
| Euler (Forward) | 1 | 1 | Conditional (hλ∈[-2,0]) | No | Monotonic | 教学演示，不可用于生产 |
| Euler-Cromer | 1 | 1 | Conditional | No | ~Bounded | 入门级物理模拟 |
| Midpoint (RK2) | 2 | 2 | Conditional | No | Monotonic | 简单改进 |
| RK4 | 4 | 4 | Conditional (hλ∈[-2.785,0]) | No | Monotonic | 通用中短时间积分 |
| **Velocity Verlet** | **2** | **2** | **hω<2** | **Yes** | **~Bounded** | **N 体、分子动力学、守恒系统** |
| Leapfrog | 2 | 1 | hω<2 | Yes | ~Bounded | 同 Verlet，少一次力评估 |
| RK45 (Dormand-Prince) | 5(4) | 7 | Conditional | No | Adaptive | 非光滑系统、碰撞、变步长需求 |

**选择流程图：**

```
系统是保守的(Hamiltonian)? ──Yes── 长时积分(>10³周期)? ──Yes── Velocity Verlet
    │                               │
    No                              No
    │                               │
含速度依赖力? ──Yes── RK4           含接触/碰撞? ──Yes── RK45 (自适应)
    │
    No
    │
Euler-Cromer (原型阶段)
```

---

## Force Laws Reference

| Force | Formula | Potential Energy | Parameters |
|-------|---------|-----------------|------------|
| Newton Gravity | F = -GMm/r² * r̂ | U = -GMm/r | G, M, m |
| Uniform Gravity (near surface) | F = mg (down) | U = mgh | g = 9.81 m/s² |
| Hooke's Law (Spring) | F = -k(r - r_eq) | U = ½k|r-r_eq|² | k, r_eq |
| Linear Drag | F = -bv | N/A (non-conservative) | b |
| Quadratic Drag | F = -c|v|v | N/A | c |
| Kinetic Friction | F = -μ_k N v̂ | N/A | μ_k, N |
| Lorentz Force | F = q(E + v×B) | U = qφ (E-field only) | q, E, B |
| Damped Harmonic | F = -kx - bv | U = ½kx² | k, b |
| Centrifugal (rotating frame) | F = -m ω×(ω×r) | U = -½m|ω×r|² | ω |
| Coriolis (rotating frame) | F = -2m ω×v | N/A | ω |

---

## Conservation Laws

| Quantity | Formula | Conserved When |
|----------|---------|----------------|
| Linear Momentum | p = Σ mᵢvᵢ | Σ F_ext = 0 |
| Angular Momentum | L = Σ rᵢ × mᵢvᵢ | Σ τ_ext = 0 (no external torque) |
| Total Mechanical Energy | E = T + U | Only conservative forces, no drag |
| Center of Mass Velocity | v_cm = p/M_total | No external forces |
| Phase Space Volume (Liouville) | ∫ dp∧dq | Hamiltonian system → symplectic integrator preserves |

### Energy Conservation Quality by Integrator

对谐振子 x'' = -ω²x, T=1000:

```
Integrator        | Max ΔE/E | Trend
------------------|----------|--------
Euler             | 10⁰      | Exponential growth
Euler-Cromer      | 10⁻²     | Bounded oscillation
RK4               | 10⁻⁶     | Linear drift (~10⁻¹² / step)
Velocity Verlet   | 10⁻⁶     | Bounded oscillation (no drift)
```

**关键洞察**：RK4 的 per-step 误差远小于 Verlet，但经过 10⁵ 步后能量单调漂移。Verlet 能量在真值附近有界振荡——这是 **conservation vs accuracy** 的区别。

---

## Coordinate Systems

| System | Position | Velocity | Acceleration |
|--------|----------|----------|-------------|
| Cartesian | (x, y, z) | (vx, vy, vz) | (ax, ay, az) |
| Cylindrical | (ρ, φ, z) | (ρ̇, ρφ̇, ż) | (ρ̈-ρφ̇², ρφ̈+2ρ̇φ̇, z̈) |
| Spherical | (r, θ, φ) | (ṙ, rθ̇, r sinθ φ̇) | (r̈-rθ̇²-r sin²θ φ̇², ...) |

---

## Common Physical Constants

| Symbol | Value | Description |
|--------|-------|-------------|
| G | 6.67430×10⁻¹¹ m³/(kg·s²) | Gravitational constant |
| g | 9.80665 m/s² | Earth surface gravity |
| AU | 1.495978707×10¹¹ m | Astronomical unit |
| M_sun | 1.98847×10³⁰ kg | Solar mass |
| M_earth | 5.9722×10²⁴ kg | Earth mass |
| R_earth | 6.371×10⁶ m | Earth radius |
| year | 3.15576×10⁷ s | Sidereal year |
| day | 86400 s | Solar day |

### Dimensionless Solar System Units

```
[Mass]  = M_sun     → G = 4π² ≈ 39.478
[Length] = AU       → Earth orbit radius = 1
[Time]   = year     → Earth orbital period = 1
[Velocity] = AU/yr  → Earth orbital speed = 2π
```

此单位制下开普勒第三定律简化为精确的 T² = a³（无比例常数）。

---

## Error Scaling Laws

```
Global Error ∝ h^p,  where p = order of accuracy

Example: reduce h by 2x →
  Euler (p=1):   error / 2
  RK2 (p=2):     error / 4
  RK4 (p=4):     error / 16
  Verlet (p=2):  error / 4  (but energy bounded!)
```

**精度 vs 代价权衡**（对谐振子 1 周期，目标误差 10⁻⁶）：

| Method | Steps needed | Force evals | CPU time (relative) |
|--------|-------------|-------------|---------------------|
| Euler | ~10⁶ | 10⁶ | 100× |
| RK2 | ~10³ | 2×10³ | 0.3× |
| RK4 | ~10² | 4×10² | 0.1× |
| Verlet | ~10² | 2×10² | **0.05×** |
| RK45 | ~50 | ~350 | 0.08× |

Verlet 不仅精度足够，且每步只需 2 次力评估（vs RK4 的 4 次）+ 辛性质。
**在保守系统的长时积分中，Verlet 是帕累托最优选择。**

---

## References

- *Geometric Numerical Integration*, Hairer, Lubich, Wanner (2006)
- *Numerical Recipes*, Press et al., Ch.16 — Integration of ODEs
- *Solving Ordinary Differential Equations I*, Hairer, Norsett, Wanner (1993)
- MIT 18.330 — Numerical Analysis, Lectures 14-20
