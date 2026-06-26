# mini-phase-explorer — 相空间探索器

交互式相空间分析工具：从 Hamilton 量直接生成相图、Poincaré 截面、能量等高线。

## Overview

给定任意 H(q, p)，本 demo 自动生成相空间全景：
- 能量等高线（1D 系统）
- 矢量场 (qdot, pdot) = (∂H/∂p, -∂H/∂q)
- 轨道积分（多初始条件采样）
- Poincaré 截面（2D 系统）
- 辛积分器对比

## Demo Program

### 输入

```julia
# 定义你的 Hamilton 量
H(q,p) = 0.5*(p[1]^2 + q[1]^2) + 0.25*q[1]^4  # Duffing 振子

# 自动分析
explore(H, q_range=(-2,2), p_range=(-2,2), E_levels=10)
```

### 输出

```
=== Phase Space Explorer ===
System: H = 0.5(p² + q²) + 0.25q⁴

Equilibrium points:
  (q=0, p=0): stable center (elliptic)
  
Energy contours (10 levels):
  E=0.1  → small-amplitude libration
  E=1.0  → nonlinear regime visible
  E=5.0  → strongly anharmonic

Orbit periods (Verlet, dt=0.01):
  E=0.1:  T=6.285 (linear T₀=6.283)
  E=1.0:  T=5.821 (amplitude-dependent)
  E=5.0:  T=4.312

Symplectic vs non-symplectic:
  Verlet:        ΔH_max = 3.2e-6  (bounded)
  Euler (expl):  ΔH_max = 0.47    (secular drift)
```

## Key Concepts

- **等能量面**：1D 系统 H(q,p)=E 是相空间中的闭合曲线 → 周期运动
- **非线性频率偏移**：Duffing 振子的周期随振幅增大而减小
- **辛积分器的优越性**：在长时间相空间结构保持上，Verlet 完胜任何非辛方法

## References

- Goldstein Ch.8 — Phase Space
- Lichtenberg & Lieberman — *Regular and Chaotic Dynamics* (1992)
