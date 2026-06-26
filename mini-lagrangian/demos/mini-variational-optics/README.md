# mini-variational-optics — 变分光学：Fermat原理
# 从最小作用量原理推导几何光学定律

## Overview

Fermat 原理是光学中的"最小作用量原理"：光在两点之间传播的路径是使光程（optical path length）取极小值的路径。这是变分原理在物理学中最优雅的应用之一——比牛顿力学的变分原理早了约 100 年。

## Physical Model

### Fermat's Principle

光程: S = ∫ n(r) ds，其中 n(r) 是折射率，ds 是路径微元。

δS = 0 → 光的实际路径。

### 从 Fermat 原理推导 Snell 定律

界面两侧折射率 n1, n2。光从 A 到 B，在界面上折射。

设界面在 y=0，A=(0,a)，B=(L,-b)，折射点 P=(x,0)。

光程: S = n1*√(x²+a²) + n2*√((L-x)²+b²)

最小化条件 dS/dx = 0:

```
n1*x/√(x²+a²) = n2*(L-x)/√((L-x)²+b²)
→ n1*sin(θ1) = n2*sin(θ2)   ← Snell定律!
```

### 变分推导

更一般地，从 δS = 0 使用 Euler-Lagrange 方程:

光线的 Lagrangian: L(x, y, y') = n(x,y)*√(1 + y'²)

EL 方程给出光线方程:

```
d/dx [n*y'/√(1+y'²)] = √(1+y'²) * ∂n/∂y
```

### 梯度折射率介质中的光线弯曲

海市蜃楼：近地面空气温度高 → 折射率小 → 光线向上弯曲。

折射率梯度: n(y) = n0 + α*y（α<0 对应热地面）

光线方程数值解 → 看到"水面"倒影。

## Verification

### 数值验证：Snell 定律

入射角 30°，n1=1.0, n2=1.5:

```
Fermat 路径数值搜索:  x_optimal = 0.387L
Snell 预测:  sin(θ2) = n1/n2*sin(30°) = 0.333 → θ2 = 19.47°
几何:  tan(θ2) = (L-x)/b → x = 0.387L  ✓ 一致
```

### 光程最小 vs 极大 vs 鞍点

球面镜反射：光程可以是极小、极大或鞍点，取决于曲率半径。
Fermat 原理是更准确的"稳定光程原理"（stationary，而非 strictly minimum）。

## Demo Program

```
=== Variational Optics ===

1. Fermat → Snell:
   n1=1.0, n2=1.5, θ1=30°
   Optimal x/L = 0.387 (theory: 0.387) ✓

2. Mirage (graded index):
   n(y) = 1.0 - 0.0001*y
   Ray traced: bends upward
   Apparent reflection at y = -50m

3. Spherical mirror:
   R = 1.0, object at 2R
   Path length analysis: local minimum

4. Gravitational lensing analogy:
   n_eff(r) = 1 + 2GM/(rc²)  (weak field)
   → light bending by mass
```

## Source Files

| File | Purpose |
|------|---------|
| `demos/mini-variational-optics/fermat.jl` | Fermat 主程序 |
| `demos/mini-variational-optics/ray_trace.jl` | 梯度折射率光线追迹 |
| `src/action.jl` | 作用量变分计算框架 |

## Key Concepts

- **变分原理的统一性**：力学的最小作用量 (∫L dt) 和光学的最小光程 (∫n ds) 是同一数学结构的两个实例。这暗示了量子力学的路径积分表述和几何光学-经典力学的类比。
- **Hamilton 的光学-力学类比**：正是这个类比启发 Schrödinger 建立了波动力学——经典力学是几何光学的极限，波动力学是波动光学的对应物。
- **Fermat 原理不是"最小"而是"稳定"**：与力学中的 Hamilton 原理一样，δS=0 可以是极小、极大或鞍点。

## References

- Born & Wolf — *Principles of Optics*, Ch.3 — Fermat's Principle
- Goldstein — *Classical Mechanics*, Ch.2 (变分原理) + Ch.9 (Hamilton-Jacobi 与光学类比)
- Feynman — *Lectures on Physics*, Vol.1 Ch.26 — Optics: The Principle of Least Time
- MIT 8.03 — Vibrations and Waves, Lecture 17 — Fermat's Principle
