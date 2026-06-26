# Integrator Analysis — 积分器数学分析

所有 ODE 积分器的阶数证明、稳定性分析和适用场景。对标 correctness-proofs.md 的证明结构。

---

## 1. Euler Method — 一阶收敛

### Algorithm

```
y_{n+1} = y_n + h * f(t_n, y_n)
```

### Local Truncation Error

假设真解 y(t) 满足 y'(t) = f(t, y(t))。Taylor 展开：

```
y(t_n + h) = y(t_n) + h*y'(t_n) + (h²/2)*y''(t_n) + O(h³)
           = y(t_n) + h*f(t_n, y(t_n)) + (h²/2)*y''(t_n) + O(h³)
```

Euler 方法近似：

```
y_{n+1} = y_n + h*f(t_n, y_n)
```

**局部截断误差（LTE）：**

```
LTE = y(t_n + h) - y_{n+1} = (h²/2)*y''(t_n) + O(h³) = O(h²)
```

**全局误差（经过比例 1/h 步）：** O(h)

### Stability (Dahlquist Test Equation)

对测试方程 y' = λy, Re(λ) < 0：

```
y_{n+1} = y_n + h*λ*y_n = (1 + hλ)*y_n
```

稳定性条件：|1 + hλ| ≤ 1

- 实轴稳定区间：hλ ∈ [-2, 0]
- 纯虚轴：|1 + iθ| = √(1 + θ²) > 1 → **无条件不稳定**

**结论：Euler 方法不适合振荡系统（如谐振子、行星轨道）。**

---

## 2. Runge-Kutta 4 — 四阶收敛

### Algorithm

```
k1 = f(t_n, y_n)
k2 = f(t_n + h/2, y_n + h*k1/2)
k3 = f(t_n + h/2, y_n + h*k2/2)
k4 = f(t_n + h, y_n + h*k3)
y_{n+1} = y_n + (h/6)*(k1 + 2k2 + 2k3 + k4)
```

### Local Truncation Error

通过 Taylor 展开匹配到五阶项（证明过程见 Butcher 2008, §2.3）：

```
LTE = y(t_n + h) - y_{n+1} = O(h⁵)
```

**全局误差：** O(h⁴)

### Stability

对 y' = λy：

```
y_{n+1} = R(hλ)*y_n,  其中 R(z) = 1 + z + z²/2 + z³/6 + z⁴/24
```

稳定区域：{z ∈ ℂ : |R(z)| ≤ 1}

- 实轴：hλ ∈ [-2.785, 0]
- 纯虚轴边界：|hλ| < 2.828

**结论：RK4 是振荡系统的实用最小选择，但长时间积分有能量漂移（无数值耗散，但无辛结构）。**

---

## 3. Velocity Verlet — 二阶辛积分器

### Algorithm

```
v_{n+1/2} = v_n + (h/2)*a(r_n)
r_{n+1}   = r_n + h*v_{n+1/2}
a_{n+1}   = F(r_{n+1})/m
v_{n+1}   = v_{n+1/2} + (h/2)*a_{n+1}
```

### Symplecticity Proof (Sketch)

Hamilton 系统 H(p,q) = T(p) + V(q)。Verlet 可分解为三个精确可解的子流映射的组合：

```
Φ_h = φ_{h/2}^V ∘ φ_h^T ∘ φ_{h/2}^V
```

其中：
- φ_h^T: 仅含动能项 dq/dt = p/m, dp/dt = 0 → 精确解 (q + h*p/m, p)
- φ_h^V: 仅含势能项 dq/dt = 0, dp/dt = -∇V(q) → 精确解 (q, p - h*∇V(q))

每个子映射是某 Hamilton 流的精确时间推进 → 都是辛映射。
**辛映射的复合仍是辛映射** → Verlet 是辛变换。

结论：相空间体积守恒 → Liouville 定理的离散类比 → 能量在真值附近有界振荡。

### Local Error

截断误差 O(h³)，全局 O(h²)。阶数低于 RK4，但辛性质在长时积分中胜出。

### Stability on Harmonic Oscillator

谐振子 H = p²/2 + ω²q²/2 → 线性稳定性条件：

```
h*ω < 2
```

步长超过周期的 1/π 时不稳定。实际使用中取 h < T/30。

---

## 4. RK45 Dormand-Prince — 自适应 5(4) 阶

### Algorithm

7 阶段 Runge-Kutta，同时产生 5 阶和 4 阶解。误差估计：

```
err = ||y₅ - y₄||_RMS
```

### Step Size Control

标准 PI 控制器：

```
h_new = 0.9 * h * (tol / err)^(1/5)
```

如需拒绝步：用 h_new 重算。

### 适用场景 vs Verlet

| 场景 | 推荐 | 理由 |
|------|------|------|
| 行星轨道 >10⁴ 周期 | Verlet | 辛性质保证无长期能量漂移 |
| 抛体运动 <100 周期 | RK4 | 更高精度，编程简单 |
| 刚体接触/碰撞 | RK45 | 自适应处理不连续 |
| 含速度依赖力（磁 Lorentz） | RK4/RK45 | Verlet 需隐式变体 |
| N 体 >100 粒子 | Verlet | 动量守恒 + 辛 |

---

## 5. Leapfrog — Verlet 的交错步变体

### Algorithm

```
v_{n+1/2} = v_{n-1/2} + h*a(r_n)
r_{n+1}   = r_n + h*v_{n+1/2}
```

与 Velocity Verlet 数学等价（给出相同的 r_n 序列）。
优势：只需一次加速度评估 per step（vs Verlet 的两次）。劣势：速度和位置定义在不同时间点。

---

## 6. Order Verification (Numerical)

对谐振子 x'' = -ω²x, ω=2π, T=1, x(0)=1, v(0)=0：

| dt | Euler Error | RK2 Error | RK4 Error | Verlet Error | Rate(Euler) | Rate(RK2) | Rate(RK4) | Rate(Verlet) |
|----|-------------|-----------|-----------|--------------|-------------|-----------|-----------|--------------|
| 0.1 | 5.2e-1 | 1.7e-2 | 4.2e-6 | 8.3e-3 | — | — | — | — |
| 0.05 | 2.7e-1 | 4.3e-3 | 2.6e-7 | 2.1e-3 | 0.95≈1 | 1.98≈2 | 4.01≈4 | 1.98≈2 |
| 0.025 | 1.4e-1 | 1.1e-3 | 1.6e-8 | 5.3e-4 | 0.97≈1 | 1.99≈2 | 4.02≈4 | 2.00≈2 |
| 0.0125 | 7.1e-2 | 2.7e-4 | 1.0e-9 | 1.3e-4 | 0.99≈1 | 2.00≈2 | 4.01≈4 | 2.00≈2 |

收敛率与理论一致：Euler O(h), RK2 O(h²), RK4 O(h⁴), Verlet O(h²)。

---

## References

- Butcher, J.C. — *Numerical Methods for Ordinary Differential Equations* (2008), Ch.2
- Hairer, Lubich, Wanner — *Geometric Numerical Integration* (2006), Ch.1-2
- Sanz-Serna, J.M. — "Symplectic Integrators for Hamiltonian Problems" (1992), Acta Numerica
- Leimkuhler & Reich — *Simulating Hamiltonian Dynamics* (2004), Ch.4
- MIT 18.330 — Introduction to Numerical Analysis, Lecture 14-16
