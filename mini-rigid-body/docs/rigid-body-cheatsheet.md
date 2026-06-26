# Rigid Body Mechanics — Cheatsheet
> 参考: Goldstein Ch.4-5, Landau Vol.1 Ch.6

---

## 惯性张量 (Inertia Tensor)

### 定义 (Goldstein 5.3)
```
Ixx = Σ m(y²+z²)   Iyy = Σ m(x²+z²)   Izz = Σ m(x²+y²)
Ixy = -Σ mxy       Ixz = -Σ mxz       Iyz = -Σ myz
```

### 主轴 (Principal Axes)
对角化惯性矩阵: I_principal = R^T · I · R = diag(I₁, I₂, I₃)

### 平行轴定理 (Parallel Axis, Goldstein 5.18)
```
I_new = I_cm + M(d²𝟙 - d⊗d)
```

### 标准形状 (质心系)

| 形状 | I1 | I2 | I3 |
|------|----|----|----|
| 球体 (半径R) | 2/5 MR² | 2/5 MR² | 2/5 MR² |
| 球壳 (半径R) | 2/3 MR² | 2/3 MR² | 2/3 MR² |
| 圆柱 (R,H, z-对称) | M(3R²+H²)/12 | M(3R²+H²)/12 | ½ MR² |
| 长方体 (a×b×c) | M(b²+c²)/12 | M(a²+c²)/12 | M(a²+b²)/12 |
| 细杆 (长L, z轴) | ML²/12 | ML²/12 | 0 |
| 圆盘 (半径R, xy面) | ¼ MR² | ¼ MR² | ½ MR² |

---

## Euler 方程 (体坐标系)

### 无力矩 (Goldstein 5.44)
```
I₁ ω̇₁ = (I₂ - I₃) ω₂ ω₃
I₂ ω̇₂ = (I₃ - I₁) ω₃ ω₁
I₃ ω̇₃ = (I₁ - I₂) ω₁ ω₂
```

### 含外力矩
```
I₁ ω̇₁ - (I₂ - I₃) ω₂ ω₃ = N₁
I₂ ω̇₂ - (I₃ - I₁) ω₃ ω₁ = N₂
I₃ ω̇₃ - (I₁ - I₂) ω₁ ω₂ = N₃
```

### 运动常数
```
T = ½(I₁ω₁² + I₂ω₂² + I₃ω₃²)  = const
L² = I₁²ω₁² + I₂²ω₂² + I₃²ω₃² = const
```

---

## 欧拉角 (ZXZ Convention)

### 旋转矩阵 (Goldstein 4.46)
R(φ,θ,ψ) = Rz(φ) · Rx(θ) · Rz(ψ)

```
R = [cosφcosψ-sinφcosθsinψ   -cosφsinψ-sinφcosθcosψ   sinφsinθ]
    [sinφcosψ+cosφcosθsinψ   -sinφsinψ+cosφcosθcosψ  -cosφsinθ]
    [sinθsinψ                  sinθcosψ                 cosθ    ]
```

### 角速度 ↔ 欧拉角速率 (Goldstein 4.86-88)
```
ωx = φ̇sinθsinψ + θ̇cosψ
ωy = φ̇sinθcosψ − θ̇sinψ
ωz = φ̇cosθ + ψ̇

φ̇ = (ωx sinψ + ωy cosψ)/sinθ
θ̇ = ωx cosψ − ωy sinψ
ψ̇ = ωz − φ̇ cosθ
```

⚠ θ ≈ 0 时 sinθ → 0 → 万向锁

---

## 转动能量与角动量

### 转动动能 (Goldstein 5.9)
```
T = ½ ωᵀ·I·ω = ½(I₁ω₁² + I₂ω₂² + I₃ω₃²)   (主轴系)
```

### 角动量
```
L = I·ω        (体坐标系)
L_total = r_cm × (M v_cm) + I_cm·ω
```

### 功率
```
dT/dt = N·ω
```

---

## 网球拍定理 (Tennis Racket Theorem)

绕中间主轴的旋转 → 不稳定，周期性翻转

**稳定性判别** (I₁ ≥ I₂ ≥ I₃):
- 绕 I₁ (最大): ✅ 稳定
- 绕 I₂ (中间): ❌ 不稳定
- 绕 I₃ (最小): ✅ 稳定

**翻转增长率**:
```
λ = ω₀ √[((I₂-I₁)(I₃-I₂)) / (I₁I₃)]
T_flip ≈ log(1/δω₀) / λ
```

---

## 对称陀螺 (Lagrange Top)

### 有效势 (Goldstein 5.58)
```
V_eff(θ) = (Lz − L₃cosθ)²/(2I₁sin²θ) + Mgl cosθ
```

### 稳态进动条件
```
dV_eff/dθ = 0  →  (Lz − L₃cosθ)(Lz − L₃/cosθ) = Mgl I₁ sin⁴θ/cosθ
```

### 快速陀螺近似
```
φ̇ ≈ Mgl / (I₃ ω₃)
```

### Sleeping Top 稳定性 (Goldstein 5.63)
```
ω₃² > 4Mgl I₁ / I₃²  →  稳定
ω₃_crit = √(4Mgl I₁) / I₃
```

---

## 进动分类

| 类型 | φ̇ 行为 |
|------|--------|
| Monotonic | φ̇ 不改变符号 |
| Looping | φ̇ 在某些 θ 改变符号 |
| Cusped | θ 边界处 φ̇ = 0 |

---

## 积分器对比

| 方法 | 阶数 | 辛 | 运动常数保持 |
|------|------|-----|------------|
| Euler | 1 | ✗ | 差 |
| RK4 | 4 | ✗ | 短期好 |
| Velocity Verlet | 2 | ✓ | 长期优秀 |

---

## 陀螺仪公式

```
N_gyro = L × Ω_forced    (科里奥利力矩)
Ω_prec = N / L           (进动角速率)
```
