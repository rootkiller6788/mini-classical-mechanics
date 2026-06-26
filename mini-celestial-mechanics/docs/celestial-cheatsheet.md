# Celestial Mechanics — Cheatsheet
> 参考: Goldstein Ch.3, Murray & Dermott, Vallado

---

## Kepler 方程

| 类型 | 方程 | 求解 |
|------|------|------|
| 椭圆 (e<1) | M = E - e sin E | Newton迭代: E_{n+1}=E_n-(E_n-e sinE_n-M)/(1-e cosE_n) |
| 双曲线 (e>1) | M_h = e sinh H - H | Newton: H_{n+1}=H_n-(e sinhH_n-H_n-M_h)/(e coshH_n-1) |
| 抛物线 (e=1) | Barker方程 | tan(ν/2) = B - 1/B, B=(3M_p/2+√((3M_p/2)²+1))^(1/3) |

## 近点角互转

```
tan(ν/2) = √((1+e)/(1-e)) tan(E/2)    真↔偏
tan(E/2) = √((1-e)/(1+e)) tan(ν/2)    偏↔真
M = E - e sin E                        偏→平
```

## 轨道根数 ↔ 状态向量

**根数→状态**: r = R_z(Ω)R_x(i)R_z(ω) · (a(1-e²)/(1+e cosν)) [cosν, sinν, 0]^T

**状态→根数**: h=r×v, e=(v×h)/μ - r̂, a=1/(2/r - v²/μ)

## 二体公式

```
Vis-viva:    v² = μ(2/r - 1/a)
周期:        T = 2π√(a³/μ)
平均角速度:  n = √(μ/a³)
半通径:      p = a(1-e²)
近心距:      r_p = a(1-e), 远心距: r_a = a(1+e)
逃逸速度:    v_esc = √(2μ/r)
```

## J2 摄动 (长期)

```
Ω̇ = -(3/2)n J₂(R/p)² cos i
ω̇ = (3/4)n J₂(R/p)² (5cos²i - 1)
```

太阳同步条件: Ω̇ = 2π/365.25 day⁻¹ → cos i = -Ω̇_sun / (1.5 n J₂(R/p)²)

## 三体问题

Lagrange 点: 共线点 L1-L3 解 x - (1-μ)(x+μ)/|x+μ|³ - μ(x-1+μ)/|x-1+μ|³ = 0
三角点 L4,L5: (1/2-μ, ±√3/2)

Jacobi 常数: C_J = 2Ω(x,y,z) - v²

## 影响球

```
SOI:   r_SOI = a (m/M)^(2/5)
Hill:  r_H = a (m/3M)^(1/3)
Roche: d = 1.26 R_1 (ρ₁/ρ₂)^(1/3) (刚体)
       d = 2.44 R_1 (ρ₁/ρ₂)^(1/3) (流体)
```

## 轨道转移

```
Hohmann: Δv₁ = |v_trans_peri - v_circ1|, Δv₂ = |v_circ2 - v_trans_apo|
双椭圆: 经由 r_b 的 3-脉冲转移
平面变化: Δv = 2v sin(Δi/2)
引力辅助: ΔV_max = 2V_planet |sin(δ/2)|
```

## 发射窗口

相位角: φ = π - n_target T_transfer
会合周期: 1/T_syn = |1/T₁ - 1/T₂|
