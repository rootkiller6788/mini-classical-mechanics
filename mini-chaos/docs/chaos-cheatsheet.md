# Chaos & Nonlinear Dynamics — Cheatsheet
> 参考: Strogatz, Goldstein Ch.11, May 1976

---

## 连续系统 (Flows)

### Lorenz
```
ẋ = σ(y-x), ẏ = x(ρ-z)-y, ż = xy-βz
σ=10, ρ=28, β=8/3 → 混沌
Lyapunov: λ₁≈0.9, λ₂≈0, λ₃≈-14.6
```

### Rössler
```
ẋ = -y-z, ẏ = x+ay, ż = b+z(x-c)
a=0.2, b=0.2, c=5.7 → 混沌
```

### Chua 电路
```
ẋ = α(y-x-f(x)), ẏ = x-y+z, ż = -βy
f(x) = m₁x + ½(m₀-m₁)(|x+1|-|x-1|)
α=15.6, β=28, m₀=-1.14, m₁=-0.714
```

### Duffing
```
ẍ + δẋ - βx + αx³ = γ cos ωt
δ=0.2, β=1, α=1, γ=0.3, ω=1.0 → 混沌
```

## 离散映射 (Maps)

| 映射 | 公式 | 混沌参数 |
|------|------|---------|
| Logistic | x_{n+1}=r x_n(1-x_n) | r≈3.57→4 |
| Hénon | x_{n+1}=1-ax²+y, y_{n+1}=bx | a=1.4, b=0.3 |
| Standard | p'=p+K sinθ, θ'=θ+p' | K>0.9716 |
| Circle | θ'=θ+Ω-(K/2π)sin(2πθ) | K>1 |
| Tent | x'=μ min(x,1-x) | μ=2 |
| Ikeda | z'=A+Bz exp(iκ-ip/(1+|z|²)) | 光学腔 |

## Lyapunov 指数

```
λ = lim_{t→∞} (1/t) ln(|δx(t)|/|δx(0)|)

λ>0: 混沌 (对初值敏感依赖)
λ=0: 周期/准周期
λ<0: 稳定不动点
```

全谱: λ₁≥λ₂≥...≥λ_n
Kaplan-Yorke 维数: D_KY = j + Σ_{i≤j} λ_i/|λ_{j+1}|

## 分岔

```
倍周期分岔 → 混沌 (Feigenbaum 路径)
准周期 → 混沌 (Ruelle-Takens 路径)
阵发性 → 混沌 (Pomeau-Manneville 路径)
```

### Feigenbaum 常数
```
δ = 4.669201609... (分岔间距比)
α = 2.502907875... (分岔宽度比)
```

## 分形

Mandelbrot: z_{n+1} = z_n² + c, z₀=0
Julia: z_{n+1} = z_n² + c (c固定, z₀可变)
Newton分形: z_{n+1} = z_n - f(z)/f'(z)

盒计数维数: D = lim_{ε→0} log N(ε)/log(1/ε)
