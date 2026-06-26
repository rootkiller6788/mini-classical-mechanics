# Calculus of Variations — Cheatsheet
> 参考: Gelfand & Fomin, Goldstein Ch.2

---

## Euler-Lagrange 方程

### 最简泛函: J[y] = ∫_a^b F(x, y, y') dx
```
∂F/∂y - d/dx(∂F/∂y') = 0
```

### Beltrami 恒等式 (F不显含x)
```
F - y' ∂F/∂y' = const
```

### 循环坐标 (F不显含y)
```
∂F/∂y' = const (第一积分)
```

### 多变量: J[y₁,...,y_n]
```
∂F/∂y_i - d/dx(∂F/∂y'_i) = 0,  i=1..n
```

### 高阶导数: F(x,y,y',y'',...,y^{(m)})
```
Σ_{k=0}^m (-1)^k d^k/dx^k (∂F/∂y^{(k)}) = 0
```

### 场论: J[u] = ∬ F(x,y,u,u_x,u_y) dxdy
```
∂F/∂u - ∂/∂x(∂F/∂u_x) - ∂/∂y(∂F/∂u_y) = 0
```

## 边界条件

| 端点条件 | 公式 |
|---------|------|
| 固定端点 | y(a)=y_a, y(b)=y_b |
| 自由端点 | ∂F/∂y'|_{end} = 0 |
| 横截条件 | F+(φ'-y')∂F/∂y' = 0 (端点沿φ移动) |
| 自然边界 | ∂F/∂y' = 0 |

## 二级变分与充分条件

### Legendre 条件 (必要条件)
```
∂²F/∂y'² ≥ 0 (弱极小)
```

### Jacobi 条件
共轭点方程: d/dx(P du/dx) + Q u = 0
无共轭点在(a,b)内

### Weierstrass 条件 (强极小)
```
E(x,y,p_opt,p) = F(x,y,p) - F(x,y,p_opt) - (p-p_opt)∂F/∂y'|_{opt} ≥ 0
```

### 充分条件: Legendre + Jacobi + Weierstrass 全部满足

## 约束变分

### 等周约束: ∫G dx = c
```
Lagrange乘子: F* = F + λG
解 EL 方程 + G-约束方程 联立
```

### 可变端点
横截条件自动确定最优端点位置

### 不等式约束 (障碍)
```
y(x) ≥ ψ(x) → 自由区域: EL方程; 接触区域: y=ψ
```

## 直接法

### Ritz 法
```
y ≈ Σ c_i φ_i(x)  (φ_i 满足边界条件)
求解: ∂J/∂c_i = 0 → 线性方程组
```

### Galerkin 法
```
∫ EL(y) φ_i dx = 0  (加权残差)
```

## 著名变分问题

| 问题 | 泛函 | 极值曲线 |
|------|------|---------|
| 最速降线 | ∫√(1+y'²)/√(2gy) dx | 摆线 |
| 悬链线 | ∫y√(1+y'²) dx | 双曲余弦 |
| 极小曲面 | ∬√(1+u_x²+u_y²) dxdy | 悬链面 |
| Dido问题 | max面积, 固定周长 | 圆 |
| 测地线 | ∫√(1+y'²) dx | 直线(平面)/大圆(球面) |
