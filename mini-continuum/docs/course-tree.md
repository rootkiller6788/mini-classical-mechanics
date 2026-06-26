# Course Tree — mini-continuum

```
连续介质力学 (Continuum Mechanics)
│
├── 运动学 (Kinematics)
│   ├── 小变形应变张量 ε = ½(∇u + ∇uᵀ)
│   ├── 变形梯度 F = ∂x/∂X
│   ├── 极分解 F = R·U = V·R
│   ├── Green-Lagrange应变 E = ½(FᵀF - I)
│   ├── Almansi应变 e = ½(I - F⁻ᵀF⁻¹)
│   └── 变形率 D = sym(ḞF⁻¹)
│
├── 应力分析 (Stress Analysis)
│   ├── Cauchy应力张量 σ (真实应力)
│   ├── 第一Piola-Kirchhoff应力 P = JσF⁻ᵀ
│   ├── 第二Piola-Kirchhoff应力 S = JF⁻¹σF⁻ᵀ
│   ├── Kirchhoff应力 τ = Jσ
│   ├── 应力不变量 I₁, J₂, J₃
│   └── 主应力与主方向
│
├── 守恒定律 (Conservation Laws)
│   ├── 质量守恒: ρ̇ + ρ∇·v = 0
│   ├── 动量守恒: ρa = ∇·σ + ρb
│   ├── 角动量守恒: σ = σᵀ
│   ├── 能量守恒: ρė = σ:D - ∇·q + ρr
│   └── Clausius-Duhem不等式
│
├── 弹性本构 (Elasticity)
│   ├── 各向同性: σ = λtr(ε)I + 2με
│   ├── 弹性常数互转 E,ν ↔ λ,μ ↔ K,G
│   ├── 平面应力/平面应变
│   ├── 正交各向异性 (9常数)
│   ├── 横观各向同性 (5常数)
│   ├── 非线性弹性 Neo-Hookean, Mooney-Rivlin
│   └── 热弹性 σ = λtr(ε)I + 2με - (3λ+2μ)αΔTI
│
├── 弹性波 (Elastic Waves)
│   ├── P波 v_p = √((λ+2μ)/ρ)
│   ├── S波 v_s = √(μ/ρ)
│   ├── Rayleigh表面波
│   ├── Love波 (层状半空间)
│   ├── Lamb波 (板)
│   ├── Snell定律, 反射/透射系数
│   └── 1D波动方程 FDTD
│
├── 梁理论 (Beam Theory)
│   ├── Euler-Bernoulli: EI·wⁱᵛ = q(x)
│   ├── Timoshenko梁 (计入剪切)
│   ├── 临界屈曲 P_cr = π²EI/(KL)²
│   ├── 梁振动 ω_n = (nπ/L)²√(EI/ρA)
│   ├── Winkler弹性地基梁
│   └── 连续梁 (三弯矩方程)
│
├── 板壳理论 (Plates & Shells)
│   ├── Kirchhoff板 D∇⁴w = q
│   ├── Mindlin板 (计入剪切)
│   ├── 圆板 (简支/固支)
│   ├── 矩形板振动
│   ├── 薄膜振动
│   └── 圆柱壳
│
├── 强度与失效 (Strength & Failure)
│   ├── von Mises, Tresca, Rankine
│   ├── Mohr-Coulomb, Drucker-Prager
│   ├── 疲劳 Basquin, Coffin-Manson
│   ├── 平均应力修正 Goodman, Gerber
│   └── Palmgren-Miner累积损伤
│
├── 断裂力学 (Fracture Mechanics)
│   ├── K因子 K_I = βσ√(πa)
│   ├── 能量释放率 G = K²/E
│   ├── J积分 (弹塑性)
│   ├── Paris裂纹扩展律
│   └── Weibull统计断裂
│
├── 塑性力学 (Plasticity)
│   ├── J₂流动理论, 径向返回
│   ├── 等向/随动硬化
│   ├── Drucker-Prager, Cam-Clay
│   └── Mohr-Coulomb塑性
│
├── 粘弹性 (Viscoelasticity)
│   ├── Maxwell, Kelvin-Voigt, SLS
│   ├── Prony级数, Boltzmann叠加
│   ├── 复模量, DMA
│   └── WLF时温叠加
│
├── 接触力学 (Contact)
│   ├── Hertz弹性接触
│   ├── JKR/DMT粘着理论
│   ├── Coulomb/Stribeck摩擦
│   └── Mindlin切向接触
│
├── 复合材料 (Composites)
│   ├── 混合律, Halpin-Tsai
│   ├── 层合板ABD矩阵
│   └── Tsai-Hill/Tsai-Wu准则
│
├── 微观力学 (Micromechanics)
│   ├── Eshelby夹杂
│   ├── Mori-Tanaka均质化
│   ├── 位错力学
│   └── Hall-Petch, Taylor硬化
│
├── 多孔弹性 (Poroelasticity)
│   ├── Biot理论
│   ├── Terzaghi固结
│   └── Gassmann方程
│
└── 结构稳定性 (Stability)
    ├── 屈曲, 后屈曲
    ├── Perry-Robertson公式
    ├── 分岔理论
    └── 动力屈曲
```

**前置依赖:**
- mini-newtonian: F=ma, 向量运算, ODE积分器

**后续路径:**
- mini-fluid-dynamics: Navier-Stokes方程
- mini-solid-state: 晶格动力学
- mini-computational-physics: FEM/FDM高级求解器
