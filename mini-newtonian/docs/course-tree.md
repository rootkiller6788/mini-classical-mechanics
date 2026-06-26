# 课程树 — mini-newtonian

> 一棵树 = 一个学生的完整学习路径。从零到轨道力学。

## 知识树

```
                      F = ma （牛顿第二定律 ─ 一切从这里开始）
                            │
        ┌───────────────────┼───────────────────┐
        │                   │                   │
   [1] types.jl        [2] kinematics.jl    [3] forces.jl
   向量运算基础           运动的语言            力的来源
   ───────────           ───────────          ──────────
   Vec3(x,y,z)           r(t) 位移            引力 F=-GMm/r²
   dot/cross/norm        v(t) 速度            弹力 F=-kx
   Base 运算重载         a(t) 加速度          阻力 -bv, -cv|v|
   State 结构体          抛体运动              摩擦 μₖN
   Trajectory 容器       圆周运动              洛伦兹 q(E+v×B)
   ParticleODE           坐标变换              N体叠加
        │                   │                   │
        │                   │    ┌──────────────┘
        │                   │    │
        │              [4] integrators.jl
        │              把力变成运动 ─ 数值方法
        │              ────────────────────────
        │              Euler (教学, O(h))
        │              RK2 中点法 (O(h²))
        │              RK4 经典 (O(h⁴))
        │              RK45 自适应 (O(h⁵))
        │              Euler-Cromer 半隐式
        │              Velocity Verlet ★ 辛
        │              Leapfrog
        │              N体Verlet
        │                   │
        ├───────────────────┤
        │                   │
   [5] energy.jl        [6] momentum.jl
   守恒律(一)            守恒律(二)
   ──────────            ────────────
   动能 T=½mv²           线动量 p=mv
   势能 U(r)             角动量 L=r×p
   引力势能              质心 R_cm
   弹性势能              弹性碰撞 1D/3D
   能量漂移监测           非弹性碰撞
   机械能守恒验证         火箭方程 Δv
        │                   │
        └─────────┬─────────┘
                  │
             [7] constraints.jl
             受约束系统
             ──────────
             斜面 θ
             Atwood 机
             滑轮系统
             圆锥摆
             环形轨道
             弯道倾斜角
             电梯表观重量
                  │
        ┌─────────┼─────────┐
        │         │         │
   mini-      mini-      mini-
   lagrangian hamiltonian rigid-body
```

## 学习顺序（精确到文件）

```
Hour 1:  types.jl          Vec3 构造、加减、点积、叉积
Hour 2:  kinematics.jl     抛体解析解、圆周运动
Hour 3:  forces.jl         万有引力、胡克定律、阻力
Hour 4:  integrators.jl    Euler 法手写、理解 dt 的影响
Hour 5:  energy.jl         验证能量守恒、发现 Euler 漂移
Hour 6:  integrators.jl    Verlet 法、对比 RK4、理解辛积分
Hour 7:  momentum.jl       碰撞模拟、动量守恒
Hour 8:  constraints.jl    斜面、Atwood 机、环形轨道
Hour 9:  examples/         逐一运行 12 个示例
Hour10:  benchmark/        看 Verlet vs RK4 定量对比
Hour11:  tests/            学习物理量不变性测试
Hour12:  → mini-lagrangian 从 F=ma 到作用量原理
```

## 每个文件教你什么

| 文件 | 核心思想 | 学生应该能回答的问题 |
|------|---------|-------------------|
| `types.jl` | 类型系统是物理的骨架 | "为什么 Vec3 用 struct 不用 class？" |
| `kinematics.jl` | 运动不需要力 | "纯运动学能告诉我们什么？" |
| `forces.jl` | 力是加速度的原因 | "同一种运动可以由不同力产生吗？" |
| `integrators.jl` | 数值近似 = 物理近似的数学实现 | "为什么 Verlet 比 RK4 更适合轨道？" |
| `energy.jl` | 能量守恒是检验数值正确性的标准 | "如何用能量漂移判断代码 bug？" |
| `momentum.jl` | 对称性 → 守恒律 | "为什么碰撞问题用动量而不用力？" |
| `constraints.jl` | 约束力不独立——由几何决定 | "约束力和主动力有什么区别？" |
