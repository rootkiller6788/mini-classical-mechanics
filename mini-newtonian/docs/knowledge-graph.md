# 知识图谱 — mini-newtonian

> 不只是罗列概念。要展示依赖关系、学习路径、以及"为什么这个顺序"。

---

## 概念依赖图（拓扑排序）

```
数学层           Vec3 (向量运算)
                  │
                  ├── norm, normalize, dot, cross
                  │   这些是所有物理量的基础类型
                  │
物理层    ┌───────┼───────┐
         │       │       │
    kinematics  forces  (无依赖，都是Vec3操作)
         │       │
         │ 运动描述   │ 力的来源
         │ (r,v,a)   │ (F=Gm₁m₂/r², F=-kx, ...)
         │       │
         └───┬───┘
             │  两者结合 → ODE问题: d²r/dt² = F/m
             │
      integrators
             │  Euler/RK2/RK4/RK45/Verlet/Leapfrog
             │  把力变成轨迹
             │
      ┌──────┼──────┐
      │      │      │
   energy  momentum  constraints
      │      │      │
   ½mv²   p=mv   斜面/滑轮/
   U(r)   L=r×p  弯道/圆锥摆
      │      │
      └──┬───┘
         │
    守恒律验证 (tests/)
    能量漂移监测
    动量守恒检查
    角动量守恒检查
```

## 为什么这个顺序？

每个文件只能依赖它上面的文件。这是刻意设计——迫使代码模块化、零循环依赖：

| 文件 | 可以 import | 不能 import | 理由 |
|------|-----------|------------|------|
| `types.jl` | (无) | 一切 | 基础类型，必须先加载 |
| `kinematics.jl` | types | forces, integrators | 运动学不依赖力 |
| `forces.jl` | types | integrators, energy | 力定律是纯函数 |
| `integrators.jl` | types, forces? (通过accel函数) | energy, momentum | 积分器只解ODE |
| `energy.jl` | types, kinematics? | momentum | 能量计算用运动状态 |
| `momentum.jl` | types | energy | 动量计算用运动状态 |
| `constraints.jl` | types, kinematics | — | 约束系统较独立 |

**实际加载顺序**（`Newtonian.jl` 中）：types → kinematics → forces → integrators → energy → momentum → constraints

## 概念精确定义

### Vec3 → 为什么手写不 import？

```julia
struct Vec3  # 不可变，栈分配，零开销
    x::Float64; y::Float64; z::Float64
end
```

**不用 `StaticArrays.SVector{3,Float64}` 的原因**：
- 教学目的：让学生理解内存布局
- 零依赖原则：复制目录即可运行
- 性能等同：Julia 的 immutable struct 同样栈分配

### State → 为什么可变？

```julia
mutable struct State
    t::Float64; r::Vec3; v::Vec3
end
```

轨迹记录需要 O(10⁴-10⁶) 个状态。不可变 struct 每次更新都分配新对象 → GC 压力。可变 struct 原地更新 → 堆分配一次。

### 力 → 为什么是函数指针不是类型？

```julia
# 力 = 函数: (Vec3, Vec3, Float64) → Vec3
# 而非 struct Force{...}
```

这样任何 Julia 函数都可以作为力传入积分器。无需提前注册、无需继承体系。函数式设计 = 最大灵活性。

### 积分器 → 为什么 Verlet 是"推荐"？

基于 `benchmark/integrator_comparison.jl` 的定量数据：

| 性质 | Euler | RK4 | Verlet |
|------|-------|-----|--------|
| 能量漂移（1000周期） | 发散 | ∝t 单调增长 | 有界振荡 |
| 相空间体积 | 不守恒 | 不守恒 | 守恒（辛） |
| 时间反演 | 不可逆 | 不可逆 | 可逆（±1e-6） |
| 每步力评估 | 1次 | 4次 | 2次 |
| 适用场景 | 仅教学 | 非保守系统 | **保守系统** |

---

## 学习路径：从零到轨道力学

```
第1小时:   types.jl         理解Vec3，写向量加减
第2小时:   kinematics.jl    理解位移→速度→加速度的关系
第3小时:   forces.jl        手写万有引力和胡克定律
第4小时:   integrators.jl   理解Euler法，亲手实现
第5小时:   energy.jl        验证能量守恒，发现Euler法的能量漂移
第6小时:   momentum.jl      碰撞模拟，动量守恒
第7小时:   constraints.jl   斜面、滑轮、圆锥摆
第8小时:   跑全部examples   从抛体到N体太阳系
```

---

## 跨模块连接：具体接口

| 从这个函数 | 到那个模块 | 接口 |
|-----------|----------|------|
| `Vec3` | mini-lagrangian | 广义坐标 q 的类型 |
| `State` | mini-hamiltonian | 相空间点 (q,p) 的类型 |
| `angular_momentum()` | mini-rigid-body | 角动量 → 惯量张量 → Euler 方程 |
| `velocity_verlet_step()` | mini-chaos | 积分器用于 Lorenz 系统 |
| `newton_gravity()` | mini-celestial | N体引力用于 Kepler 轨道 |
| `energy_drift()` | mini-continuum | 能量监测用于弹性波 |
