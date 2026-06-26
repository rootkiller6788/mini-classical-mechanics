# API 函数图谱 — mini-newtonian

> 不是罗列函数。是展示数据如何流动，以及"如果我想要 X，应该调用哪个函数，它依赖什么"。

## 数据流全景

```
用户输入:
  r₀, v₀ (初始条件)
  force_func (加速度函数)
  t_end, dt (积分参数)
    │
    ▼
ParticleODE(r₀, v₀, t_end, force_func)    ← types.jl
    │
    ├──→ solve_fixed_step(ode, dt; method)   ← integrators.jl
    │       │
    │       ├── method=:euler    → euler_step()
    │       ├── method=:verlet   → velocity_verlet_step()  ★推荐
    │       ├── method=:rk4      → rk4_step()
    │       └── method=:rk45     → rk45_step()  自适应
    │
    ▼
Trajectory (轨迹容器)
    │
    ├──→ energy_drift(traj, m, pot_func)      ← energy.jl
    │       验证数值质量: ΔE/E₀ 是否过大?
    │
    ├──→ angular_momentum(r, v, m)            ← momentum.jl
    │       验证: L 是否守恒? (中心力场下)
    │
    └──→ 绘图/分析
```

## 按"我想做 X"索引

### 我想计算抛体轨迹
```
projectile_motion(r₀, v₀, t; g)  →  解析解 (无空气阻力)
  或
ParticleODE(gravity_accel, ...)  →  solve_fixed_step  →  数值解 (可加阻力)
```

### 我想知道最大射程对应的抛射角
```
projectile_range(v₀, θ; g) 遍历 θ 找最大值
  (有空气阻力时不是 45°，需数值求解)
```

### 我想模拟弹簧-质量系统
```
hooke_force(r, r_eq, k) 作为加速度传入 ParticleODE
然后用 velocity_verlet_step 积分
最后 energy_drift 验证能量守恒
```

### 我想模拟太阳系
```
newton_gravity() + gravity_nbody() 作为加速度
nbody_verlet_step() 积分 (N体专用)
total_kinetic_energy() + total_gravitational_potential() 验证总能量
```

### 我想处理碰撞
```
elastic_collision_1d(m₁,v₁,m₂,v₂)  →  (v₁', v₂')
elastic_collision_3d(m₁,v₁,m₂,v₂,n)  →  (v₁', v₂')
  n 是碰撞法向
```

### 我想检查数值积分器是否可靠
```
energy_drift(traj, m, pot_func)  →  ΔE/E₀
  < 1e-4: 优秀
  < 1e-2: 可接受
  > 1e-2: 减小 dt 或换 Verlet
```

### 我想对比不同积分器
```
运行 benchmark/integrator_comparison.jl
  输出: Euler/Verlet/RK4 的能量漂移 % + 相位误差 % + 耗时
```

### 我想求斜面加速度
```
incline_acceleration(θ; mu, g)  →  a = g(sinθ - μ·cosθ)
angle_of_repose(μ_s)  →  θ_crit = arctan(μ_s)
```

### 我想求环形轨道最低速度
```
loop_the_loop_min_speed(R)  →  v = √(gR)
  低于此速 → 脱离轨道
```

### 我想计算火箭 Δv 预算
```
tsiolkovsky_delta_v(v_e, m₀, m_f)  →  Δv = v_e·ln(m₀/m_f)
  一级火箭: 算 v₁
  二级火箭: m₀₂ = m₁ - m_propellant₁, 算 v₂
  三级: ...总 Δv = v₁+v₂+v₃
```

## 函数调用链（精确）

### 积分一步 Verlet
```
velocity_verlet_step(accel_func, state, dt)
  │
  ├── accel_func(state.r, state.v, state.t)    ← 用户提供
  │      │
  │      └── e.g., newton_gravity(r, r_src, M, m)  ← forces.jl
  │                  │
  │                  └── normalize(dr)  ← types.jl
  │                       dot, cross, norm
  │
  ├── r_new = r + v·dt + ½·a_old·dt²
  ├── a_new = accel_func(r_new, v_half, t+dt)
  └── v_new = v + ½·(a_old + a_new)·dt
```

### 能量监测
```
energy_drift(traj, m, potential_func)
  │
  ├── for each state in traj.states:
  │     T = kinetic_energy(m, state.v)         ← energy.jl
  │     U = potential_func(state.r)            ← 用户提供
  │     E = T + U
  │
  └── return max(|E - E₀|/E₀)
```

## 类型层次

```
Vec3                    ← 基础向量 (不可变)
  ├── State             ← 运动状态 (可变): t, r::Vec3, v::Vec3
  ├── NBodyState        ← N体状态: t, masses, pos::Vector{Vec3}, vel::Vector{Vec3}
  └── Trajectory        ← 轨迹容器: states::Vector{State}

ParticleODE            ← ODE定义: accel::Function, t₀, t_end, r₀::Vec3, v₀::Vec3
```

```
Vec3(x,y,z) ──┬── dot(a,b) ── 标量
              ├── cross(a,b) ── Vec3
              ├── norm(v) ── Float64
              ├── normalize(v) ── Vec3
              └── distance(a,b) ── Float64
```

## 运动学 (`kinematics.jl`)

```
projectile_motion(r0,v0,t) ──→ 轨迹
projectile_range(v0,theta) ──→ 射程
projectile_flight_time(v0z,z0) ──→ 飞行时间
circular_position(center,r,omega,t) ──→ 位置
circular_velocity(r,omega,t) ──→ 速度
centripetal_acceleration(v,r) ──→ 加速度
galilean_position/velocity ──→ 坐标变换
```

## 力定律 (`forces.jl`)

```
newton_gravity(r_obj,r_src,M,m) ──→ 引力 Vec3
gravity_nbody(r_obj,m,sources) ──→ 多体引力
hooke_force(r,r_eq,k) ──→ 弹力
linear_drag(v,b) ──→ 线性阻力
quadratic_drag(v,c) ──→ 平方阻力
kinetic_friction(v,N,mu_k) ──→ 摩擦力
lorentz_force(q,E,v,B) ──→ 洛伦兹力
```

## ODE 积分器 (`integrators.jl`)

```
solve_fixed_step(ode,dt) ──┬── euler_step
                           ├── rk2_step
                           ├── rk4_step
                           ├── velocity_verlet_step  ← 辛, 推荐
                           └── euler_cromer_step

solve_adaptive(ode,tol) ──── rk45_step (Dormand-Prince)
```

## 能量 (`energy.jl`)

```
kinetic_energy(m,v) ──→ T = ½mv²
gravitational_potential(r_obj,r_src,M,m) ──→ U
elastic_potential(r,r_eq,k) ──→ U
energy_drift(traj,m,pot_func) ──→ ΔE 监测
```

## 动量与碰撞 (`momentum.jl`)

```
linear_momentum(m,v) ──→ p
angular_momentum(r,v,m) ──→ L
elastic_collision_1d(m1,v1,m2,v2) ──→ (v1',v2')
elastic_collision_3d(m1,v1,m2,v2,n) ──→ (v1',v2')
tsiolkovsky_delta_v(ve,m0,mf) ──→ Δv
center_of_mass(masses,pos) ──→ R_cm
```

## 约束系统 (`constraints.jl`)

```
incline_acceleration(theta,mu) ──→ a
atwood_machine(m1,m2) ──→ a
conical_pendulum_parameters(m,L,omega) ──→ (theta,T)
loop_the_loop_min_speed(R) ──→ v_min
banked_curve_angle(v,R) ──→ theta
```

## 数据流

```
初始条件 (r0,v0)
    │
    ├──→ 加速度函数 a(r,v,t)
    │       │
    │       ├── 力定律 → F
    │       └── F/m → a
    │
    ├──→ ODE 积分器
    │       │
    │       └──→ 轨迹 Trajectory
    │               │
    │               ├──→ 能量监测
    │               ├──→ 动量监测
    │               └──→ 可视化
    │
    └──→ 解析解（对照验证）
```
