# Course Alignment — 课程对照

## MIT 8.012 Physics I: Classical Mechanics

| 8.012 Topic | mini-newtonian | 实现 |
|-------------|---------------|------|
| Ch.1-2 Kinematics | `kinematics.jl` | 运动学解析解 + 轨迹计算 |
| Ch.3-4 Newton's Laws | `forces.jl` | 所有力定律的手写实现 |
| Ch.5 Work & Energy | `energy.jl` | 动能/势能/功, 能量守恒验证 |
| Ch.6 Momentum & Collisions | `momentum.jl` | 线动量/角动量/碰撞/火箭方程 |
| Ch.7 Numerical Methods | `integrators.jl` | Euler, RK2, RK4, Verlet, RK45 |
| Ch.8 Central Force | `examples/central_force.jl` | 开普勒问题, 轨道力学 |
| Ch.9 Harmonic Oscillator | `examples/oscillator.jl` | SHO/阻尼/受迫/共振/耦合 |
| Ch.10 Rigid Body | (mini-rigid-body) | 独立子模块 |
| Ch.11 Lagrangian | (mini-lagrangian) | 独立子模块 |

## Goldstein Classical Mechanics

| Goldstein | mini-newtonian |
|-----------|---------------|
| Ch.1 Survey | `kinematics.jl` |
| Ch.2 Variational Principles | (mini-variational-principle) |
| Ch.3 Central Force | `examples/central_force.jl` |
| Ch.4 Kinematics of Rigid Body | (mini-rigid-body) |
| Ch.5 Rigid Body Equations | (mini-rigid-body) |
| Ch.6 Oscillations | `examples/oscillator.jl` |
| Ch.8 Hamiltonian Mechanics | (mini-hamiltonian) |

## Numerical Methods

| Method | Implementation | Order | Properties |
|--------|---------------|-------|------------|
| Euler | `euler_step` | 1 | 显式, 不稳 |
| Euler-Cromer | `euler_cromer_step` | 1 | 半隐式, 能量好 |
| Midpoint (RK2) | `rk2_step` | 2 | 改进精度 |
| RK4 | `rk4_step` | 4 | 经典通用 |
| RK45 (Dormand-Prince) | `rk45_step` | 5(4) | 自适应 |
| Velocity Verlet | `velocity_verlet_step` | 2 | 辛, 守恒系统首选 |
| Leapfrog | `leapfrog_step` | 2 | 辛, 交错步 |
| N-body Verlet | `nbody_verlet_step` | 2 | 多体辛积分 |
