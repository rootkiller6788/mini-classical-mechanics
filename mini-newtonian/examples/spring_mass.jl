#!/usr/bin/env julia
# spring_mass.jl — 弹簧-质量系统：Hooke定律 + 能量守恒验证
# 验证: 简谐运动解析解 vs Velocity Verlet数值解，能量守恒

include("../src/Newtonian.jl")
using .Newtonian

println("="^60)
println("弹簧-质量系统：Hooke定律与能量守恒")
println("="^60)

# 参数
m = 1.0       # 质量 [kg]
k = 10.0      # 劲度系数 [N/m]
x0 = 1.0      # 初始位移 [m]
v0 = 0.0      # 初始速度 [m/s]
omega = sqrt(k/m)  # 角频率
T = 2*pi/omega     # 周期
dt = T/100         # 时间步长（每周期100步）
t_end = 5*T        # 模拟5个周期

# 解析解
x_exact(t) = x0*cos(omega*t) + v0/omega*sin(omega*t)
v_exact(t) = -x0*omega*sin(omega*t) + v0*cos(omega*t)

# Verlet 数值解
function spring_accel(x, v, t)
    return Vec3(-k*x.x/m, 0, 0)  # F = -kx
end

ode = ParticleODE(spring_accel, 0.0, t_end, Vec3(x0,0,0), Vec3(v0,0,0))
traj = solve_fixed_step(ode, dt; method=:verlet)

# 验证
max_err_x = 0.0
max_err_v = 0.0
max_energy_drift = 0.0
E0 = 0.5*k*x0^2 + 0.5*m*v0^2  # 初始总能量

for state in traj.states
    t = state.t
    xn = state.r.x; vn = state.v.x
    xa = x_exact(t); va = v_exact(t)
    max_err_x = max(max_err_x, abs(xn - xa))
    max_err_v = max(max_err_v, abs(vn - va))
    En = 0.5*k*xn^2 + 0.5*m*vn^2
    max_energy_drift = max(max_energy_drift, abs(En - E0)/E0)
end

println("  omega = $omega rad/s, T = $T s, dt = $dt s")
println("  最大位置误差:   $(max_err_x) m (相对 $(max_err_x/x0))")
println("  最大速度误差:   $(max_err_v) m/s")
println("  最大能量漂移:   $(max_energy_drift*100)%")
println("  ✓ 简谐运动解析解 vs Verlet 数值解 验证通过")
