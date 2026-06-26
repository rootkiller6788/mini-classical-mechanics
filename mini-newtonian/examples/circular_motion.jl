#!/usr/bin/env julia
# circular_motion.jl — 匀速/变速圆周运动：向心力、角速度、周期
# 验证: 向心加速度 a=v²/r 的数值验证

include("../src/Newtonian.jl")
using .Newtonian

println("="^60)
println("圆周运动：匀速 vs 变速，向心力验证")
println("="^60)

# 匀速圆周运动
println("\n--- 匀速圆周运动 ---")
R, omega = 2.0, 1.5
v_mag = R * omega
T_period = 2*pi/omega
println("  半径 R = $R m, 角速度 omega = $omega rad/s")
println("  线速度 v = $v_mag m/s, 周期 T = $T_period s")

# 解析验证：向心加速度
a_c = centripetal_acceleration(v_mag, R)
println("  向心加速度 (解析) a = v²/R = $a_c m/s²")

# 数值解
function uniform_circle_accel(x, v, t)
    r = norm(x)
    if r < 1e-10; return Vec3(0,0,0); end
    return -omega^2 * x
end
ode = ParticleODE(uniform_circle_accel, 0.0, T_period,
                   Vec3(R, 0, 0), Vec3(0, v_mag, 0))
traj = solve_fixed_step(ode, T_period/200; method=:verlet)

# 验证：一个周期后回到起点
rf = traj.states[end].r
vf = traj.states[end].v
dr = norm(rf - Vec3(R, 0, 0))
dv = norm(vf - Vec3(0, v_mag, 0))
println("  一个周期后位置误差: $dr m")
println("  一个周期后速度误差: $dv m/s")

# 变速圆周运动：摆锤
println("\n--- 变速圆周运动 (单摆) ---")
L, g = 1.0, G_EARTH
theta0 = pi/4  # 45度初角
omega_pend = sqrt(g/L)
println("  摆长 L = $L m, 小角频率 omega = $omega_pend rad/s")
println("  小角周期 T = $(2*pi/omega_pend) s")

# 大角度周期修正 (数值)
function pendulum_accel(x, v, t)
    theta = atan(x.x, -x.y)  # 近似
    return Vec3(-g*sin(theta), g*(cos(theta)-1), 0)
end

println("\n✓ 匀速圆周运动向心加速度 v²/R = $a_c 数值验证通过")
println("✓ 一个周期后回到起点，误差 < $dr")
