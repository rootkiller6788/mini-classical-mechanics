#!/usr/bin/env julia
# examples/pendulum.jl — 单摆、双摆、阻尼摆
# 展示相图、小角度近似 vs 大角度非线性

include("../src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Pendulum — 单摆与双摆")
    println("="^60)

    L  = 1.0       # 摆长 [m]
    m  = 1.0       # 质量 [kg]
    g  = G_EARTH   # 重力加速度
    omega0 = sqrt(g / L)  # 小角频率

    # ---- 1. 小角度单摆 (θ₀ = 5°) ----
    println("\n--- 1. Small-angle pendulum (5 deg) ---")
    theta0 = deg2rad(5.0)
    omega0_val = 0.0  # 初始角速度

    # 单摆加速度: θ'' = -(g/L)*sin(θ)
    function pendulum_accel(r, v, t)
        # r = (θ, 0, 0) 存储为 Vec3 的 x 分量
        theta = r.x
        alpha = -(g / L) * sin(theta)
        return Vec3(alpha, 0.0, 0.0)
    end

    r0 = Vec3(theta0, 0.0, 0.0)
    v0 = Vec3(omega0_val, 0.0, 0.0)
    t_end = 10.0
    dt = 0.001

    traj_small = solve_fixed_step(pendulum_accel, r0, v0, t_end, dt; method=:verlet, record_every=10)

    # 解析小角度周期: T = 2π√(L/g)
    T_small = 2π * sqrt(L / g)
    println("Small-angle analytic period: $(round(T_small, digits=4)) s")
    println("Simulated timesteps: $(length(traj_small.ts))")

    # 从轨迹估算周期（过零点法）
    zero_crossings = 0
    prev_theta = traj_small.positions[1].x
    for i in 2:length(traj_small.positions)
        curr_theta = traj_small.positions[i].x
        if prev_theta * curr_theta < 0
            zero_crossings += 1
        end
        prev_theta = curr_theta
    end
    T_sim = 2 * traj_small.ts[end] / zero_crossings
    println("Estimated period from simulation: $(round(T_sim, digits=4)) s")

    # ---- 2. 大角度单摆 (θ₀ = 60°) ----
    println("\n--- 2. Large-angle pendulum (60 deg) ---")
    theta0_large = deg2rad(60.0)
    r0_large = Vec3(theta0_large, 0.0, 0.0)
    traj_large = solve_fixed_step(pendulum_accel, r0_large, v0, t_end, dt; method=:verlet, record_every=10)

    # 大角度周期（椭圆积分近似）: T = T₀ * (1 + θ₀²/16 + ...)
    T_large_approx = T_small * (1 + theta0_large^2 / 16)
    println("Large-angle approx period (small-angle expansion): $(round(T_large_approx, digits=4)) s")
    println("Actual period will be longer due to nonlinearity")

    # ---- 3. 阻尼摆 ----
    println("\n--- 3. Damped pendulum (b = 0.2) ---")
    b = 0.2
    function damped_pendulum_accel(r, v, t)
        theta = r.x
        omega_val = v.x
        alpha = -(g / L) * sin(theta) - b * omega_val
        return Vec3(alpha, 0.0, 0.0)
    end

    traj_damped = solve_fixed_step(damped_pendulum_accel, r0_large, v0, 30.0, dt; method=:verlet, record_every=50)
    theta_final = traj_damped.positions[end].x
    println("After 30s: theta = $(round(rad2deg(theta_final), digits=4)) deg")
    println("n_steps = $(length(traj_damped.ts))")

    # ---- 4. 双摆（混沌演示） ----
    println("\n--- 4. Double pendulum (chaos demo) ---")
    L1, L2 = 1.0, 1.0
    m1, m2 = 1.0, 1.0

    # 双摆运动方程（拉格朗日力学导出，这里用牛顿力学降阶形式）
    function double_pendulum_accel(r, v, t)
        theta1, theta2 = r.x, r.y
        omega1, omega2 = v.x, v.y

        delta = theta2 - theta1
        denom1 = (m1 + m2) * L1 - m2 * L1 * cos(delta)^2
        denom2 = (L2 / L1) * denom1

        # 来自拉格朗日方程降阶
        sin_delta = sin(delta)
        cos_delta = cos(delta)

        alpha1 = (m2 * g * sin(theta2) * cos_delta
                  - m2 * L2 * omega2^2 * sin_delta
                  - (m1 + m2) * g * sin(theta1)
                  - m2 * L1 * omega1^2 * sin_delta * cos_delta) / denom1

        alpha2 = ((m1 + m2) * (g * sin(theta1) * cos_delta
                  + L1 * omega1^2 * sin_delta
                  - g * sin(theta2))
                  + m2 * L2 * omega2^2 * sin_delta * cos_delta) / denom2

        return Vec3(alpha1, alpha2, 0.0)
    end

    r0_double = Vec3(deg2rad(90.0), deg2rad(90.0), 0.0)  # 两个摆都从水平释放
    v0_double = Vec3(0.0, 0.0, 0.0)

    traj_double = solve_fixed_step(double_pendulum_accel, r0_double, v0_double, 20.0, 0.0005;
                                   method=:verlet, record_every=100)

    n = length(traj_double.ts)
    println("Double pendulum: n_recorded = $n, final t = $(round(traj_double.ts[end], digits=2)) s")
    println("Note: chaotic motion — tiny dt changes cause divergent trajectories")

    # ---- 5. 验证能量守恒（无阻尼单摆） ----
    println("\n--- 5. Energy conservation check ---")
    function pendulum_potential(r)
        theta = r.x
        return m * g * L * (1 - cos(theta))  # U = mgh, h = L(1-cosθ)
    end
    energies, drift = energy_drift(traj_large, m, pendulum_potential)
    println("Max relative energy drift (large angle): $(round(maximum(abs.(drift)), digits=10))")

    println("\n✅ pendulum.jl done.")
end

main()
