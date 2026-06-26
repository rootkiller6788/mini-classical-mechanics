#!/usr/bin/env julia
# examples/projectile.jl — 抛体运动（含空气阻力）
# 对比解析解与数值解，验证积分器精度

include("../src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Projectile Motion — 抛体运动")
    println("="^60)

    # ---- 参数 ----
    v0_mag = 50.0         # 初速 [m/s]
    theta  = deg2rad(45)  # 发射角 45°
    z0     = 0.0          # 初始高度 [m]
    m      = 1.0          # 质量 [kg]
    t_end  = 10.0         # 模拟时间 [s]
    dt     = 0.001        # 步长 [s]

    v0 = Vec3(v0_mag * cos(theta), 0.0, v0_mag * sin(theta))
    r0 = Vec3(0.0, 0.0, z0)

    # ---- 无阻力（解析解存在） ----
    println("\n--- 1. No Drag (Analytic + Numerical) ---")

    g_val = G_EARTH
    function accel_no_drag(r, v, t)
        return Vec3(0.0, 0.0, -g_val)
    end

    # 解析解
    t_flight = projectile_flight_time(v0.z, z0; g=g_val)
    range_analytic = projectile_range(v0_mag, theta; g=g_val)
    println("Analytic: flight_time = $(round(t_flight, digits=3)) s, range = $(round(range_analytic, digits=3)) m")

    # 数值解 (Velocity Verlet)
    traj = solve_fixed_step(accel_no_drag, r0, v0, t_end, dt; method=:verlet, record_every=10)
    println("Numerical: n_steps = $(length(traj.ts)), final pos = $(traj.positions[end])")

    # 验证能量守恒
    function pot_no_drag(r) ; return g_val * m * r.z ; end
    energies, drift = energy_drift(traj, m, pot_no_drag)
    println("Energy drift (max relative): $(round(maximum(abs.(drift)), digits=10))")

    # ---- 含线性阻力 ----
    println("\n--- 2. Linear Drag (b = 0.5) ---")
    b = 0.5  # 阻力系数
    function accel_linear_drag(r, v, t)
        gravity = Vec3(0.0, 0.0, -g_val)
        drag = -b * v / m
        return gravity + drag
    end

    traj_drag = solve_fixed_step(accel_linear_drag, r0, v0, t_end, dt; method=:verlet, record_every=10)
    r_final = traj_drag.positions[end]
    println("With linear drag: final pos = $(r_final), range ≈ $(round(r_final.x, digits=3)) m")

    # ---- 含平方阻力 ----
    println("\n--- 3. Quadratic Drag (c = 0.01) ---")
    c = 0.01
    function accel_quad_drag(r, v, t)
        gravity = Vec3(0.0, 0.0, -g_val)
        speed = norm(v)
        drag = -c * speed * v / m
        return gravity + drag
    end

    traj_qdrag = solve_fixed_step(accel_quad_drag, r0, v0, t_end, dt; method=:verlet, record_every=10)
    r_final_q = traj_qdrag.positions[end]
    println("With quadratic drag: final pos = $(r_final_q), range ≈ $(round(r_final_q.x, digits=3)) m")

    # ---- 对比不同积分器 ----
    println("\n--- 4. Integrator Comparison (RK4 vs Verlet, no drag) ---")
    traj_rk4 = solve_fixed_step(accel_no_drag, r0, v0, t_end, dt; method=:rk4, record_every=10)
    diff = distance(traj.positions[end], traj_rk4.positions[end])
    println("|r_verlet - r_rk4| at t=$t_end: $diff")

    # ---- 自适应步长 ----
    println("\n--- 5. Adaptive RK45 ---")
    traj_adapt, accepted, rejected = solve_adaptive(accel_no_drag, r0, v0, t_end;
                                                    dt_init=0.01, tol=1e-6)
    println("Adaptive: accepted=$accepted, rejected=$rejected")
    println("Final: $(traj_adapt.positions[end])")

    println("\n✅ projectile.jl done.")
end

main()
