#!/usr/bin/env julia
# examples/constrained.jl — 约束系统
# 斜面、滑轮、Atwood机、弯道、圆锥摆

include("../src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Constrained Systems — 约束系统")
    println("="^60)

    g = G_EARTH

    # ---- 1. 斜面 ----
    println("\n--- 1. Inclined Plane ---")
    theta30 = deg2rad(30)
    a_no_friction = incline_acceleration(theta30)
    a_with_friction = incline_acceleration(theta30; mu=0.3)
    println("30 deg incline: a_no_friction=$(round(a_no_friction,digits=2)) m/s²")
    println("30 deg + μ=0.3: a=$(round(a_with_friction,digits=2)) m/s²")

    # 静止角
    mu_s_test = 0.577  # ≈ tan(30°)
    repose = rad2deg(angle_of_repose(mu_s_test))
    println("μ_s=$mu_s_test → angle of repose = $(round(repose,digits=1))° (expect ~30°)")

    # ---- 2. Atwood 机 ----
    println("\n--- 2. Atwood Machine ---")
    for (m1, m2) in [(5.0, 3.0), (2.0, 2.0), (1.0, 4.0)]
        a1, a2, T = atwood_machine(m1, m2)
        direction = a1 > 0 ? "m₁ down" : a1 < 0 ? "m₂ down" : "balanced"
        println("  m₁=$m1 kg, m₂=$m2 kg → a₁=$(round(a1,digits=2)) m/s², T=$(round(T,digits=1)) N ($direction)")
    end

    # ---- 3. 圆锥摆 ----
    println("\n--- 3. Conical Pendulum ---")
    L = 2.0
    m = 1.0
    omega_crit = conical_pendulum_critical_omega(L)
    println("Critical ω for L=$L m: $(round(omega_crit,digits=3)) rad/s = $(round(omega_crit*60/(2π),digits=1)) RPM")

    for omega in [omega_crit*1.5, omega_crit*2.0, omega_crit*3.0]
        theta, T = conical_pendulum_parameters(m, L, omega)
        if !isnan(theta)
            println("  ω=$(round(omega,digits=2)) rad/s → θ=$(round(rad2deg(theta),digits=2))°, T=$(round(T,digits=2)) N")
        end
    end

    # ---- 4. 环形轨道 ----
    println("\n--- 4. Loop-the-Loop ---")
    R_loop = 10.0  # 环半径 10m
    v_min = loop_the_loop_min_speed(R_loop)
    println("Loop radius: $R_loop m")
    println("Minimum speed at top: $(round(v_min,digits=2)) m/s = $(round(v_min*3.6,digits=1)) km/h")

    # 不同速度下的法向力
    for frac in [0.8, 1.0, 1.5, 2.0]
        v = frac * v_min
        is_detached, N = detachment_condition(70.0, v, R_loop, 0.0)  # 70kg 的人/车，顶点 θ=0
        status = is_detached ? "DETACHED!" : "safe"
        println("  v=$(round(v,digits=2)) m/s ($(round(frac*100))% v_min) → N=$(round(N,digits=1)) N ($status)")
    end

    # ---- 5. 斜面+滑轮 ----
    println("\n--- 5. Incline-Pulley System ---")
    for (m1, m2, theta_deg) in [(2.0, 1.0, 30), (1.0, 2.0, 30), (3.0, 1.5, 45)]
        theta = deg2rad(theta_deg)
        a, T = incline_pulley_system(m1, m2, theta)
        dir = incline_pulley_direction(m1, m2, theta)
        dir_str = dir == :m1_up ? "m₁ up incline" : dir == :m1_down ? "m₁ down incline" : "static"
        println("  m₁=$m1 on $(theta_deg)°, m₂=$m2 → a=$(round(a,digits=2)) m/s², T=$(round(T,digits=1)) N ($dir_str)")
    end

    # ---- 6. 弯道 ----
    println("\n--- 6. Banked Curves ---")
    R_curve = 200.0  # 弯道半径 [m]
    for v_kmh in [60.0, 80.0, 100.0, 120.0]
        v = v_kmh / 3.6
        theta_bank = rad2deg(banked_curve_angle(v, R_curve))
        v_max = max_safe_speed_banked_curve(R_curve, deg2rad(10), 0.7) * 3.6  # 10° 倾斜 + μ=0.7
        println("  R=$R_curve m: v=$(round(v_kmh)) km/h → ideal bank=$(round(theta_bank,digits=1))° | max safe (10° bank, μ=0.7)=$(round(v_max,digits=0)) km/h")
    end

    # ---- 7. 电梯 ----
    println("\n--- 7. Elevator Apparent Weight ---")
    m_person = 70.0
    for a_elev in [-g*0.3, 0.0, g*0.2, g*0.5, g*1.0]
        apparent = apparent_weight_in_elevator(m_person, a_elev) / (m_person * g) * 100
        desc = a_elev > 0 ? "accelerating UP" : a_elev < 0 ? "accelerating DOWN" : "constant velocity"
        println("  a=$(round(a_elev,digits=1)) m/s² ($desc) → weight = $(round(apparent,digits=1))% normal")
    end

    println("\n✅ constrained.jl done.")
end

main()
