#!/usr/bin/env julia
# spinning_top.jl — 重力场中对称陀螺的进动与章动
# 参考: Goldstein Ch.5.7, Landau Ch.36

include("../src/RigidBody.jl"); using .RigidBody

function main()
    println("="^60)
    println("  Spinning Top — Precession & Nutation")
    println("="^60)

    # ==============================================
    # 陀螺参数: 均匀圆盘 + 细杆
    # ==============================================
    M_disk = 0.5   # 圆盘质量
    R_disk = 0.1   # 圆盘半径
    M_rod = 0.1    # 杆质量
    L_rod = 0.3    # 杆长（质心到支点距离）

    g = 9.81       # 重力加速度

    # 圆盘绕对称轴的转动惯量
    I3 = 0.5 * M_disk * R_disk^2           # 圆盘轴向
    I1_disk = 0.25 * M_disk * R_disk^2     # 圆盘横向

    # 杆对支点的转动惯量 (平行轴)
    I1_rod = (1/3) * M_rod * L_rod^2       # 杆绕端点

    # 总横向惯量 (绕支点)
    I1 = I1_disk + M_disk * L_rod^2 + I1_rod

    # 总质量 × 质心高度
    M_total = M_disk + M_rod
    l_cm = (M_disk * L_rod + M_rod * L_rod/2) / M_total

    println("\n陀螺参数:")
    println("  圆盘: M=$(M_disk) kg, R=$(R_disk) m")
    println("  杆:   M=$(M_rod) kg, L=$(L_rod) m")
    println("  I₁ (横向) = $(round(I1, 8)) kg·m²")
    println("  I₃ (轴向) = $(round(I3, 8)) kg·m²")
    println("  总质量 M = $(M_total) kg, 质心高 l = $(round(l_cm, 4)) m")

    # ==============================================
    # Sleeping Top 条件
    # ==============================================
    println("\n" * "-"^40)
    println("Sleeping Top (直立陀螺) 稳定性:")

    ω_crit = critical_spin_rate(M_total, g, l_cm, I1, I3)
    println("  临界自转角速度 ω₃_crit = $(round(ω_crit, 2)) rad/s")
    println("  即 ≈ $(round(ω_crit * 60 / (2π), 1)) RPM")

    # 测试不同自转速率
    for ω3_test in [0.5*ω_crit, ω_crit, 2*ω_crit]
        status, _ = sleeping_top_stability(M_total, g, l_cm, I1, I3, ω3_test)
        ratio = ω3_test / ω_crit
        println("  ω₃=$(round(ω3_test,1)) ($(round(ratio,2))× ω_crit) → $(status)")
    end

    # ==============================================
    # 稳态进动
    # ==============================================
    println("\n" * "-"^40)
    println("稳态进动分析:")

    ω3 = 2 * ω_crit  # 快速自转
    L3 = I3 * ω3
    println("  自转 ω₃ = $(round(ω3, 1)) rad/s, L₃ = $(round(L3, 4))")

    # 不同章动角下的进动速率比较
    for theta_deg in [15.0, 30.0, 45.0, 60.0]
        theta = theta_deg * π / 180
        Lz_needed = L3 * cos(theta)  # 稳态条件: Lz = L3 cosθ

        φdot_fast = fast_top_precession(M_total, g, l_cm, I3, ω3)
        φdot_exact = precession_rate_top(I1, Lz_needed, L3, theta)

        println("  θ=$(rpad(theta_deg,4))°: φ̇_fast=$(round(φdot_fast,4)), φ̇_exact=$(round(φdot_exact,4))")
    end

    # ==============================================
    # 有效势与章动
    # ==============================================
    println("\n" * "-"^40)
    println("有效势与章动:")

    Lz_val = L3 * cos(π/6) + 0.02  # 略偏离稳态值
    println("  Lz = $(round(Lz_val, 4)), L₃ = $(round(L3, 4))")

    thetas, Veff = symmetric_top_effective_potential(I1, I3, M_total, g, l_cm, Lz_val, L3)

    # 找到势能最小点和转折点
    Vmin, imin = findmin(Veff)
    theta_min = thetas[imin]
    println("  势能最小处 θ = $(round(theta_min * 180/π, 2))°")

    # 章动振幅: 给定能量略高于最小值
    E_level = Vmin + 0.02
    amp_range = nutation_amplitude_range(I1, I3, M_total, g, l_cm, E_level, Lz_val, L3)
    println("  章动振幅转折点: $(length(amp_range)) 个")

    if length(amp_range) >= 2
        theta_low = minimum(amp_range)
        theta_high = maximum(amp_range)
        println("  θ ∈ [$(round(theta_low*180/π,2))°, $(round(theta_high*180/π,2))°]")

        ptype = precession_type(I1, Lz_val, L3, amp_range)
        println("  进动类型: $ptype")
    end

    # 小振幅章动周期
    T_nut_small = nutation_period_small(I1, M_total, g, l_cm)
    println("  小振幅章动周期 ≈ $(round(T_nut_small, 4)) s")

    # ==============================================
    # Lagrange Top 完整模拟
    # ==============================================
    println("\n" * "-"^40)
    println("Lagrange Top 完整模拟:")

    theta0 = π / 6  # 初始倾斜 30°
    phi0 = 0.0
    psi0 = 0.0
    theta_dot0 = 0.05  # 微小程序动

    t_sim = 2.0
    dt_sim = 0.005
    times, thetas, phis, psis = simulate_lagrange_top(
        I1, I3, M_total, g, l_cm, theta0, phi0, psi0,
        theta_dot0, Lz_val, L3, t_sim, dt_sim)

    println("  模拟步数: $(length(times))")
    println("  θ 范围: [$(round(minimum(thetas)*180/π,2))°, $(round(maximum(thetas)*180/π,2))°]")
    println("  总进动角 Δφ = $(round(phis[end] * 180/π, 1))°")
    println("  总自转角 Δψ = $(round(psis[end] * 180/π, 1))°")

    # 每 50 步输出
    println("\n  状态快照 (每 100 步):")
    for i in 1:100:length(times)
        println("    t=$(rpad(round(times[i],3),6)): θ=$(rpad(round(thetas[i]*180/π,2),6))°  φ=$(rpad(round(phis[i]*180/π,1),6))°  ψ=$(rpad(round(psis[i]*180/π,1),6))°")
    end

    println("\n✅ spinning_top.jl 完成.")
end

main()
