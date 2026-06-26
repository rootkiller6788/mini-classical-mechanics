#!/usr/bin/env julia
# tennis_racket_theorem.jl — 网球拍定理 / Dzhanibekov 效应
# 参考: Goldstein Ch.5.6, Landau Ch.37
# 演示: 绕中间主轴的旋转不稳定，导致周期性180°翻转

include("../src/RigidBody.jl"); using .RigidBody

function main()
    println("="^60)
    println("  Tennis Racket Theorem — Dzhanibekov Effect")
    println("="^60)

    # 构造非对称刚体: I₁ > I₂ > I₃
    # 长方体 a=3, b=2, c=1 (质量=1)
    I = cuboid_inertia(1.0, 3.0, 2.0, 1.0)
    pm = principal_moments(I)
    println("\n惯性张量 (长方体 3×2×1):")
    println("  I₁ = $(round(pm[1],4))  (最大 — 稳定)")
    println("  I₂ = $(round(pm[2],4))  (中间 — 不稳定!)")
    println("  I₃ = $(round(pm[3],4))  (最小 — 稳定)")

    # ==============================================
    # 实验 1: 绕中间轴旋转 ≈ 稳定分析
    # ==============================================
    println("\n" * "-"^40)
    println("实验 1: 绕中间主轴 (I₂) 旋转 + 微小扰动")

    ω_mag = 1.0
    δω = 0.01  # 微小扰动

    # 绕中间轴的稳定性
    _, growth_rate = axis_stability(I, 2, ω_mag)
    println("  线性增长率: λ = $(round(growth_rate, 4))")
    println("  e-折叠时间: τ = $(round(1/growth_rate, 4)) s")

    # 翻转周期估计
    T_flip = flipping_period_estimate(I, ω_mag, δω)
    println("  翻转周期估计: T ≈ $(round(T_flip, 3)) s")

    # 模拟: 绕 I₂, 给微小扰动
    w0 = [δω, ω_mag, δω]  # 主要绕 y 轴 (I₂)
    T0_init, L2_init = motion_constants(I, w0)
    println("\n  初始 ω = $(round.(w0, 4))")
    println("  初始 T = $(round(T0_init, 6)), L² = $(round(L2_init, 6))")

    t_end = 15.0
    dt = 0.01
    traj_mid = simulate_free_rigid_body(I, w0, t_end, dt; method=:rk4)

    # 每隔一定步数输出 ω
    println("\n  角速度演化 (每 100 步):")
    for i in 1:100:min(length(traj_mid), 1501)
        t = (i - 1) * dt
        w = traj_mid[i]
        println("    t=$(rpad(round(t,2),5)): ω = $(rpad(join(round.(w,digits=4)," "),30))  |ω|=$(round(sqrt(sum(x->x*x,w)),4))")
    end

    drift_mid = constant_drift_report(I, traj_mid)
    println("\n  运动常数漂移: T=$(round(drift_mid.T_drift,10)), L²=$(round(drift_mid.L2_drift,10))")

    # ==============================================
    # 实验 2: 绕最大主轴 → 稳定
    # ==============================================
    println("\n" * "-"^40)
    println("实验 2: 绕最大主轴 (I₁) 旋转 — 稳定")

    w0_I1 = [ω_mag, δω, δω]  # 主要绕 x 轴 (I₁)
    T0_I1, _ = motion_constants(I, w0_I1)
    println("  初始 ω = $(round.(w0_I1, 4))")
    println("  初始 T = $(round(T0_I1, 6))")

    traj_I1 = simulate_free_rigid_body(I, w0_I1, 5.0, dt; method=:rk4)
    w_final_I1 = traj_I1[end]
    dot_change = abs(dot(w0_I1, w_final_I1) / (norm(w0_I1) * norm(w_final_I1)))
    println("  5s 后 ω 方向余弦 = $(round(dot_change, 6))")
    println("  → ω 保持接近原始方向 (稳定)")

    # ==============================================
    # 实验 3: 绕最小主轴 → 稳定
    # ==============================================
    println("\n" * "-"^40)
    println("实验 3: 绕最小主轴 (I₃) 旋转 — 稳定")

    w0_I3 = [δω, δω, ω_mag]  # 主要绕 z 轴 (I₃)
    T0_I3, _ = motion_constants(I, w0_I3)
    println("  初始 ω = $(round.(w0_I3, 4))")
    println("  初始 T = $(round(T0_I3, 6))")

    traj_I3 = simulate_free_rigid_body(I, w0_I3, 5.0, dt; method=:rk4)
    w_final_I3 = traj_I3[end]
    dot_change3 = abs(dot(w0_I3, w_final_I3) / (norm(w0_I3) * norm(w_final_I3)))
    println("  5s 后 ω 方向余弦 = $(round(dot_change3, 6))")
    println("  → ω 保持接近原始方向 (稳定)")

    # ==============================================
    # 总结
    # ==============================================
    println("\n" * "="^60)
    println("结论:")
    println("  1. 绕最大惯量轴 (I₁) → 稳定 ✓")
    println("  2. 绕中间惯量轴 (I₂) → 不稳定 ✗ (翻转)")
    println("  3. 绕最小惯量轴 (I₃) → 稳定 ✓")
    println("  → 这就是为什么网球拍绕手柄长轴转是稳定的,")
    println("    绕拍面法向轴转也是稳定的, 但绕中间轴会翻转!")
    println("="^60)
    println("\n✅ tennis_racket_theorem.jl 完成.")
end

main()
