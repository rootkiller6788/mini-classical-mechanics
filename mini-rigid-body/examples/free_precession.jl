#!/usr/bin/env julia
# free_precession.jl — 无力矩刚体自由进动
# 参考: Goldstein Ch.5.5, Landau Ch.37

include("../src/RigidBody.jl"); using .RigidBody

function main()
    println("="^60)
    println("  Rigid Body — Free Precession (Poinsot)")
    println("="^60)

    # 构造非对称刚体: 4 个点质量
    m = [1.0, 1.0, 1.0, 1.0]
    r = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
         [0.0, 1.0, 0.0], [0.0, 0.0, 0.5]]

    I = inertia_tensor(m, r)
    println("\n惯性张量 (质心系):")
    display(round.(inertia_matrix(I), digits=4))

    pm = principal_moments(I)
    pa = principal_axes(I)
    println("\n主转动惯量: I₁=$(round(pm[1],4)), I₂=$(round(pm[2],4)), I₃=$(round(pm[3],4))")

    # 稳定性分析
    stab = stability_analysis(I)
    println("\n主轴稳定性 (网球拍定理):")
    for (k, v) in stab
        println("  $k: $v")
    end

    # ==============================================
    # 模拟: 初始角速度靠近中间轴 → 不稳定翻转
    # ==============================================
    w0 = [0.1, 1.0, 0.05]  # 主要绕 I₂ (中间轴)
    T0 = rotational_kinetic_energy(I, w0)
    L0 = angular_momentum_rigid(I, w0)
    println("\n初始条件:")
    println("  ω₀ = $(round.(w0, digits=4))")
    println("  T  = $(round(T0, 6))")
    println("  |L| = $(round(norm(L0), 6))")

    # RK4 积分
    dt = 0.01
    t_end = 10.0
    traj = simulate_free_rigid_body(I, w0, t_end, dt; method=:rk4)
    n = length(traj)
    println("\n模拟: $(n-1) 步, Δt=$dt, t_end=$t_end (RK4)")

    # 运动常数漂移
    drift = constant_drift_report(I, traj)
    println("\n运动常数相对漂移:")
    println("  T  drift: $(round(drift.T_drift, 12))")
    println("  L² drift: $(round(drift.L2_drift, 12))")

    # 输出几个时刻的角速度
    println("\n角速度演化 (体坐标系):")
    for frac in [0.0, 0.25, 0.5, 0.75, 1.0]
        idx = max(1, Int(round(frac * (n - 1))) + 1)
        t_val = (idx - 1) * dt
        w = traj[idx]
        println("  t=$(round(t_val,2)): ω = $(round.(w, digits=5))")
    end

    # ==============================================
    # 对比: 积分器精度对比
    # ==============================================
    println("\n--- 积分器对比 (运动常数保真度) ---")
    for method in [:euler, :rk4, :verlet]
        traj_m = simulate_free_rigid_body(I, w0, 5.0, 0.01; method=method)
        drift_m = constant_drift_report(I, traj_m)
        println("  $(rpad(method, 8)): T_drift=$(round(drift_m.T_drift,10)), L²_drift=$(round(drift_m.L2_drift,10))")
    end

    # ==============================================
    # Polhode 曲线数据 (可导出可视化)
    # ==============================================
    println("\n--- Polhode 曲线采样 ---")
    w_samples = polhode_curve(I, w0, 5.0, 0.05; method=:rk4)
    ω1_vals = [w[1] for w in w_samples]
    ω2_vals = [w[2] for w in w_samples]
    ω3_vals = [w[3] for w in w_samples]
    println("  样本数: $(length(w_samples))")
    println("  ω₁ ∈ [$(round(minimum(ω1_vals),4)), $(round(maximum(ω1_vals),4))]")
    println("  ω₂ ∈ [$(round(minimum(ω2_vals),4)), $(round(maximum(ω2_vals),4))]")
    println("  ω₃ ∈ [$(round(minimum(ω3_vals),4)), $(round(maximum(ω3_vals),4))]")

    # 惯性椭球半轴
    a, b, c = inertia_ellipsoid(I)
    println("\n惯性椭球半轴: a=$(round(a,4)), b=$(round(b,4)), c=$(round(c,4))")

    println("\n✅ free_precession.jl 完成.")
end

main()
