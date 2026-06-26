#!/usr/bin/env julia
# examples/bead_rotating_hoop.jl — 旋转圆环上的珠子
# 经典拉格朗日力学例题：约束系统 + 有效势 + 分岔

include("../src/Lagrangian.jl")
using .Lagrangian
include("../../mini-newtonian/src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Bead on a Rotating Hoop")
    println("="^60)

    R = 1.0; m = 1.0; g = G_EARTH

    # 圆环绕竖直轴以角速度 ω 旋转
    # 广义坐标: θ（珠子与最低点的夹角）
    # 动能来自两个分量：沿环运动 + 随环旋转
    # v_θ = R*θdot（切向）
    # v_φ = R*sin(θ)*ω（环的旋转带给珠子的速度）
    # T = 0.5*m*(R²θdot² + R²sin²(θ)ω²)
    # U = m*g*R*(1-cos(θ))
    # L = T - U
    # ∂L/∂θ = m*R²*sin(θ)*cos(θ)*ω² - m*g*R*sin(θ)
    # ∂L/∂θdot = m*R²*θdot
    # EL: θddot = sin(θ)*cos(θ)*ω² - (g/R)*sin(θ) = sin(θ)*(ω²*cos(θ) - g/R)

    println("R = $R m, m = $m kg, g = $g m/s²")
    println()

    # 平衡点分析：θddot = 0 → sin(θ)=0 或 ω²*cos(θ) = g/R
    println("--- Equilibrium Points ---")
    # θ=0 和 θ=π 始终是平衡点
    println("  θ=0 (bottom): always an equilibrium")
    println("  θ=π (top):    always an equilibrium")

    # 当 ω > ω_crit 时，出现新的平衡点
    omega_crit = sqrt(g / R)
    println("  Critical ω = sqrt(g/R) = $(round(omega_crit, digits=3)) rad/s")

    for omega in [0.0, 0.5*omega_crit, omega_crit, 1.5*omega_crit, 3.0*omega_crit]
        if omega > omega_crit
            theta_eq = acos(g / (R * omega^2))
            println("  ω=$(round(omega,digits=2)): new equilibrium at θ=±$(round(rad2deg(theta_eq),digits=1))°")
        end
    end

    # ---- 模拟：从底部微扰，观察不同 ω 下的运动 ----
    println("\n--- Simulation ---")

    for (label, omega) in [("Subcritical (ω=0.7ωc)", 0.7*omega_crit),
                           ("Critical (ω=ωc)", omega_crit),
                           ("Supercritical (ω=2ωc)", 2.0*omega_crit)]

        function bead_accel(r, v, t)
            theta = r.x
            theta_dot = v.x
            # 有效加速度（含离心效应）
            alpha = sin(theta) * (omega^2 * cos(theta) - g/R)
            return Vec3(alpha, 0, 0)
        end

        theta0 = deg2rad(5)  # 从接近底部开始
        r0 = Vec3(theta0, 0, 0)
        v0 = Vec3(0, 0, 0)
        dt = 0.001; t_end = 10.0

        traj = solve_fixed_step(bead_accel, r0, v0, t_end, dt; method=:verlet, record_every=50)

        # 计算平均位置（最后 20% 的时间）
        n = length(traj.ts)
        n_avg = max(1, n ÷ 5)
        avg_theta = mean([traj.positions[i].x for i in n-n_avg+1:n])
        println("  $label:")
        println("    Final θ ≈ $(round(rad2deg(avg_theta), digits=2))° (avg last 20%)")

        if omega > omega_crit
            theta_eq_theory = acos(g / (R * omega^2))
            println("    Theory equilibrium: ±$(round(rad2deg(theta_eq_theory), digits=1))°")
        end
    end

    # ---- 有效势分析 ----
    println("\n--- Effective Potential ---")
    # U_eff(θ) = U(θ) - 0.5*m*R²*sin²(θ)*ω²  (离心势)
    #           = m*g*R*(1-cos(θ)) - 0.5*m*R²*sin²(θ)*ω²

    for (label, omega) in [("ω=0.5ωc", 0.5*omega_crit), ("ω=ωc", omega_crit), ("ω=2ωc", 2.0*omega_crit)]
        println("  $label:")
        for theta_deg in [0, 30, 60, 90, 120, 150, 180]
            theta = deg2rad(theta_deg)
            U_eff = m*g*R*(1-cos(theta)) - 0.5*m*R^2*sin(theta)^2*omega^2
            println("    θ=$(theta_deg)°: U_eff = $(round(U_eff, digits=3))")
        end
    end

    println("\n✅ bead_rotating_hoop.jl done.")
end

main()
