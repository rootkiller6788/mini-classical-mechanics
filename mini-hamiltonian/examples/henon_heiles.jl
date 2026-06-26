#!/usr/bin/env julia
# examples/henon_heiles.jl — Henon-Heiles系统：哈密顿混沌
include("../src/Hamiltonian.jl")
using .Hamiltonian

function main()
    println("="^60)
    println("  Henon-Heiles System — Hamiltonian Chaos")
    println("="^60)

    # H = 0.5*(px²+py²) + 0.5*(x²+y²) + x²*y - y³/3
    # 天体力学中星系势的简化模型
    # E < 1/6 → 规则运动; E > 1/6 → 混沌出现

    function HH(q, p)
        x,y = q[1], q[2]; px,py = p[1], p[2]
        return 0.5*(px^2+py^2) + 0.5*(x^2+y^2) + x^2*y - y^3/3
    end
    sys = numerical_hamiltonian_system(HH, 2)
    dt = 0.01

    println("E_crit = 1/6 ≈ 0.1667 — above this, chaos emerges\n")

    # 不同能量下的 Poincare 截面
    for E in [0.05, 0.10, 0.125, 0.15, 0.1667]
        # 从 y=0, py>0 开始，x 由 E 决定
        # H(x,0,0,py) = 0.5*py² + 0.5*x² = E
        # 取 x=0 → py=√(2E)
        q0=[0.0, 0.0]; p0=[0.0, sqrt(2E)]
        t_end = 500.0
        traj = solve_hamiltonian(sys, q0, p0, t_end, dt; method=:verlet, record_every=5)

        # Poincare 截面: y=0, py>0
        surface(q,p) = q[2]
        pts = poincare_section(traj, surface)
        # 只保留 py>0 的点
        valid_pts = [pt for pt in pts if pt.p[2] > 0]
        n_pts = length(valid_pts)

        # 计算截面点的散布程度（规则→闭合曲线，混沌→弥散）
        if n_pts > 2
            x_vals = [pt.q[1] for pt in valid_pts]
            spread = maximum(x_vals) - minimum(x_vals)
            println("E=$E: $(n_pts) Poincare points, x_spread=$(round(spread,digits=4))")
            if E < 0.12
                println("  → Regular (points lie on invariant curves)")
            elseif E < 0.16
                println("  → Mixed (some regular islands, some chaotic sea)")
            else
                println("  → Mostly chaotic (points fill area)")
            end
        end
    end

    # 最大 Lyapunov 指数近似（两条邻近轨道的分离）
    println("\n--- Lyapunov Exponent Estimate ---")
    d0 = 1e-8
    for E in [0.0833, 0.125, 0.1667]
        py0 = sqrt(2E)
        qA=[0.0, 0.0]; pA=[0.0, py0]
        qB=[d0, 0.0]; pB=[0.0, py0]

        t_end = 100.0; n_steps = Int(ceil(t_end/dt))
        d = d0; log_sum = 0.0; count = 0
        for _ in 1:n_steps
            qA,pA = stormer_verlet_step(sys, qA, pA, dt)
            qB,pB = stormer_verlet_step(sys, qB, pB, dt)
            d_new = sqrt(sum((qA-qB).^2) + sum((pA-pB).^2))
            if d_new > 1e-15 && d > 1e-15
                log_sum += log(d_new / d)
                count += 1
                # 重归一化
                factor = d0 / d_new
                qB = qA + factor*(qB - qA)
                pB = pA + factor*(pB - pA)
            end
            d = d0
        end
        lyap = count > 0 ? log_sum / (count*dt) : 0.0
        println("  E=$E: λ ≈ $(round(lyap,digits=6)) ($(lyap > 0.005 ? "CHAOTIC" : "regular"))")
    end

    println("\n✅ henon_heiles.jl done.")
end
main()
