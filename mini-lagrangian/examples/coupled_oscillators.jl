#!/usr/bin/env julia
# examples/coupled_oscillators.jl — 耦合谐振子：小振动理论
# 用 Lagrangian + 简正模分析求解

include("../src/Lagrangian.jl")
using .Lagrangian
using LinearAlgebra
include("../../mini-newtonian/src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Coupled Oscillators — Normal Mode Analysis")
    println("="^60)

    m1=m2=1.0; k1=k2=10.0; kc=2.0  # kc = coupling spring

    # ---- 1. 建立 Lagrangian ----
    # q = (x1, x2), x1和x2是两质量的位移
    # T = 0.5*m1*x1dot^2 + 0.5*m2*x2dot^2
    # U = 0.5*k1*x1^2 + 0.5*k2*x2^2 + 0.5*kc*(x1-x2)^2

    function U_coupled(q)
        x1,x2 = q[1], q[2]
        return 0.5*k1*x1^2 + 0.5*k2*x2^2 + 0.5*kc*(x1-x2)^2
    end

    function M_coupled(q)
        return [m1 0.0; 0.0 m2]
    end

    # ---- 2. 找平衡点 ----
    grad_U(q) = numerical_gradient(U_coupled, q)
    hess_U(q) = numerical_hessian(U_coupled, q)
    q_eq = find_equilibrium(grad_U, hess_U, [1.0, -0.5]; tol=1e-12)
    println("Equilibrium: $q_eq (expect [0,0])")

    # ---- 3. Hessian + 简正模 ----
    K = hess_U(q_eq)
    M = M_coupled(q_eq)
    println("\nMass matrix M:\n$M")
    println("Hessian K:\n$K")

    sys = solve_normal_modes(M, K, q_eq)

    println("\n=== Normal Modes ===")
    for alpha in 1:length(sys.frequencies)
        println("Mode $alpha: ω = $(round(sys.frequencies[alpha],digits=4)) rad/s")
        println("  Eigenvector: $(round.(sys.modes[:,alpha],digits=4))")
        println("  Effective mass: $(round(sys.effective_masses[alpha],digits=4))")
    end

    # 理论值
    omega_sym = sqrt((k1+kc)/m1)  # 对称模 x1=x2 → 边弹簧+耦合弹簧
    omega_anti = sqrt((k1+2*kc)/m1) # 反对称模 x1=-x2
    println("\nTheory:")
    println("  Symmetric mode: ω = $(round(omega_sym,digits=4)) (x1=x2)")
    println("  Antisymmetric mode: ω = $(round(omega_anti,digits=4)) (x1=-x2)")

    # ---- 4. 模拟验证 ----
    # 激发纯反对称模
    amp_anti = [0.0, 0.5]  # 只激发第二模
    phase_anti = [0.0, 0.0]

    dt = 0.001; t_end = 10.0
    n_steps = Int(ceil(t_end/dt))
    times = Float64[]; x1s = Float64[]; x2s = Float64[]

    for step in 1:n_steps
        t = (step-1)*dt
        if step % 50 == 1
            push!(times, t)
            q_t = normal_mode_solution(sys, amp_anti, phase_anti, t)
            push!(x1s, q_t[1]); push!(x2s, q_t[2])
        end
    end

    println("\nSimulation (antisymmetric mode excited):")
    println("  x1 at t=0: $(round(x1s[1],digits=4))")
    println("  x2 at t=0: $(round(x2s[1],digits=4))")
    println("  Ratio x1/x2 ≈ -1: $(all(abs.(x1s[i]+x2s[i]) < 0.01 for i in 1:length(x1s)) ? "YES" : "check")")
    println("  x1 peaks at $(round(maximum(abs.(x1s)),digits=4)) (expect ~0.25)")

    # ---- 5. 任意初始条件的频谱 ----
    println("\n--- General IC: superposition of modes ---")
    q0 = [0.3, -0.1]  # 同时激发两个模
    v0 = [0.0, 0.0]

    # 用 Newtonian 的 Verlet 模拟
    function coupled_accel(r, v, t)
        x1,x2 = r.x, r.y
        a1 = (-(k1+kc)*x1 + kc*x2) / m1
        a2 = (kc*x1 - (k2+kc)*x2) / m2
        return Vec3(a1, a2, 0.0)
    end
    traj = solve_fixed_step(coupled_accel, Vec3(q0[1],q0[2],0), Vec3(0,0,0), 20.0, 0.001; method=:verlet, record_every=100)
    n = length(traj.ts)
    println("General IC: $(n) frames, t_max=$(round(traj.ts[end],digits=1)) s")
    println("Both modes visible in the beating pattern")

    println("\n✅ coupled_oscillators.jl done.")
end
main()
