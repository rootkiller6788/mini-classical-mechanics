#!/usr/bin/env julia
# examples/coupled_oscillators_hamiltonian.jl — 耦合振子：相空间中的简正模
include("../src/Hamiltonian.jl")
using .Hamiltonian

function main()
    println("="^60)
    println("  Coupled Oscillators in Phase Space")
    println("="^60)

    m1=m2=1.0; k1=k2=10.0; kc=2.0

    # H = p1²/(2m1) + p2²/(2m2) + ½k1 q1² + ½k2 q2² + ½kc (q1-q2)²
    function H_coupled(q, p)
        return p[1]^2/(2m1) + p[2]^2/(2m2) + 0.5*k1*q[1]^2 + 0.5*k2*q[2]^2 + 0.5*kc*(q[1]-q[2])^2
    end
    
    sys = numerical_hamiltonian_system(H_coupled, 2)

    # ---- 1. 对称模 ----
    println("--- Symmetric Mode (q1=q2) ---")
    q0_sym=[0.3, 0.3]; p0_sym=[0.0, 0.0]
    traj_sym = solve_hamiltonian(sys, q0_sym, p0_sym, 10.0, 0.01; method=:verlet, record_every=20)
    omega_sym = sqrt((k1+kc)/m1)
    println("Expected ω_sym = $(round(omega_sym,digits=4))")
    println("H drift: $(round(maximum(abs.(traj_sym.H_vals .- traj_sym.H_vals[1])),digits=10))")
    # q1 ≈ q2 在全程保持
    max_diff = maximum([abs(q[1]-q[2]) for q in traj_sym.qs])
    println("max|q1-q2| = $(round(max_diff,digits=8)) (should be ~0)")

    # ---- 2. 反对称模 ----
    println("\n--- Antisymmetric Mode (q1=-q2) ---")
    q0_anti=[0.3, -0.3]; p0_anti=[0.0, 0.0]
    traj_anti = solve_hamiltonian(sys, q0_anti, p0_anti, 10.0, 0.01; method=:verlet, record_every=20)
    omega_anti = sqrt((k1+2*kc)/m1)
    println("Expected ω_anti = $(round(omega_anti,digits=4))")
    max_sum = maximum([abs(q[1]+q[2]) for q in traj_anti.qs])
    println("max|q1+q2| = $(round(max_sum,digits=8)) (should be ~0)")

    # ---- 3. 一般初始条件 = 两模叠加 → 拍频 ----
    println("\n--- General IC: Beating ---")
    q0_gen=[0.5, 0.0]; p0_gen=[0.0, 0.0]
    traj_gen = solve_hamiltonian(sys, q0_gen, p0_gen, 20.0, 0.01; method=:verlet, record_every=10)
    n=length(traj_gen.ts)
    println("$n frames recorded, t_max=$(round(traj_gen.ts[end],digits=1))")
    omega_beat = abs(omega_sym - omega_anti)
    T_beat = 2π/omega_beat
    println("Beat frequency: $(round(omega_beat,digits=4)), period: $(round(T_beat,digits=2)) s")

    # 检测拍频：q1 的振幅包络
    q1_vals = [q[1] for q in traj_gen.qs]
    q1_abs = abs.(q1_vals)
    # 找极大值
    peaks = Float64[]
    for i in 2:n-1
        if q1_abs[i] > q1_abs[i-1] && q1_abs[i] > q1_abs[i+1] && q1_abs[i] > 0.1
            push!(peaks, traj_gen.ts[i])
        end
    end
    if length(peaks) >= 2
        beat_period_sim = 2*(peaks[end] - peaks[1])/(length(peaks)-1)
        println("Observed beat period: $(round(beat_period_sim,digits=2)) s")
        println("Theoretical: $(round(T_beat,digits=2)) s")
    end

    # ---- 4. Poincaré 截面 (q1, p1) at q2=0 ----
    println("\n--- Poincare Section at q2=0 ---")
    surface(q,p) = q[2]  # q2 = 0
    pts = poincare_section(traj_gen, surface)
    println("$(length(pts)) crossings at q2=0")
    if length(pts) > 0
        println("  First crossing: q1=$(round(pts[1].q[1],digits=4)), p1=$(round(pts[1].p[1],digits=4))")
    end

    println("\n✅ coupled_oscillators_hamiltonian.jl done.")
end
main()
