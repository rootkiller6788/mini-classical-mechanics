#!/usr/bin/env julia
# examples/pendulum_phase.jl — 单摆相图：分离线、振动区与旋转区
include("../src/Hamiltonian.jl")
using .Hamiltonian

function main()
    println("="^60)
    println("  Pendulum Phase Portrait — Libration & Rotation")
    println("="^60)

    m=1.0; L=1.0; g=9.81
    H(q,p) = p[1]^2/(2*m*L^2) - m*g*L*cos(q[1])
    grad_q(q,p) = [m*g*L*sin(q[1])]
    grad_p(q,p) = [p[1]/(m*L^2)]
    sys = analytic_hamiltonian_system(H, grad_q, grad_p, 1)
    dt=0.01
    sep_energy = m*g*L  # 分离线能量（恰好到达顶点 θ=π）

    println("Separatrix energy: $(round(sep_energy,digits=3)) J")

    # ---- 1. 振动区 (E < sep_energy) ----
    println("\n--- Libration (E < separatrix) ---")
    for E_frac in [0.1, 0.5, 0.9]
        E = E_frac * sep_energy
        theta_max = acos(-E/(m*g*L))
        q0=[theta_max]; p0=[0.0]
        traj = solve_hamiltonian(sys, q0, p0, 4*2π*sqrt(L/g), dt; method=:verlet, record_every=20)
        println("  E=$(round(E,digits=2)): θ_max=±$(round(rad2deg(theta_max),digits=1))°, n=$(length(traj.ts)) frames")
    end

    # ---- 2. 分离线 (E = sep_energy) ----
    println("\n--- Separatrix (E = separatrix) ---")
    q0=[deg2rad(179.0)]; p0=[0.0]  # 近乎恰好到顶点
    traj_sep = solve_hamiltonian(sys, q0, p0, 10.0, dt; method=:verlet, record_every=20)
    println("  E≈sep: $(length(traj_sep.ts)) frames, final θ=$(round(rad2deg(traj_sep.qs[end][1]),digits=1))°")

    # ---- 3. 旋转区 (E > sep_energy) ----
    println("\n--- Rotation (E > separatrix) ---")
    for E_mul in [1.5, 3.0]
        E = E_mul * sep_energy
        # 从 θ=0 出发，动量足够大做完整旋转
        p0_val = sqrt(2*m*L^2*(E + m*g*L))
        q0=[0.0]; p0=[p0_val]
        traj = solve_hamiltonian(sys, q0, p0, 10.0, dt; method=:verlet, record_every=20)
        theta_range = maximum(traj.qs)[1] - minimum(traj.qs)[1]
        println("  E=$(round(E,digits=1)): θ_range=$(round(rad2deg(theta_range),digits=1))°, n=$(length(traj.ts))")
    end

    # ---- 4. Poincare 截面 (θ mod 2π, p) ----
    println("\n--- Poincare Section (stroboscopic at θ=0) ---")
    surface(q,p) = q[1] - 0.0  # θ=0 截面
    points = poincare_section(traj_sep, surface)
    println("  $(length(points)) crossings at θ=0 for near-separatrix orbit")

    println("\n✅ pendulum_phase.jl done.")
end
main()
