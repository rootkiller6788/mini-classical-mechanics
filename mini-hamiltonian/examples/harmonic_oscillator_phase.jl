#!/usr/bin/env julia
# examples/harmonic_oscillator_phase.jl — 谐振子相空间
include("../src/Hamiltonian.jl")
using .Hamiltonian

function main()
    println("="^60)
    println("  Harmonic Oscillator — Phase Space Portrait")
    println("="^60)

    m=1.0; k=1.0; omega=sqrt(k/m)
    H(q,p) = p[1]^2/(2m) + 0.5*k*q[1]^2
    grad_q(q,p) = [k*q[1]]
    grad_p(q,p) = [p[1]/m]
    sys = analytic_hamiltonian_system(H, grad_q, grad_p, 1)

    # 多条能量等高线
    println("Energy contours (E = 0.5, 1, 2, 4, 8):")
    for E in [0.5, 1.0, 2.0, 4.0, 8.0]
        amp = sqrt(2E/k)
        pmax = sqrt(2m*E)
        println("  E=$E: q_amp=±$(round(amp,digits=3)), p_max=±$(round(pmax,digits=3))")
    end

    # 对比 Hamilton vs analytic
    q0=[1.0]; p0=[0.0]; dt=0.01; t_end=2π
    traj = solve_hamiltonian(sys, q0, p0, t_end, dt; method=:verlet, record_every=10)
    n = length(traj.ts)
    println("\nVerlet: period ~$(round(t_end,digits=3)), n=$(n)")
    println("  q_final=$(round(traj.qs[end][1],digits=6)) (expect 1.0)")
    println("  H_drift=$(round(maximum(abs.(traj.H_vals .- traj.H_vals[1])),digits=10))")

    # 作用量-角变量
    aa = HarmonicOscillatorActionAngle(m, omega)
    J, theta = qp_to_action_angle_ho(aa, 1.0, 0.0)
    println("\nAction-angle: J=$(round(J,digits=4)), θ=$(round(theta,digits=4))")
    println("E from J: $(round(energy_from_action_ho(aa,J),digits=4))")

    # 相空间椭圆验证
    println("\nPhase ellipse check: p² + ω²q² = const?")
    for i in [1, n÷4, n÷2, 3n÷4, n]
        q=traj.qs[i][1]; p=traj.ps[i][1]
        invar = p^2 + omega^2*q^2
        println("  t=$(round(traj.ts[i],digits=2)): $(round(invar,digits=6))")
    end

    println("\n✅ harmonic_oscillator_phase.jl done.")
end
main()
