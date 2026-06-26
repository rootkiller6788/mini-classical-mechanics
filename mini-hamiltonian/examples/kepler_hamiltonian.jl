#!/usr/bin/env julia
# examples/kepler_hamiltonian.jl — Kepler问题：Runge-Lenz矢量与轨道参数
include("../src/Hamiltonian.jl")
using .Hamiltonian

function main()
    println("="^60)
    println("  Kepler Problem — Hamiltonian + Runge-Lenz Vector")
    println("="^60)

    m=1.0; k=1.0  # GMm
    kp = KeplerProblem(m, k)

    # 初始条件（椭圆轨道 e=0.5）
    a = 1.0; e = 0.5
    r_peri = a*(1-e)
    v_peri = sqrt(k/m * (1+e)/(a*(1-e)))  # 近地点速度

    r0 = [r_peri, 0.0, 0.0]
    p0 = [0.0, m*v_peri, 0.0]

    # 从守恒量推导轨道参数
    params = kepler_orbit_params(kp, r0, p0)
    println("Initial conditions: r=($(r0[1]),0,0), p=(0,$(round(p0[2],digits=3)),0)")
    println("\nDerived from integrals of motion:")
    println("  E = $(round(params.E,digits=4))")
    println("  L² = $(round(params.L2,digits=4))")
    println("  e = $(round(params.eccentricity,digits=4)) (target: $e)")
    println("  a = $(round(params.semi_major,digits=4)) (target: $a)")
    println("  A = [$(round(params.A[1],digits=4)), $(round(params.A[2],digits=4)), $(round(params.A[3],digits=4))]")

    # 验证积分对合
    verify_kepler_integrals(kp, r0, p0)

    # Hamilton 演化
    H_func(q,p) = kepler_energy(kp, q, p)
    sys = numerical_hamiltonian_system(H_func, 3)
    dt = 0.005; t_end = 2π  # 一个周期
    traj = solve_hamiltonian(sys, r0, p0, t_end, dt; method=:verlet, record_every=50)

    # 验证 Runge-Lenz 矢量在演化中守恒
    n = length(traj.ts)
    A_mags = Float64[]
    for i in 1:n
        A = runge_lenz_vector(kp, traj.qs[i], traj.ps[i])
        push!(A_mags, sqrt(sum(A.^2)))
    end
    A_drift = (maximum(A_mags) - minimum(A_mags)) / A_mags[1]
    println("\nRunge-Lenz magnitude drift: $(round(A_drift,digits=10))")
    println("H drift: $(round(maximum(abs.(traj.H_vals .- traj.H_vals[1])),digits=10))")

    # 轨道几何验证
    r_vals = [sqrt(sum(q.^2)) for q in traj.qs]
    r_min, r_max = minimum(r_vals), maximum(r_vals)
    e_numerical = (r_max - r_min) / (r_max + r_min)
    a_numerical = (r_max + r_min) / 2
    println("\nNumerical orbit:")
    println("  r_peri=$(round(r_min,digits=4)), r_apo=$(round(r_max,digits=4))")
    println("  e_numerical=$(round(e_numerical,digits=4)) (target: $e)")
    println("  a_numerical=$(round(a_numerical,digits=4)) (target: $a)")

    println("\n✅ kepler_hamiltonian.jl done.")
end
main()
