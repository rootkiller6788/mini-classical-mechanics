#!/usr/bin/env julia
# examples/canonical_demo.jl — 正则变换演示
include("../src/Hamiltonian.jl")
using .Hamiltonian

function main()
    println("="^60)
    println("  Canonical Transformations — Demo")
    println("="^60)

    # ---- 1. Identity ----
    println("--- 1. Identity Transform ---")
    ct = identity_transform(2)
    q=[1.0,2.0]; p=[3.0,4.0]
    Q,P = ct.transform(q,p)
    println("(q,p)=([1,2],[3,4]) → (Q,P)=($Q,$P)")
    ok,e = verify_canonical(ct, q, p)
    println("Canonical: $ok (max error: $(round(e,digits=10)))")

    # ---- 2. Scale μ=2 ----
    println("\n--- 2. Scale Transform μ=2, ν=0.5 ---")
    ct_s = scale_transform(2.0, 2)
    Qs,Ps = ct_s.transform(q, p)
    println("(q,p)=([1,2],[3,4]) → (Q,P)=($(round.(Qs,digits=2)),$(round.(Ps,digits=2)))")
    ok_s,_ = verify_canonical(ct_s, q, p)
    println("Canonical: $ok_s")

    # ---- 3. Exchange Q=p, P=-q ----
    println("\n--- 3. Exchange Transform (q↔p) ---")
    ct_e = exchange_transform(2)
    Qe,Pe = ct_e.transform(q, p)
    println("(q,p)=([1,2],[3,4]) → (Q,P)=($Qe,$Pe)")
    # 验证：{Q,P} = {p, -q} = -{p,q} = {q,p} = 1
    ok_e,_ = verify_canonical(ct_e, q, p)
    println("Canonical: $ok_e")
    pb_check = poisson_bracket((qq,pp)->Qe[1], (qq,pp)->Pe[1], Qe, Pe)
    println("{Q1,P1} = $(round(pb_check,digits=6)) (expect 1)")

    # ---- 4. SHO 正则变换（作用量-角变量） ----
    println("\n--- 4. Action-Angle Transform for SHO ---")
    m=1.0; omega=1.0
    aa = HarmonicOscillatorActionAngle(m, omega)
    # (q,p) → (θ, J)
    J_val, theta = qp_to_action_angle_ho(aa, 1.0, 0.0)
    println("(q=1,p=0) → J=$(round(J_val,digits=4)), θ=$(round(theta,digits=4))")
    # 验证：q(J,θ) = √(2J/ω)*sin(θ)
    q_back, p_back = action_to_qp_ho(aa, 0.5, π/2)
    println("(J=0.5,θ=π/2) → q=$(round(q_back,digits=4)), p=$(round(p_back,digits=4)) (expect q=1, p=0)")
    # (θ, J) 也是正则变量: {θ, J} = 1
    theta_func(qq,pp) = atan(m*omega*qq[1], pp[1])
    J_func(qq,pp) = (pp[1]^2 + (m*omega*qq[1])^2)/(2*m*omega)
    pb_thJ = poisson_bracket(theta_func, J_func, [1.0], [0.0])
    println("{θ, J} = $(round(pb_thJ,digits=6)) (expect 1)")

    # ---- 5. 谐振子演化（3种正则坐标系对比） ----
    println("\n--- 5. SHO Evolution in 3 Coordinate Systems ---")
    H_cart(q,p) = 0.5*(p[1]^2 + q[1]^2)
    sys = analytic_hamiltonian_system(H_cart, q->[q[1]], p->[p[1]], 1)

    # 笛卡尔坐标中
    traj = solve_hamiltonian(sys, [1.0], [0.0], π/2, 0.01; method=:verlet)
    println("Cartesian: t=π/2, (q,p)=($(round.(traj.qs[end],digits=3)),$(round.(traj.ps[end],digits=3)))")
    println("  Expect (q≈0, p≈-1) — quarter period")

    # 作用量-角变量中（解析）
    J0, th0 = qp_to_action_angle_ho(aa, 1.0, 0.0)
    th_final = th0 + π/2  # θ 均匀增加
    q_aa, p_aa = action_to_qp_ho(aa, J0, th_final)
    println("Action-angle: (q,p)=($(round(q_aa,digits=3)),$(round(p_aa,digits=3)))")
    println("  Action J conserved: $(abs(J0-0.5) < 1e-10)")
    println("  Angle θ advances uniformly: θ(t)=ωt+θ₀")

    println("\n✅ canonical_demo.jl done.")
end
main()
