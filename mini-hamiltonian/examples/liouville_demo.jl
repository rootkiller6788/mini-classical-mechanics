#!/usr/bin/env julia
# examples/liouville_demo.jl — Liouville 定理数值验证
include("../src/Hamiltonian.jl")
using .Hamiltonian

function main()
    println("="^60)
    println("  Liouville Theorem Verification")
    println("="^60)

    # SHO: H = p²/2 + q²/2
    H(q,p) = 0.5*(p[1]^2 + q[1]^2)
    grad_q(q,p) = [q[1]]
    grad_p(q,p) = [p[1]]
    sys = analytic_hamiltonian_system(H, grad_q, grad_p, 1)

    println("System: SHO, H = 0.5(p² + q²)")
    println("Hamiltonian flow should preserve phase space area.\n")

    # 测试不同积分器对 Liouville 定理的保持
    for method in [:verlet, :symplectic_euler]
        println("--- $method ---")
        ensemble = generate_ensemble([1.0], [0.0], 50; sigma_q=0.05, sigma_p=0.05)
        V0 = ensemble_volume_proxy(ensemble)
        dt = 0.05; t_end = 2π  # 一个周期

        # 推进
        for _ in 1:Int(ceil(t_end/dt))
            for k in 1:50
                q,p = ensemble[k]
                if method == :verlet
                    q_new, p_new = stormer_verlet_step(sys, q, p, dt)
                else
                    q_new, p_new = symplectic_euler_step(sys, q, p, dt)
                end
                ensemble[k] = (q_new, p_new)
            end
        end

        V1 = ensemble_volume_proxy(ensemble)
        drift = (V1 - V0) / (abs(V0) + 1e-300)
        println("  Initial volume proxy: $(round(V0,digits=6))")
        println("  After 1 period:       $(round(V1,digits=6))")
        println("  Relative drift:       $(round(drift,digits=8))")
        println("  Liouville preserved:  $(abs(drift) < 0.05 ? "YES" : "WARNING")")
    end

    # 非 Hamilton 流的对照：显式 Euler（非辛）
    println("\n--- Explicit Euler (non-symplectic) ---")
    ensemble = generate_ensemble([1.0], [0.0], 50; sigma_q=0.05, sigma_p=0.05)
    V0 = ensemble_volume_proxy(ensemble)
    dt = 0.01; t_end = 2π
    for _ in 1:Int(ceil(t_end/dt))
        for k in 1:50
            q,p = ensemble[k]
            dH_dq, dH_dp = sys.grad_H(q, p)
            q_new = q + dt*dH_dp
            p_new = p - dt*dH_dq
            ensemble[k] = (q_new, p_new)
        end
    end
    V1 = ensemble_volume_proxy(ensemble)
    drift = (V1 - V0) / (abs(V0) + 1e-300)
    println("  Initial: $(round(V0,digits=6)), After 1 period: $(round(V1,digits=6))")
    println("  Drift: $(round(drift,digits=4)) — phase space area GROWS (dissipative!)")
    println("  → Non-symplectic integrators violate Liouville's theorem")

    println("\n✅ liouville_demo.jl done.")
end
main()
