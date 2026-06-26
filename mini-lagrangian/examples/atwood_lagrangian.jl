#!/usr/bin/env julia
# examples/atwood_lagrangian.jl — Atwood机：拉格朗日 vs 牛顿对比
# 展示约束系统的拉格朗日处理 vs 牛顿力学中的张力分析

include("../src/Lagrangian.jl")
using .Lagrangian
include("../../mini-newtonian/src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Atwood Machine — Lagrangian vs Newtonian")
    println("="^60)

    m1,m2 = 5.0, 3.0; g = G_EARTH

    # ---- 1. Lagrangian formulation ----
    # 约束: x1 + x2 = const (绳子不可伸长)
    # 广义坐标 q = x (m1向下的位移)
    # x1 = x, x2 = -x (滑轮另一侧向上)
    # T = 0.5*m1*xdot^2 + 0.5*m2*(-xdot)^2 = 0.5*(m1+m2)*xdot^2
    # U = -m1*g*x - m2*g*(-x) = (m2-m1)*g*x (设滑轮处U=0)
    # L = 0.5*(m1+m2)*xdot^2 + (m1-m2)*g*x
    # EL: (m1+m2)*xddot = (m1-m2)*g → a = (m1-m2)/(m1+m2)*g

    function atwood_L(q, qdot)
        x, xdot = q[1], qdot[1]
        return 0.5*(m1+m2)*xdot^2 + (m1-m2)*g*x
    end
    atwood_grad_q(q,qdot) = [(m1-m2)*g]
    atwood_grad_qdot(q,qdot) = [(m1+m2)*qdot[1]]
    atwood_M(q) = [m1+m2]
    sys = EulerLagrangeSystem(1, atwood_L, atwood_grad_q, atwood_grad_qdot, atwood_M)
    f_ode = el_to_first_order(sys)

    y0 = [0.0, 0.0]; dt=0.001; t_end=2.0
    n_steps = Int(ceil(t_end/dt)); y=copy(y0)
    for _ in 1:n_steps; y = rk4_step(f_ode, 0.0, y, dt); end
    a_lag = y[2]/t_end
    println("Lagrangian: a = $(round(a_lag,digits=4)) m/s²")

    # ---- 2. Newtonian formula ----
    a_newt, _, T = atwood_machine(m1, m2)
    println("Newtonian: a = $(round(a_newt,digits=4)) m/s², T = $(round(T,digits=2)) N")
    println("Match: $(abs(a_lag - a_newt) < 1e-4 ? "YES" : "NO")")

    # ---- 3. Parameter sweep ----
    println("\n--- Parameter Sweep ---")
    for (a,b) in [(2,1),(4,1),(10,1),(10,9),(100,99)]
        a_val, _, T_val = atwood_machine(Float64(a), Float64(b))
        println("  m1=$a, m2=$b → a=$(round(a_val,digits=3)) m/s², T=$(round(T_val,digits=1)) N")
    end
    println("  (m1≈m2 → a≈0 → 可精确测g的实验设计)")

    # ---- 4. Energy conservation ----
    E0 = energy_from_lagrangian(atwood_L, [0.0], [0.0], atwood_grad_qdot)
    E_f = energy_from_lagrangian(atwood_L, [y[1]], [y[2]], atwood_grad_qdot)
    println("\nEnergy: E0=$E0, E_final=$E_f, drift=$(round((E_f-E0)/max(abs(E0),1e-300),digits=10))")

    println("\n✅ atwood_lagrangian.jl done.")
end
main()
