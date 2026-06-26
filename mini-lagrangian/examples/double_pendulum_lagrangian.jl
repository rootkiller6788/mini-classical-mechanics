#!/usr/bin/env julia
# examples/double_pendulum_lagrangian.jl — 双摆：拉格朗日方法的威力
# 牛顿力学做双摆极其繁琐，拉格朗日方法几步到位

include("../src/Lagrangian.jl")
using .Lagrangian
include("../../mini-newtonian/src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Double Pendulum via Lagrangian")
    println("="^60)

    L1=1.0; L2=1.0; m1=1.0; m2=1.0; g=G_EARTH

    # ---- 1. 推导拉格朗日量 ----
    # q = (θ1, θ2)
    # x1 = L1*sin(θ1),              z1 = -L1*cos(θ1)
    # x2 = L1*sin(θ1)+L2*sin(θ2),   z2 = -L1*cos(θ1)-L2*cos(θ2)
    # T = 0.5*m1*(x1dot²+z1dot²) + 0.5*m2*(x2dot²+z2dot²)
    # 展开后:
    # T = 0.5*(m1+m2)*L1²*θ1dot² + 0.5*m2*L2²*θ2dot² + m2*L1*L2*θ1dot*θ2dot*cos(θ1-θ2)
    # U = -(m1+m2)*g*L1*cos(θ1) - m2*g*L2*cos(θ2)
    # L = T - U

    function double_pendulum_L(q, qdot)
        th1, th2 = q[1], q[2]
        dth1, dth2 = qdot[1], qdot[2]
        cos_d = cos(th1 - th2)
        T = 0.5*(m1+m2)*L1^2*dth1^2 + 0.5*m2*L2^2*dth2^2 + m2*L1*L2*dth1*dth2*cos_d
        U = -(m1+m2)*g*L1*cos(th1) - m2*g*L2*cos(th2)
        return T - U
    end

    function double_pendulum_grad_L_q(q, qdot)
        th1, th2 = q[1], q[2]
        dth1, dth2 = qdot[1], qdot[2]
        cos_d = cos(th1 - th2); sin_d = sin(th1 - th2)
        # ∂L/∂θ1
        dL_dth1 = -m2*L1*L2*dth1*dth2*sin_d - (m1+m2)*g*L1*sin(th1)
        # ∂L/∂θ2
        dL_dth2 =  m2*L1*L2*dth1*dth2*sin_d - m2*g*L2*sin(th2)
        return [dL_dth1, dL_dth2]
    end

    function double_pendulum_grad_L_qdot(q, qdot)
        th1, th2 = q[1], q[2]
        dth1, dth2 = qdot[1], qdot[2]
        cos_d = cos(th1 - th2)
        p1 = (m1+m2)*L1^2*dth1 + m2*L1*L2*dth2*cos_d
        p2 = m2*L2^2*dth2 + m2*L1*L2*dth1*cos_d
        return [p1, p2]
    end

    function double_pendulum_M(q)
        th1, th2 = q[1], q[2]
        cos_d = cos(th1 - th2)
        M11 = (m1+m2)*L1^2
        M12 = m2*L1*L2*cos_d
        M22 = m2*L2^2
        return [M11 M12; M12 M22]
    end

    sys = EulerLagrangeSystem(2, double_pendulum_L, double_pendulum_grad_L_q, double_pendulum_grad_L_qdot, double_pendulum_M)
    f_ode = el_to_first_order(sys)

    # ---- 2. 与牛顿法对比 ----
    theta1_0 = deg2rad(90)
    theta2_0 = deg2rad(90)
    y0 = [theta1_0, theta2_0, 0.0, 0.0]
    dt = 0.0005; t_end = 5.0

    n_steps = Int(ceil(t_end/dt))
    y = copy(y0)
    record_every = 50
    traj_theta1 = Float64[]; traj_theta2 = Float64[]
    times = Float64[]

    for step in 1:n_steps
        if step % record_every == 1
            push!(times, (step-1)*dt)
            push!(traj_theta1, y[1]); push!(traj_theta2, y[2])
        end
        y = rk4_step(f_ode, (step-1)*dt, y, dt)
    end

    # 牛顿法双摆（mini-newtonian 的版本）
    function double_pend_newton_accel(r, v, t)
        th1, th2 = r.x, r.y
        dth1, dth2 = v.x, v.y
        delta = th2 - th1
        sd, cd = sin(delta), cos(delta)
        denom = (m1+m2)*L1 - m2*L1*cd^2
        a1 = (m2*g*sin(th2)*cd - m2*L2*dth2^2*sd - (m1+m2)*g*sin(th1) - m2*L1*dth1^2*sd*cd) / denom
        a2 = ((m1+m2)*(g*sin(th1)*cd + L1*dth1^2*sd - g*sin(th2)) + m2*L2*dth2^2*sd*cd) / (L2/L1*denom)
        return Vec3(a1, a2, 0)
    end

    r0_dp = Vec3(theta1_0, theta2_0, 0); v0_dp = Vec3(0, 0, 0)
    traj_newton = solve_fixed_step(double_pend_newton_accel, r0_dp, v0_dp, t_end, dt; method=:rk4, record_every=record_every)

    # 对比
    println("Lagrangian vs Newtonian (at t=$t_end s):")
    n_newton = length(traj_newton.ts)
    n_lag = length(times)
    n_compare = min(n_lag, n_newton)
    println("  Lagrangian θ1=$(round(rad2deg(traj_theta1[end]),digits=4))°, θ2=$(round(rad2deg(traj_theta2[end]),digits=4))°")
    println("  Newtonian  θ1=$(round(rad2deg(traj_newton.positions[end].x),digits=4))°, θ2=$(round(rad2deg(traj_newton.positions[end].y),digits=4))°")

    # 能量守恒
    E0 = energy_from_lagrangian(double_pendulum_L, [theta1_0, theta2_0], [0.0, 0.0], double_pendulum_grad_L_qdot)
    E_final = energy_from_lagrangian(double_pendulum_L, [y[1], y[2]], [y[3], y[4]], double_pendulum_grad_L_qdot)
    println("\nEnergy: E0=$(round(E0,digits=4)), E_final=$(round(E_final,digits=4)), drift=$(round((E_final-E0)/abs(E0),digits=8))")

    println("\n✅ double_pendulum_lagrangian.jl done.")
end

main()
