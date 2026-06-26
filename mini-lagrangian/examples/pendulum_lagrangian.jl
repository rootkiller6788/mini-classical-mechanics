#!/usr/bin/env julia
# examples/pendulum_lagrangian.jl — 用拉格朗日方法重新推导单摆
# 验证：L = T - U → EL方程 → 与牛顿力学一致

include("../src/Lagrangian.jl")
using .Lagrangian
include("../../mini-newtonian/src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Pendulum via Lagrangian Mechanics")
    println("="^60)

    L_val = 1.0; m = 1.0; g = G_EARTH

    # ---- 1. 建立拉格朗日量 ----
    # 广义坐标 q = [θ]
    # 笛卡尔坐标: x = L*sin(θ), z = -L*cos(θ)  (z轴向上)
    # 动能: T = 0.5*m*(xdot^2 + zdot^2) = 0.5*m*L^2*θdot^2
    # 势能: U = m*g*z = -m*g*L*cos(θ)
    # L = 0.5*m*L^2*θdot^2 + m*g*L*cos(θ)

    function lagrangian(q, qdot)
        theta = q[1]
        theta_dot = qdot[1]
        T = 0.5 * m * L_val^2 * theta_dot^2
        U = -m * g * L_val * cos(theta)  # z=0 在悬挂点, 摆锤在下方时 z<0
        return T - U
    end

    # ∂L/∂θ = -m*g*L*sin(θ)
    function grad_L_q(q, qdot)
        return [-m * g * L_val * sin(q[1])]
    end

    # ∂L/∂θdot = m*L²*θdot
    function grad_L_qdot(q, qdot)
        return [m * L_val^2 * qdot[1]]
    end

    function M_func(q)
        return [m * L_val^2]
    end

    sys = EulerLagrangeSystem(1, lagrangian, grad_L_q, grad_L_qdot, M_func)
    f_ode = el_to_first_order(sys)

    # ---- 2. 与牛顿法对比 ----
    theta0 = deg2rad(30)
    omega0 = 0.0

    # 拉格朗日法
    y0 = [theta0, omega0]
    dt = 0.001; t_end = 5.0
    n_steps = Int(ceil(t_end/dt))
    y = copy(y0)
    traj_lag = Float64[]
    times = Float64[]
    for step in 1:n_steps
        if step % 10 == 1
            push!(times, (step-1)*dt)
            push!(traj_lag, y[1])
        end
        y = rk4_step(f_ode, (step-1)*dt, y, dt)
    end

    # 牛顿法（用 mini-newtonian 的 Verlet）
    function pend_accel(r, v, t)
        theta = r.x
        return Vec3(-g/L_val*sin(theta), 0, 0)
    end
    r0 = Vec3(theta0, 0, 0); v0 = Vec3(omega0, 0, 0)
    traj_newton = solve_fixed_step(pend_accel, r0, v0, t_end, dt; method=:verlet, record_every=10)

    # 对比终态
    theta_lag_final = traj_lag[end]
    theta_newton_final = traj_newton.positions[end].x
    println("θ_final (Lagrangian RK4): $(round(rad2deg(theta_lag_final),digits=4))°")
    println("θ_final (Newtonian Verlet): $(round(rad2deg(theta_newton_final),digits=4))°")
    println("Difference: $(round(abs(theta_lag_final - theta_newton_final), digits=8)) rad")

    # ---- 3. 能量守恒验证 ----
    E0 = energy_from_lagrangian(lagrangian, [theta0], [omega0], grad_L_qdot)
    E_final = energy_from_lagrangian(lagrangian, [y[1]], [y[2]], grad_L_qdot)
    println("\nEnergy conservation:")
    println("  E_initial = $(round(E0, digits=6))")
    println("  E_final   = $(round(E_final, digits=6))")
    println("  ΔE/E      = $(round((E_final-E0)/abs(E0), digits=10))")

    # ---- 4. 小角度周期验证 ----
    T_small = 2π * sqrt(L_val / g)
    println("\nSmall-angle period: $(round(T_small, digits=4)) s")

    println("\n✅ pendulum_lagrangian.jl done.")
end

main()
