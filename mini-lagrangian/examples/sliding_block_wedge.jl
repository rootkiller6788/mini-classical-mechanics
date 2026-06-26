#!/usr/bin/env julia
# examples/sliding_block_wedge.jl — 滑块在可动斜面上
# 经典约束系统：用拉格朗日乘子或广义坐标两种方法

include("../src/Lagrangian.jl")
using .Lagrangian
include("../../mini-newtonian/src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Sliding Block on Movable Wedge")
    println("="^60)

    M=5.0; m=1.0; theta=deg2rad(30); g=G_EARTH; mu=0.0

    # 系统：质量M的楔子可在地面上无摩擦滑动，质量m的滑块在楔子上滑动
    # 广义坐标: q1 = X (楔子水平位移), q2 = s (滑块沿斜面的位移)
    # 滑块笛卡尔坐标: x = X + s*cosθ, y = s*sinθ
    # T = 0.5*M*Xdot² + 0.5*m*((Xdot+sdot*cosθ)² + (sdot*sinθ)²)
    # U = m*g*s*sinθ
    # L = T - U

    function wedge_L(q, qdot)
        X, s = q[1], q[2]; Xd, sd = qdot[1], qdot[2]
        T = 0.5*M*Xd^2 + 0.5*m*((Xd+sd*cos(theta))^2 + (sd*sin(theta))^2)
        U = m*g*s*sin(theta)
        return T - U
    end

    function wedge_grad_q(q, qdot)
        # ∂L/∂X = 0 (X 是循环坐标!) → 守恒量
        # ∂L/∂s = -m*g*sinθ
        return [0.0, -m*g*sin(theta)]
    end

    function wedge_grad_qdot(q, qdot)
        X, s = q[1], q[2]; Xd, sd = qdot[1], qdot[2]
        pX = M*Xd + m*(Xd + sd*cos(theta))  # 总水平动量！
        ps = m*(Xd*cos(theta) + sd)
        return [pX, ps]
    end

    function wedge_M(q)
        return [(M+m) m*cos(theta); m*cos(theta) m]
    end

    sys = EulerLagrangeSystem(2, wedge_L, wedge_grad_q, wedge_grad_qdot, wedge_M)
    f = el_to_first_order(sys)

    # 初始条件：静止释放
    y0 = [0.0, 0.0, 0.0, 0.0]
    dt=0.001; t_end=2.0
    n=Int(ceil(t_end/dt)); y=copy(y0)
    for _ in 1:n; y=rk4_step(f,0.0,y,dt); end

    a_s = y[4]/t_end  # 滑块沿斜面的平均加速度
    a_X = y[3]/t_end  # 楔子的平均加速度

    # 理论值（从 EL 方程推导）:
    # sddot = (M+m)*g*sinθ / (M + m*sin²θ)
    # Xddot = -m*g*sinθ*cosθ / (M + m*sin²θ)
    denom = M + m*sin(theta)^2
    sddot_theory = (M+m)*g*sin(theta) / denom
    Xddot_theory = -m*g*sin(theta)*cos(theta) / denom

    println("Block acceleration along wedge:")
    println("  Numerical: $(round(a_s,digits=4)) m/s²")
    println("  Theory:    $(round(sddot_theory,digits=4)) m/s²")
    println("\nWedge acceleration (recoil):")
    println("  Numerical: $(round(a_X,digits=4)) m/s²")
    println("  Theory:    $(round(Xddot_theory,digits=4)) m/s²")

    # 验证水平动量守恒 (X 是循环坐标)
    pX_initial = wedge_grad_qdot([0,0],[0,0])[1]
    pX_final = wedge_grad_qdot([y[1],y[2]],[y[3],y[4]])[1]
    println("\nHorizontal momentum conservation:")
    println("  pX_initial: $pX_initial, pX_final: $(round(pX_final,digits=6))")
    println("  Conserved: $(abs(pX_final - pX_initial) < 1e-6)")

    # 能量守恒
    E0 = wedge_L([0,0],[0,0])
    E_f = wedge_L([y[1],y[2]],[y[3],y[4]])
    println("\nEnergy: drift = $(round((E_f-E0)/max(abs(E0),1e-300),digits=10))")

    println("\n✅ sliding_block_wedge.jl done.")
end
main()
