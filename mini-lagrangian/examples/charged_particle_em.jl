#!/usr/bin/env julia
# examples/charged_particle_em.jl — 电磁场中的带电粒子（速度依赖势）
# 拉格朗日方法的标志性应用：L = 0.5*m*v^2 - q*(φ - v·A)
# ∂L/∂v = mv + qA = p_canonical（正则动量 ≠ 机械动量）

include("../src/Lagrangian.jl")
using .Lagrangian
include("../../mini-newtonian/src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Charged Particle in EM Field (Lagrangian)")
    println("="^60)

    m=1.0; q=1.0

    # 均匀磁场 B = B0*ẑ, 无电场 E=0
    # 矢势: A = (-B0*y/2, B0*x/2, 0)（对称规范）
    # 标势: φ = 0
    B0 = 1.0

    function A_vec(r)
        return Vec3(-B0*r.y/2, B0*r.x/2, 0.0)
    end

    function phi(r)
        return 0.0
    end

    # L = 0.5*m*v² + q*v·A - q*φ
    function em_lagrangian(q_coords, qdot)
        r = Vec3(q_coords[1], q_coords[2], q_coords[3])
        v = Vec3(qdot[1], qdot[2], qdot[3])
        A = A_vec(r)
        T = 0.5 * m * norm2(v)
        U_em = q * (phi(r) - dot(v, A))  # 广义势
        return T - U_em
    end

    function em_grad_L_q(q_coords, qdot)
        r = Vec3(q_coords[1], q_coords[2], q_coords[3])
        v = Vec3(qdot[1], qdot[2], qdot[3])
        # ∂L/∂r = q*∂(v·A)/∂r - q*∂φ/∂r
        # 对 A = (-B0*y/2, B0*x/2, 0): v·A = B0/2*(-vx*y + vy*x)
        # ∂(v·A)/∂x = B0*vy/2, ∂(v·A)/∂y = -B0*vx/2, ∂(v·A)/∂z = 0
        return [q*B0*qdot[2]/2, -q*B0*qdot[1]/2, 0.0]
    end

    function em_grad_L_qdot(q_coords, qdot)
        r = Vec3(q_coords[1], q_coords[2], q_coords[3])
        v = Vec3(qdot[1], qdot[2], qdot[3])
        A = A_vec(r)
        # p_canonical = m*v + q*A
        p = m * v + q * A
        return [p.x, p.y, p.z]
    end

    function em_M(q_coords)
        return [m 0 0; 0 m 0; 0 0 m]
    end

    sys = EulerLagrangeSystem(3, em_lagrangian, em_grad_L_q, em_grad_L_qdot, em_M)
    f_ode = el_to_first_order(sys)

    # 初始条件：v 在 xy 平面，B 沿 z → 圆周运动
    v0_mag = 5.0
    y0 = [1.0, 0.0, 0.0, 0.0, v0_mag, 0.0]  # [x,y,z, vx,vy,vz]
    dt = 0.001; t_end = 10.0
    n_steps = Int(ceil(t_end/dt))

    y = copy(y0)
    record_every = 50
    xs, ys = Float64[], Float64[]

    for step in 1:n_steps
        if step % record_every == 1
            push!(xs, y[1]); push!(ys, y[2])
        end
        y = rk4_step(f_ode, (step-1)*dt, y, dt)
    end

    # 分析
    println("Initial: r=(1,0,0), v=(0,$v0_mag,0)")
    println("Final: r=($(round(y[1],digits=3)),$(round(y[2],digits=3)),$(round(y[3],digits=3)))")

    # 回旋半径理论值: r_cycl = m*v / (|q|*B)
    r_cycl_theory = m * v0_mag / (abs(q) * B0)
    # 从轨迹估算半径
    r_avg = mean([sqrt(x^2 + y^2) for (x,y) in zip(xs, ys)])
    println("Cyclotron radius: theory=$(round(r_cycl_theory,digits=3)), sim avg=$(round(r_avg,digits=3))")

    # 回旋频率: ω_c = |q|*B / m
    omega_c = abs(q) * B0 / m
    println("Cyclotron frequency: $(round(omega_c,digits=3)) rad/s")

    # 验证正则动量守恒（磁场沿 z，p_z 守恒）
    r_final = Vec3(y[1], y[2], y[3])
    v_final = Vec3(y[4], y[5], y[6])
    A_final = A_vec(r_final)
    p0 = m * Vec3(0, v0_mag, 0) + q * A_vec(Vec3(1,0,0))
    p_final = m * v_final + q * A_final
    println("\nCanonical momentum p_z: initial=$(round(p0.z,digits=6)), final=$(round(p_final.z,digits=6))")
    println("(p_z conserved: $(abs(p0.z - p_final.z) < 1e-3 ? "YES" : "NO"))")

    # 机械动量不守恒！（因为磁场力做... 不对，磁场力不做功）
    # 但 mv 方向不断改变
    println("Kinetic energy: initial=$(round(0.5*m*v0_mag^2,digits=2)), final=$(round(0.5*m*norm2(v_final),digits=2))")
    println("(KE conserved: $(abs(0.5*m*v0_mag^2 - 0.5*m*norm2(v_final)) < 1e-2 ? "YES" : "NO"))")

    println("\n✅ charged_particle_em.jl done.")
end

main()
