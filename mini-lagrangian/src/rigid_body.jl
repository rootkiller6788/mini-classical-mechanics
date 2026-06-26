# rigid_body.jl — 刚体拉格朗日动力学
# 参考：Goldstein Ch.4-5, Landau Vol.1 Ch.6, MIT 8.012 Lecture 26-27
# L3 Mathematical Structures: SO(3) rotation group, Euler angles
# L4 Fundamental Laws: Euler equations
# L6 Canonical Systems: torque-free top, heavy top

"""
Euler 角 (z-x-z 约定):
  φ: 进动角 (绕空间 z 轴)   ∈ [0, 2π)
  θ: 章动角 (绕节线)        ∈ [0, π]
  ψ: 自转角 (绕物体 z 轴)   ∈ [0, 2π)

旋转矩阵: R(φ,θ,ψ) = R_z(φ) R_x(θ) R_z(ψ)
"""
struct EulerAngles
    phi::Float64
    theta::Float64
    psi::Float64
end

"""
绕 x 轴的旋转矩阵.
"""
function rot_x(angle::Float64)
    c,s = cos(angle), sin(angle)
    return [1 0 0; 0 c -s; 0 s c]
end

"""
绕 y 轴的旋转矩阵.
"""
function rot_y(angle::Float64)
    c,s = cos(angle), sin(angle)
    return [c 0 s; 0 1 0; -s 0 c]
end

"""
绕 z 轴的旋转矩阵.
"""
function rot_z(angle::Float64)
    c,s = cos(angle), sin(angle)
    return [c -s 0; s c 0; 0 0 1]
end

"""
从 Euler 角构造完整旋转矩阵 (z-x-z 约定).
R(φ,θ,ψ) = R_z(φ) · R_x(θ) · R_z(ψ)
"""
function rotation_matrix_euler(angles::EulerAngles)
    return rot_z(angles.phi) * rot_x(angles.theta) * rot_z(angles.psi)
end

"""
从空间坐标转换到物体坐标: r_body = R^T · r_space.
"""
function body_to_space(R::Matrix{Float64}, r_body::Vector{Float64})
    return R * r_body
end

function space_to_body(R::Matrix{Float64}, r_space::Vector{Float64})
    return R' * r_space
end

"""
刚体的转动惯量张量 (3×3 矩阵):
I_{ij} = ∫ ρ(r) (r²δ_{ij} - r_i r_j) d³r

对于离散点质量: I = Σ m_α (|r_α|² I - r_α r_α^T)
"""
function inertia_tensor_point(masses::Vector{Float64}, positions::Vector{Vector{Float64}})
    n = length(masses)
    I = zeros(3, 3)
    for alpha in 1:n
        m = masses[alpha]; r = positions[alpha]
        r_sq = dot(r, r)
        for i in 1:3, j in 1:3
            I[i,j] += m * ((i==j ? r_sq : 0.0) - r[i]*r[j])
        end
    end
    return I
end

"""
平行轴定理: I_{new} = I_cm + M (d²I - d d^T)
其中 d 是从质心到新参考点的位移矢量.
"""
function inertia_tensor_parallel_axis(I_cm::Matrix{Float64}, M::Float64, d::Vector{Float64})
    d_sq = dot(d, d)
    I_shift = M * (d_sq * Matrix{Float64}(I,3,3) - d * d')
    return I_cm + I_shift
end

"""
主转动惯量: 对角化 I 得到 I₁, I₂, I₃ 和主轴方向.
"""
function principal_inertias(I::Matrix{Float64})
    eigen_vals = eigen(Symmetric(I)).values
    return eigen_vals
end

"""
物体角速度与 Euler 角导数的关系:
ω_body = [φ̇ sinθ sinψ + θ̇ cosψ,
          φ̇ sinθ cosψ - θ̇ sinψ,
          φ̇ cosθ + ψ̇]
"""
function angular_velocity_body(angles::EulerAngles, dphi::Float64, dtheta::Float64, dpsi::Float64)
    sθ, cθ = sin(angles.theta), cos(angles.theta)
    sψ, cψ = sin(angles.psi), cos(angles.psi)
    wx = dphi*sθ*sψ + dtheta*cψ
    wy = dphi*sθ*cψ - dtheta*sψ
    wz = dphi*cθ + dpsi
    return [wx, wy, wz]
end

"""
刚体动能 (在主轴坐标系中): T = ½(I₁ ω₁² + I₂ ω₂² + I₃ ω₃²).
在一般坐标系: T = ½ ω^T I ω.
"""
function rigid_body_kinetic_energy(I_body::Vector{Float64}, omega::Vector{Float64})
    return 0.5 * (I_body[1]*omega[1]^2 + I_body[2]*omega[2]^2 + I_body[3]*omega[3]^2)
end

"""
Euler 方程 (无力矩情况): 
I₁ ω̇₁ = (I₂ - I₃) ω₂ ω₃
I₂ ω̇₂ = (I₃ - I₁) ω₃ ω₁
I₃ ω̇₃ = (I₁ - I₂) ω₁ ω₂
"""
function euler_equations_ode(I_body::Vector{Float64})
    I1, I2, I3 = I_body[1], I_body[2], I_body[3]
    function f(t, omega)
        return [
            (I2 - I3)/I1 * omega[2] * omega[3],
            (I3 - I1)/I2 * omega[3] * omega[1],
            (I1 - I2)/I3 * omega[1] * omega[2],
        ]
    end
    return f
end

"""
Euler 方程含外力矩 N:
I₁ ω̇₁ - (I₂ - I₃) ω₂ ω₃ = N₁
I₂ ω̇₂ - (I₃ - I₁) ω₃ ω₁ = N₂
I₃ ω̇₃ - (I₁ - I₂) ω₁ ω₂ = N₃
"""
function euler_equations_with_torque(I_body::Vector{Float64}, torque_func::Function)
    I1, I2, I3 = I_body[1], I_body[2], I_body[3]
    function f(t, omega)
        N = torque_func(t, omega)
        return [
            ((I2 - I3)*omega[2]*omega[3] + N[1]) / I1,
            ((I3 - I1)*omega[3]*omega[1] + N[2]) / I2,
            ((I1 - I2)*omega[1]*omega[2] + N[3]) / I3,
        ]
    end
    return f
end

"""
无力矩进动 (Poinsot 构造).

对于对称陀螺 I₁ = I₂ ≠ I₃:
  ω₃ = const
  ω₁ = A cos(Ω t + δ), ω₂ = A sin(Ω t + δ)
  进动频率 Ω = (I₃ - I₁)/I₁ ω₃

数值模拟无力矩自由进动.
"""
function torque_free_precession(I_body::Vector{Float64}, omega0::Vector{Float64}, t_end::Float64, h::Float64)
    f_ode = euler_equations_ode(I_body)
    y0 = copy(omega0)
    return rk4_integrate(f_ode, y0, (0.0, t_end), h)
end

"""
对称陀螺的解析进动频率 Ω_prec = (I₃ - I₁)/I₁ ω₃.
"""
function precession_frequency(I1::Float64, I3::Float64, omega3::Float64)
    return (I3 - I1) / I1 * omega3
end

"""
重力陀螺 (heavy top) 拉格朗日量.

固定点位于原点, 质心在物体 z 轴上距离 ℓ 处.
Euler 角: (φ, θ, ψ) — 进动, 章动, 自转.
L = ½I₁(φ̇² sin²θ + θ̇²) + ½I₃(φ̇ cosθ + ψ̇)² - Mgℓ cosθ

守恒量: p_φ, p_ψ (循环坐标), 能量 E.
"""
function heavy_top_lagrangian(I1::Float64, I3::Float64, M::Float64, g::Float64, l::Float64)
    function L_func(q, qdot)
        phi, theta, psi = q[1], q[2], q[3]
        dphi, dtheta, dpsi = qdot[1], qdot[2], qdot[3]
        T = 0.5*I1*(dphi^2*sin(theta)^2 + dtheta^2) + 0.5*I3*(dphi*cos(theta) + dpsi)^2
        V = M*g*l*cos(theta)
        return T - V
    end
    function grad_q(q, qdot)
        phi, theta, psi = q[1], q[2], q[3]
        dphi, dtheta, dpsi = qdot[1], qdot[2], qdot[3]
        st, ct = sin(theta), cos(theta)
        # ∂L/∂φ = 0 (循环)
        dL_dphi = 0.0
        # ∂L/∂θ = I₁ φ̇² sinθ cosθ - I₃(φ̇ cosθ+ψ̇) φ̇ sinθ + Mgℓ sinθ
        dL_dtheta = I1*dphi^2*st*ct - I3*(dphi*ct + dpsi)*dphi*st + M*g*l*st
        # ∂L/∂ψ = 0 (循环)
        dL_dpsi = 0.0
        return [dL_dphi, dL_dtheta, dL_dpsi]
    end
    function grad_qdot(q, qdot)
        phi, theta, psi = q[1], q[2], q[3]
        dphi, dtheta, dpsi = qdot[1], qdot[2], qdot[3]
        st, ct = sin(theta), cos(theta)
        p_phi = I1*dphi*st^2 + I3*(dphi*ct + dpsi)*ct
        p_theta = I1*dtheta
        p_psi = I3*(dphi*ct + dpsi)
        return [p_phi, p_theta, p_psi]
    end
    function M_func(q)
        st, ct = sin(q[2]), cos(q[2])
        return [I1*st^2+I3*ct^2 0 I3*ct; 0 I1 0; I3*ct 0 I3]
    end
    return L_func, grad_q, grad_qdot, M_func
end

"""
有效势 for heavy top (消除 φ, ψ 两个循环坐标):
U_eff(θ) = (p_φ - p_ψ cosθ)²/(2I₁ sin²θ) + Mgℓ cosθ

平衡点 dU_eff/dθ = 0 给出稳态进动条件.
"""
function heavy_top_effective_potential(theta::Float64, p_phi::Float64, p_psi::Float64, I1::Float64, I3::Float64, M::Float64, g::Float64, l::Float64)
    st, ct = sin(theta), cos(theta)
    centrifugal = (p_phi - p_psi*ct)^2 / (2*I1*st^2)
    gravitational = M*g*l*ct
    return centrifugal + gravitational
end

"""
欧拉角的运动学方程: dR/dt = R · [ω]_×
其中 [ω]_× 是角速度的反对称矩阵.
"""
function rotation_matrix_derivative(R::Matrix{Float64}, omega::Vector{Float64})
    omega_cross = [0 -omega[3] omega[2]; omega[3] 0 -omega[1]; -omega[2] omega[1] 0]
    return R * omega_cross
end

"""
四元数旋转表示 (避免万向节锁).
q = (w, x, y, z): w = cos(α/2), (x,y,z) = sin(α/2) n̂

四元数乘法: q1·q2.
"""
function quaternion_multiply(q1::Vector{Float64}, q2::Vector{Float64})
    w1,x1,y1,z1 = q1; w2,x2,y2,z2 = q2
    return [
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    ]
end

"""
四元数到旋转矩阵.
"""
function quaternion_to_matrix(q::Vector{Float64})
    w,x,y,z = q
    return [
        w^2+x^2-y^2-z^2  2*(x*y-w*z)      2*(x*z+w*y);
        2*(x*y+w*z)      w^2-x^2+y^2-z^2  2*(y*z-w*x);
        2*(x*z-w*y)      2*(y*z+w*x)      w^2-x^2-y^2+z^2
    ]
end

"""
四元数运动学: dq/dt = ½ q · (0, ω).
"""
function quaternion_kinematics(q::Vector{Float64}, omega::Vector{Float64})
    omega_q = [0.0, omega[1], omega[2], omega[3]]
    return 0.5 .* quaternion_multiply(q, omega_q)
end
