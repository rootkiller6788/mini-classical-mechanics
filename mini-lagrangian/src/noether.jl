# noether.jl — Noether 定理：对称性与守恒量
# 参考：Goldstein Ch.13, MIT 8.012 Lecture 23

"""
Noether 定理：
每个连续对称变换对应一个守恒量（Noether 荷）

若坐标变换 q -> q + ε*Q(q) 使 L 不变（到全导数），
则守恒量: C = (∂L/∂qdot) · Q(q)

常见对称性与守恒量：
  时间平移不变  → 能量守恒（Hamiltonian）
  空间平移不变  → 动量守恒
  旋转不变      → 角动量守恒
"""

"""
计算 Noether 荷
Q(q): 对称变换的生成元（矢量场）
grad_L_qdot(q, qdot): ∂L/∂qdot
返回守恒量 C = p · Q
"""
function noether_charge(grad_L_qdot::Vector{Float64}, Q::Vector{Float64})
    return dot(grad_L_qdot, Q)
end

"""
时间平移对称性 → 能量守恒
对于 L(q, qdot) 不明显含 t（∂L/∂t = 0），能量:
  E = qdot · ∂L/∂qdot - L
"""
function energy_from_lagrangian(L::Function, q::Vector{Float64}, qdot::Vector{Float64}, grad_L_qdot::Function)
    L_val = L(q, qdot)
    p = grad_L_qdot(q, qdot)
    return dot(qdot, p) - L_val
end

"""
空间平移对称性 → 动量守恒
变换: q_i -> q_i + ε （对所有粒子沿同一方向平移）
生成元 Q = (1,0,0,...) 表示沿第一个坐标平移
"""
function translation_symmetry_Q(n_particles::Int, direction::Int=1)
    # 生成沿 direction 方向（1=x, 2=y, 3=z）的整体平移
    n_dof = 3 * n_particles
    Q = zeros(n_dof)
    for i in 1:n_particles
        Q[3*(i-1) + direction] = 1.0
    end
    return Q
end

"""
旋转对称性 → 角动量守恒
绕 z 轴旋转: x -> x - ε*y, y -> y + ε*x
生成元: Q = (-y1, x1, 0, -y2, x2, 0, ..., -yN, xN, 0)
"""
function rotation_symmetry_Q_z(positions::Vector{Vector{Float64}})
    n_particles = length(positions)
    Q = zeros(3 * n_particles)
    for i in 1:n_particles
        x, y = positions[i][1], positions[i][2]
        Q[3*(i-1)+1] = -y
        Q[3*(i-1)+2] =  x
    end
    return Q
end

"""
验证守恒量是否在轨迹上保持不变
trajectory: 包含 q(t), qdot(t) 的数组
compute_charge: (q, qdot) -> Float64
"""
function verify_conservation(trajectory_q::Vector{Vector{Float64}}, trajectory_qdot::Vector{Vector{Float64}}, compute_charge::Function)
    n = length(trajectory_q)
    charges = [compute_charge(trajectory_q[i], trajectory_qdot[i]) for i in 1:n]
    C0 = charges[1]
    drift = [(c - C0) / (abs(C0) + 1e-300) for c in charges]
    return charges, drift
end

"""
洛伦兹 boost 对称性生成元 (Galilean boost 的相对论推广)
对于单粒子: Q = (t, x)
"""
function boost_symmetry_Q(t::Float64, x::Float64)
    [t, x]
end

"""
尺度变换对称性: q → e^ε q
若 L 是齐次的 (例如 Kepler 问题), 则存在 Runge-Lenz 守恒量
"""
function scale_symmetry_Q(positions::Vector{Vector{Float64}})
    n = length(positions)
    Q = zeros(3*n)
    for i in 1:n
        Q[3*(i-1)+1] = positions[i][1]
        Q[3*(i-1)+2] = positions[i][2]
        Q[3*(i-1)+3] = positions[i][3]
    end
    return Q
end

"""
场论 Noether 定理 (1+1D 标量场)
对于作用量 S[φ] = ∫L(φ, ∂_t φ, ∂_x φ) dtdx
守恒流 j^μ: ∂_μ j^μ = 0

空间平移 → 应力-能量张量 T^μ_ν
时间平移 → 能量守恒: ∂_t T^00 + ∂_x T^10 = 0
"""
function field_noether_current(L::Function, dL_dphit::Function, dL_dphix::Function, phi::Function, x::Float64, t::Float64; dxy=1e-3)
    phit = (phi(x, t+dxy) - phi(x, t-dxy))/(2*dxy)
    phix = (phi(x+dxy, t) - phi(x-dxy, t))/(2*dxy)
    T00 = dL_dphit(phit, phix)*phit - L(phit, phix)
    T01 = dL_dphit(phit, phix)*phix
    T10 = dL_dphix(phit, phix)*phit
    T11 = dL_dphix(phit, phix)*phix - L(phit, phix)
    return (T00=T00, T01=T01, T10=T10, T11=T11)
end

"""
对称性代数闭合性检查: [Q_a, Q_b] = f_{ab}^c Q_c
对于旋转生成元: [L_x, L_y] = L_z
"""
function lie_algebra_closure(Q1::Function, Q2::Function, q::Vector{Float64}; eps_val=1e-5)
    n = length(q)
    Q1_val = Q1(q); Q2_val = Q2(q)
    dQ2_dq = zeros(n, n)
    for i in 1:n
        qp = copy(q); qp[i] += eps_val; qm = copy(q); qm[i] -= eps_val
        dQ2_dq[:, i] = (Q2(qp) - Q2(qm))/(2*eps_val)
    end
    commutator = dQ2_dq * Q1_val - zeros(n,n) * Q2_val  # [Q1, Q2] = dQ2·Q1
    return commutator
end

"""
验证拉格朗日在对称变换下的不变性 (至全导数)
L(q+εQ, q̇+εQ̇) - L(q, q̇) = ε dF/dt
"""
function verify_symmetry_invariance(L::Function, grad_L_q::Function, grad_L_qdot::Function,
                                     Q::Vector{Float64}, Qdot::Vector{Float64}, q::Vector{Float64}, qdot::Vector{Float64})
    dL = dot(grad_L_q(q, qdot), Q) + dot(grad_L_qdot(q, qdot), Qdot)
    return abs(dL) < 1e-10
end

## ============================================================
## 电磁场 Noether 定理与守恒流
## ============================================================

"""
电磁场的守恒流来自规范不变性和时空对称性.
自由电磁场拉格朗日密度: L = -1/4 F_{mu nu} F^{mu nu}

时间平移 -> 能量守恒 -> Poynting 矢量
空间平移 -> 动量守恒 -> Maxwell 应力张量
旋转变换 -> 角动量守恒
"""
function poynting_vector(E::Vector{Float64}, B::Vector{Float64}; mu0=1.25663706212e-6)
    return (1/mu0) .* cross(E, B)
end

function maxwell_stress_tensor(E::Vector{Float64}, B::Vector{Float64}; eps0=8.854187817e-12, mu0=1.25663706212e-6)
    # T_{ij} = eps0 (E_i E_j - 1/2 delta_{ij} E^2) + 1/mu0 (B_i B_j - 1/2 delta_{ij} B^2)
    T = zeros(3, 3)
    E_sq = dot(E, E); B_sq = dot(B, B)
    for i in 1:3, j in 1:3
        T[i,j] = eps0*(E[i]*E[j] - 0.5*(i==j ? 1 : 0)*E_sq) + (1/mu0)*(B[i]*B[j] - 0.5*(i==j ? 1 : 0)*B_sq)
    end
    return T
end

"""
质量守恒与连续方程: partial_t rho + div J_mass = 0
拉格朗日框架下: 物质守恒是 Noether 定理在 U(1) 相位对称性下的结果.
"""
function continuity_check(rho_func::Function, J_func::Function, x::Float64, t::Float64; h=1e-6)
    drho_dt = (rho_func(x, t+h) - rho_func(x, t-h))/(2*h)
    dJ_dx = (J_func(x+h, t) - J_func(x-h, t))/(2*h)
    residual = drho_dt + dJ_dx
    return residual
end

"""
共形对称性: 尺度变换 x -> lambda x, t -> lambda t.
对于 T^{mu}_{mu} = 0 的理论 (如电磁场), 有标度不变性.
对于有质量粒子, 迹反常 = m^2 phi^2 打破标度不变性.
"""
function trace_anomaly(T00::Float64, T11::Float64, T22::Float64, T33::Float64; signature=:minkowski)
    if signature == :minkowski
        return -T00 + T11 + T22 + T33  # (-,+,+,+) signature
    else
        return T00 + T11 + T22 + T33   # (+,+,+,+) signature
    end
end

"""
螺旋度守恒 (理想流体): H = int v · omega d^3x
其中 omega = curl(v). 在非粘性正压流体中, 螺旋度是 Noether 不变量.
"""
function helicity_density(v::Vector{Float64}, omega::Vector{Float64})
    return dot(v, omega)
end

"""
磁螺旋守恒 (理想 MHD): H_M = int A · B d^3x
其中 B = curl(A). 在理想磁流体力学中守恒.
"""
function magnetic_helicity(A::Function, B::Function, bounds::Vector{Float64}; n_grid=30)
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    dx = (x_max-x_min)/n_grid; dy = (y_max-y_min)/n_grid; dz = (z_max-z_min)/n_grid
    H = 0.0
    r = zeros(3)
    for i in 1:n_grid, j in 1:n_grid, k in 1:n_grid
        r[1] = x_min + (i-0.5)*dx; r[2] = y_min + (j-0.5)*dy; r[3] = z_min + (k-0.5)*dz
        A_val = A(r, 0.0); B_val = B(r, 0.0)
        H += dot(A_val, B_val) * dx * dy * dz
    end
    return H
end
