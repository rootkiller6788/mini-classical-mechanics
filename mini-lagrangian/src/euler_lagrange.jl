# euler_lagrange.jl — Euler-Lagrange 方程求解器
# d/dt(∂L/∂qdot) - ∂L/∂q = 0

using LinearAlgebra

"""
Euler-Lagrange 方程数值求解（通用形式）

输入：
  lagrangian(q, qdot) -> Float64             拉格朗日量 L(q, qdot)
  grad_L_q(q, qdot) -> Vector{Float64}       ∂L/∂q
  grad_L_qdot(q, qdot) -> Vector{Float64}    ∂L/∂qdot
  M_inv_func(q) -> Matrix{Float64}           质量矩阵的逆

EL 方程转化为一阶 ODE:
  dq/dt = qdot
  d(qdot)/dt = M^{-1} * (∂L/∂q - dM/dt * qdot)
"""
struct EulerLagrangeSystem
    n_dof::Int
    lagrangian::Function           # L(q, qdot) -> Float64
    grad_L_q::Function             # ∂L/∂q (q, qdot) -> Vector
    grad_L_qdot::Function          # ∂L/∂qdot (q, qdot) -> Vector
    mass_matrix::Function          # M(q) -> Matrix
end

"""
将 EL 方程转换为一阶 ODE: [q, qdot] -> [qdot, qdotdot]
"""
function el_to_first_order(sys::EulerLagrangeSystem)
    function f(t, y)
        n = sys.n_dof
        q = y[1:n]
        qdot = y[n+1:2n]

        dL_dq = sys.grad_L_q(q, qdot)
        dL_dqdot = sys.grad_L_qdot(q, qdot)

        # 需要 d/dt(∂L/∂qdot) = M*qddot + Mdot*qdot
        # 数值上我们直接用有限差分或解析 Mdot
        # 对于大多数力学系统 Mdot*qdot 可以用 Christoffel 符号表示
        # 简化：用数值微分估算 d/dt(∂L/∂qdot)
        # 更实际的做法：直接解 M*qddot = dL/dq - Mdot*qdot

        M = sys.mass_matrix(q)
        # 数值估算 Mdot（如果 M 依赖 q）
        # Mdot_{ij} = Σ_k ∂M_{ij}/∂q_k * qdot_k
        Mdot = zeros(n, n)
        eps_val = 1e-6
        for i in 1:n, j in 1:n
            for k in 1:n
                q_plus = copy(q); q_plus[k] += eps_val
                q_minus = copy(q); q_minus[k] -= eps_val
                M_plus = sys.mass_matrix(q_plus)
                M_minus = sys.mass_matrix(q_minus)
                Mdot[i,j] += (M_plus[i,j] - M_minus[i,j]) / (2*eps_val) * qdot[k]
            end
        end

        rhs = dL_dq - Mdot * qdot
        qddot = M \ rhs

        return [qdot; qddot]
    end
    return f
end

"""
简化版 EL 求解器：用于 L = T(qdot,q) - U(q) 的标准力学系统
T = 0.5*qdot'*M(q)*qdot, U = U(q)

此时:
  ∂L/∂qdot = M*qdot
  ∂L/∂q = 0.5*∂(qdot'*M*qdot)/∂q - ∂U/∂q
  d/dt(∂L/∂qdot) = M*qddot + Mdot*qdot

EL: M*qddot = 0.5*grad(qdot'*M*qdot) - Mdot*qdot - grad(U)
"""
struct StandardLagrangian
    n_dof::Int
    M_func::Function              # M(q) -> Matrix
    grad_M_func::Function         # grad of each M_ij component
    grad_U::Function              # grad U(q) -> Vector
end

function standard_el_to_ode(sys::StandardLagrangian)
    function f(t, y)
        n = sys.n_dof
        q = y[1:n]
        qdot = y[n+1:2n]

        M = sys.M_func(q)
        dU = sys.grad_U(q)

        # 计算 Mdot
        Mdot = zeros(n, n)
        eps_val = 1e-6
        for i in 1:n, j in 1:n
            for k in 1:n
                qp = copy(q); qp[k] += eps_val
                qm = copy(q); qm[k] -= eps_val
                Mp = sys.M_func(qp)
                Mm = sys.M_func(qm)
                Mdot[i,j] += (Mp[i,j] - Mm[i,j]) / (2*eps_val) * qdot[k]
            end
        end

        # 动能对 q 的梯度：∂T/∂q_i = 0.5*qdot'*(∂M/∂q_i)*qdot
        dT_dq = zeros(n)
        for i in 1:n
            qp = copy(q); qp[i] += eps_val
            qm = copy(q); qm[i] -= eps_val
            Mp = sys.M_func(qp)
            Mm = sys.M_func(qm)
            dM_dqi = (Mp - Mm) / (2*eps_val)
            dT_dq[i] = 0.5 * dot(qdot, dM_dqi * qdot)
        end

        rhs = dT_dq - dU - Mdot * qdot
        qddot = M \ rhs

        return [qdot; qddot]
    end
    return f
end

## ============================================================
## 含非保守力的 Lagrange 方程
## ============================================================

"""
含耗散函数的 Lagrange 方程 (Rayleigh 耗散函数).
对于速度依赖的耗散力 F_i^diss = -partial R/partial qdot_i:
  d/dt(partial L/partial qdot_i) - partial L/partial q_i = -partial R/partial qdot_i
R = 1/2 c qdot^2 给出粘性阻尼力 -c*qdot.
"""
function rayleigh_dissipation_lagrange(R_func::Function, grad_R_qdot::Function, sys::EulerLagrangeSystem)
    function f(t, y)
        n = sys.n_dof
        q = y[1:n]; qdot = y[n+1:2n]
        dL_dq = sys.grad_L_q(q, qdot)
        dL_dqdot = sys.grad_L_qdot(q, qdot)
        M = sys.mass_matrix(q)
        # Mdot numerical
        Mdot = zeros(n, n); eps_val = 1e-6
        for i in 1:n, j in 1:n
            for k in 1:n
                qp = copy(q); qp[k] += eps_val; qm = copy(q); qm[k] -= eps_val
                Mdot[i,j] += (sys.mass_matrix(qp)[i,j] - sys.mass_matrix(qm)[i,j])/(2*eps_val)*qdot[k]
            end
        end
        diss = grad_R_qdot(q, qdot)
        rhs = dL_dq - diss - Mdot*qdot
        qddot = M \ rhs
        return [qdot; qddot]
    end
    return f
end

"""
含时拉格朗日量的 EL 方程.
当 L = L(q, qdot, t) 显含时间时, 总能量不守恒:
  dE/dt = -partial L/partial t
"""
function time_dependent_el_to_ode(L_func::Function, grad_L_q::Function, grad_L_qdot::Function, M_func::Function)
    function f(t, y)
        n = length(y) ÷ 2
        q = y[1:n]; qdot = y[n+1:2n]
        M = M_func(q)
        dL_dq = grad_L_q(q, qdot, t)
        dL_dqdot = grad_L_qdot(q, qdot, t)
        # Mdot via finite differences
        Mdot = zeros(n, n); eps_val = 1e-6
        for i in 1:n, j in 1:n
            for k in 1:n
                qp = copy(q); qp[k] += eps_val; qm = copy(q); qm[k] -= eps_val
                Mdot[i,j] += (M_func(qp)[i,j] - M_func(qm)[i,j])/(2*eps_val)*qdot[k]
            end
        end
        rhs = dL_dq - Mdot*qdot
        qddot = M \ rhs
        return [qdot; qddot]
    end
    return f
end

"""
能量变化率: dE/dt = -partial L/partial t (显含时系统的能量不守恒定理).
"""
function energy_change_rate(L_func::Function, q::Vector{Float64}, qdot::Vector{Float64}, t::Float64; h=1e-6)
    Lp = L_func(q, qdot, t+h)
    Lm = L_func(q, qdot, t-h)
    return -(Lp - Lm)/(2*h)
end
