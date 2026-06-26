# legendre.jl — Legendre 变换：L → H
# 参考：Goldstein Ch.8, MIT 8.012 Lecture 23

"""
Legendre 变换: 从拉格朗日量到哈密顿量
给定 L(q, qdot)，定义 p = ∂L/∂qdot
解出 qdot = qdot(q, p)（需要 ∂²L/∂qdot² 可逆）
H(q, p) = p·qdot - L(q, qdot)

对于标准力学系统 T = 1/2*qdot'*M*qdot:
  p = M*qdot → qdot = M^{-1}*p
  H = 1/2*p'*M^{-1}*p + U(q) = T + U = E
"""

struct LegendreTransform
    n_dof::Int
    M_func::Function              # M(q) 质量矩阵
    M_inv_func::Function          # M^{-1}(q)
    U_func::Function              # U(q) 势能
    grad_U::Function              # ∇U(q)
end

"""
构造 L = T - U 的标准 Legendre 变换
"""
function standard_legendre(M_func::Function, U_func::Function, grad_U::Function, n::Int)
    function M_inv(q)
        M = M_func(q)
        return inv(M)
    end
    return LegendreTransform(n, M_func, M_inv, U_func, grad_U)
end

"""
qdot → p: p = M(q)*qdot
"""
function qdot_to_momentum(lt::LegendreTransform, q::Vector{Float64}, qdot::Vector{Float64})
    M = lt.M_func(q)
    return M * qdot
end

"""
p → qdot: qdot = M^{-1}(q)*p
"""
function momentum_to_qdot(lt::LegendreTransform, q::Vector{Float64}, p::Vector{Float64})
    Minv = lt.M_inv_func(q)
    return Minv * p
end

"""
计算 H(q, p)
"""
function hamiltonian(lt::LegendreTransform, q::Vector{Float64}, p::Vector{Float64})
    Minv = lt.M_inv_func(q)
    qdot = Minv * p
    T = 0.5 * dot(p, qdot)  # 0.5*p'*M^{-1}*p
    U = lt.U_func(q)
    return T + U
end

"""
Hamilton 方程: dq/dt = ∂H/∂p, dp/dt = -∂H/∂q
转换为一阶 ODE: y = [q; p]
"""
function hamiltons_equations(lt::LegendreTransform)
    function f(t, y)
        n = lt.n_dof
        q = y[1:n]
        p = y[n+1:2n]
        Minv = lt.M_inv_func(q)
        qdot = Minv * p

        # dp/dt = -∂H/∂q = -∂T/∂q - ∂U/∂q
        # ∂T/∂q_i = -0.5*p'*M^{-1}*(∂M/∂q_i)*M^{-1}*p
        dT_dq = zeros(n)
        eps_val = 1e-6
        for i in 1:n
            qp = copy(q); qp[i] += eps_val
            qm = copy(q); qm[i] -= eps_val
            Mp_inv = lt.M_inv_func(qp)
            Mm_inv = lt.M_inv_func(qm)
            Tp = 0.5 * dot(p, Mp_inv * p)
            Tm = 0.5 * dot(p, Mm_inv * p)
            dT_dq[i] = (Tp - Tm) / (2*eps_val)
        end

        dU_dq = lt.grad_U(q)
        pdot = -dT_dq - dU_dq

        return [qdot; pdot]
    end
    return f
end

"""
Poisson 括号: {f, g} = Σ (∂f/∂q_i)*(∂g/∂p_i) - (∂f/∂p_i)*(∂g/∂q_i)
数值计算（有限差分）
"""
function poisson_bracket(f::Function, g::Function, q::Vector{Float64}, p::Vector{Float64}; eps_val=1e-6)
    n = length(q)
    result = 0.0
    for i in 1:n
        qp = copy(q); qp[i] += eps_val
        qm = copy(q); qm[i] -= eps_val
        pp = copy(p); pp[i] += eps_val
        pm = copy(p); pm[i] -= eps_val

        df_dqi = (f(qp, p) - f(qm, p)) / (2*eps_val)
        dg_dpi = (g(q, pp) - g(q, pm)) / (2*eps_val)
        df_dpi = (f(q, pp) - f(q, pm)) / (2*eps_val)
        dg_dqi = (g(qp, p) - g(qm, p)) / (2*eps_val)

        result += df_dqi * dg_dpi - df_dpi * dg_dqi
    end
    return result
end

"""
验证 {q_i, p_j} = δ_ij（基本 Poisson 括号）
"""
function verify_canonical_poisson(q::Vector{Float64}, p::Vector{Float64})
    n = length(q)
    errors = zeros(n, n)
    for i in 1:n, j in 1:n
        qi_func(qq, pp) = qq[i]
        pj_func(qq, pp) = pp[j]
        pb = poisson_bracket(qi_func, pj_func, q, p)
        expected = i == j ? 1.0 : 0.0
        errors[i,j] = abs(pb - expected)
    end
    return errors
end

"""
Routhian: 混合 L 和 H
对于循环坐标（∂L/∂q_k = 0），用 p_k 替代 qdot_k
R(q_1..q_m, qdot_1..qdot_m, p_{m+1}..p_n, t) = Σ_{k>m} p_k*qdot_k - L
"""
function routhian(L::Function, grad_L_qdot::Function, q::Vector{Float64}, qdot::Vector{Float64}, cyclic_indices::Vector{Int})
    n = length(q)
    p_full = grad_L_qdot(q, qdot)
    L_val = L(q, qdot)
    R = L_val  # 符号约定因作者而异，此处用 Goldstein 的符号
    for k in cyclic_indices
        R -= p_full[k] * qdot[k]  # R = L - Σ p_k*qdot_k for cyclic coords
    end
    return R
end

"""
Routhian 运动方程: 对非循环坐标用 Lagrange, 对循环坐标用 Hamilton.

d/dt(∂R/∂q̇_i) = ∂R/∂q_i  (i ≤ m, 非循环)
ṗ_k = -∂R/∂q_k = 0        (k > m, 循环→p_k=const)
q̇_k = ∂R/∂p_k             (k > m)
"""
function routhian_equations(R::Function, grad_R_q::Function, grad_R_qdot::Function, grad_R_p::Function, state::Vector{Float64}, m::Int)
    n_cyclic = (length(state) - m) ÷ 2
    q = state[1:m]; qdot = state[m+1:2m]
    p_cyclic = state[2m+1:2m+n_cyclic]; q_cyclic = state[2m+n_cyclic+1:end]
    dq = copy(qdot)
    dqdot = grad_R_q(q, qdot, p_cyclic) - zeros(length(q))  # simplified
    dp_cyclic = -grad_R_p(q, qdot, p_cyclic)
    dq_cyclic = grad_R_p(q, qdot, p_cyclic)
    return [dq; dqdot; dp_cyclic; dq_cyclic]
end

"""
Legendre 变换的正则性: det(∂²L/∂q̇²) ≠ 0.
"""
function is_legendre_regular(L::Function, q::Vector{Float64}, qdot::Vector{Float64}; eps_val=1e-5)
    n = length(qdot)
    H = zeros(n, n)
    for i in 1:n
        for j in 1:n
            qpp = copy(qdot); qpp[i]+=eps_val; qpp[j]+=eps_val
            qpm = copy(qdot); qpm[i]+=eps_val; qpm[j]-=eps_val
            qmp = copy(qdot); qmp[i]-=eps_val; qmp[j]+=eps_val
            qmm = copy(qdot); qmm[i]-=eps_val; qmm[j]-=eps_val
            H[i,j] = (L(q,qpp) - L(q,qpm) - L(q,qmp) + L(q,qmm))/(4*eps_val^2)
        end
    end
    return det(H) > 1e-10
end

"""
Fenchel 共轭: f*(p) = sup_x (p·x - f(x)).
Legendre 变换是 Fenchel 共轭的特殊情况 (当 f 是凸且光滑).
"""
function fenchel_conjugate(f::Function, p::Float64, x_range::AbstractRange)
    maximum(p*x - f(x) for x in x_range)
end
