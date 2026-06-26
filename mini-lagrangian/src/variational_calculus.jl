# variational_calculus.jl — 变分法基础与 Hamilton 原理
# 参考：Gelfand & Fomin "Calculus of Variations", Goldstein Ch.2, Landau Vol.1 Ch.1
# L2 Core Concepts: Hamilton's principle, functional derivatives
# L3 Mathematical Structures: Gateaux derivative, first/second variation

"""
Gateaux 导数 (泛函的方向导数):
  δF[y; η] = lim_{ε→0} (F[y + εη] - F[y]) / ε
对于作用量泛函 S[q], 这导出 Euler-Lagrange 方程.
"""
function gateaux_derivative(F::Function, y::Function, eta::Function; eps_val=1e-6)
    return (F(t -> y(t) + eps_val*eta(t)) - F(y)) / eps_val
end

"""
泛函导数 δF/δy(t): 满足 δF = ∫ (δF/δy(t)) η(t) dt.

对于 F[y] = ∫ L(t, y, y') dt:
  δF/δy = ∂L/∂y - d/dt (∂L/∂y')
此即 Euler-Lagrange 算子.
"""
function functional_derivative(L::Function, dL_dy::Function, dL_dydot::Function, y::Function, t::Float64)
    return dL_dy(t, y(t))
end

"""
Euler-Lagrange 方程 (变分法直接形式):
给定泛函 F[y] = ∫_a^b f(x, y, y') dx,
极值曲线满足: ∂f/∂y - d/dx (∂f/∂y') = 0.

此函数验证 y(x) 在 x 处的 EL 残差.
"""
function euler_lagrange_residual(f::Function, df_dy::Function, df_dydot::Function, y::Function, x::Float64; h=1e-6)
    dydot = df_dydot(x, y(x), (y(x+h)-y(x-h))/(2h))
    dydot_fwd = df_dydot(x+h, y(x+h), (y(x+2h)-y(x))/(2h))
    dydot_bwd = df_dydot(x-h, y(x-h), (y(x)-y(x-2h))/(2h))
    d_dx = (dydot_fwd - dydot_bwd) / (2h)
    dydx = df_dy(x, y(x), (y(x+h)-y(x-h))/(2h))
    return dydx - d_dx
end

"""
Hamilton 原理 (最小作用量原理):
在固定端点条件下, 真实路径使作用量 S[q] = ∫ L dt 取极值.

δS = 0, δq(t₁) = δq(t₂) = 0 ⇒ d/dt(∂L/∂q̇) - ∂L/∂q = 0
"""
function hamiltons_principle_action(L::Function, q_path::Function, qdot_path::Function, t_span::Tuple{Float64,Float64}; n=1000)
    return action(L, q_path, qdot_path, t_span, n)
end

"""
Maupertuis 原理 (最简作用量原理):
对于能量守恒系统, 约化作用量 S0 = ∫ p·dq = ∫ 2T dt 沿真实路径取极值.

与 Hamilton 原理的区别:
- Hamilton: 固定时间间隔 (t₁,t₂), δS = 0
- Maupertuis: 固定能量 E, δS₀ = 0 (路径可改变总时间)

S0 = ∫_{q₁}^{q₂} √(2m(E - U(q))) ds
其中 ds = √(dq·dq) 是位形空间的线元.
"""
function maupertuis_action(E::Float64, U::Function, q_path::AbstractVector{Vector{Float64}}, h::Float64)
    m = 1.0
    S0 = 0.0
    for i in 1:length(q_path)-1
        dq = q_path[i+1] - q_path[i]
        ds = norm(dq)
        q_mid = (q_path[i] + q_path[i+1]) / 2
        v = sqrt(2*max(0, E - U(q_mid)))
        S0 += v * ds
    end
    return S0
end

"""
Jacobi 度规 (Jacobi metric):
将 Maupertuis 原理表述为位形空间中的测地线原理.
度规张量: g_{ij}(q) = 2(E - U(q)) m_{ij}(q)
真实路径是此度规下的测地线.
"""
function jacobi_metric(E::Float64, U::Function, M::Matrix{Float64}, q::Vector{Float64})
    factor = 2 * max(0.0, E - U(q))
    return factor * M
end

"""
在 Jacobi 度规下两点间的测地线距离.
"""
function jacobi_geodesic_distance(E::Float64, U::Function, M_func::Function, q_path::AbstractVector{Vector{Float64}})
    dist = 0.0
    for i in 1:length(q_path)-1
        dq = q_path[i+1] - q_path[i]
        q_mid = (q_path[i] + q_path[i+1]) / 2
        g = jacobi_metric(E, U, M_func(q_mid), q_mid)
        ds_sq = dot(dq, g * dq)
        dist += sqrt(max(0, ds_sq))
    end
    return dist
end

"""
Beltrami 恒等式 (不显含自变量的情况):
若 L = L(y, y') 不显含 x, 则以下量为常数:
  L - y' ∂L/∂y' = const
"""
function beltrami_identity(L::Function, dL_dydot::Function, y::Float64, ydot::Float64)
    return L(y, ydot) - ydot * dL_dydot(y, ydot)
end

"""
验证 Beltrami 恒等式沿 Euler-Lagrange 解是否守恒.
"""
function verify_beltrami(L::Function, dL_dydot::Function, y_vals::Vector{Float64}, ydot_vals::Vector{Float64})
    n = length(y_vals)
    C = zeros(n)
    for i in 1:n
        C[i] = beltrami_identity(L, dL_dydot, y_vals[i], ydot_vals[i])
    end
    C0 = C[1]
    drift = (C .- C0) ./ max(abs(C0), 1e-300)
    return C, drift
end

"""
自然边界条件 (transversality condition):
若端点 y(a) 或 y(b) 自由, 则:
  ∂L/∂y'|_{x=a} = 0  或  ∂L/∂y'|_{x=b} = 0
"""
function transversality_condition(L::Function, dL_dydot::Function, y::Float64, ydot::Float64)
    return dL_dydot(y, ydot)
end

"""
等周问题 (Isoperimetric problem):
在约束 ∫ G(y, y') dx = const 下极小化 ∫ F(y, y') dx.
修正泛函: H = F + λ G
"""
function isoperimetric_el(F::Function, G::Function, dF_dy::Function, dG_dy::Function,
                          dF_dydot::Function, dG_dydot::Function, lambda::Float64, y::Function, x::Float64; h=1e-6)
    yp_fwd = (y(x+h) - y(x-h))/(2h)
    dydot_cur = dF_dydot(y(x), yp_fwd) + lambda*dG_dydot(y(x), yp_fwd)
    dydot_fwd = dF_dydot(y(x+h), (y(x+2h)-y(x))/(2h)) + lambda*dG_dydot(y(x+h), (y(x+2h)-y(x))/(2h))
    dydot_bwd = dF_dydot(y(x-h), (y(x)-y(x-2h))/(2h)) + lambda*dG_dydot(y(x-h), (y(x)-y(x-2h))/(2h))
    d_dx = (dydot_fwd - dydot_bwd) / (2h)
    dydx = dF_dy(y(x), yp_fwd) + lambda*dG_dy(y(x), yp_fwd)
    return abs(dydx - d_dx)
end

"""
含高阶导数的 Euler-Lagrange 方程:
F[y] = ∫ f(x, y, y', y'') dx
EL: ∂f/∂y - d/dx(∂f/∂y') + d²/dx²(∂f/∂y'') = 0
"""
function euler_lagrange_2nd_order(f::Function, df_dy::Function, df_dydot::Function, df_dy2dot::Function, y::Function, x::Float64; h=1e-6)
    yp = (y(x+h) - y(x-h))/(2h)
    ypp = (y(x+h) - 2*y(x) + y(x-h))/h^2
    term0 = df_dy(y(x), yp, ypp)
    term1_cur = df_dydot(y(x), yp, ypp)
    term1_fwd = df_dydot(y(x+h), (y(x+2h)-y(x))/(2h), (y(x+2h)-2*y(x+h)+y(x))/h^2)
    term1_bwd = df_dydot(y(x-h), (y(x)-y(x-2h))/(2h), (y(x)-2*y(x-h)+y(x-2h))/h^2)
    d_term1 = (term1_fwd - term1_bwd)/(2h)
    term2 = df_dy2dot(y(x), yp, ypp)
    term2_fwd = df_dy2dot(y(x+h), (y(x+2h)-y(x))/(2h), (y(x+2h)-2*y(x+h)+y(x))/h^2)
    term2_bwd = df_dy2dot(y(x-h), (y(x)-y(x-2h))/(2h), (y(x)-2*y(x-h)+y(x-2h))/h^2)
    d2_term2 = (term2_fwd - 2*term2 + term2_bwd)/h^2
    return term0 - d_term1 + d2_term2
end

"""
多变量变分问题: 多个因变量的 EL 方程组.
"""
function multi_variable_el_residual(L::Function, grad_q::Function, grad_qdot::Function,
                                     q_func::Function, qdot_func::Function, t::Float64, n_vars::Int; h=1e-6)
    residuals = zeros(n_vars)
    q = [q_func(t)[i] for i in 1:n_vars]
    qdot = [qdot_func(t)[i] for i in 1:n_vars]
    dL_dq = grad_q(q, qdot)
    dL_dqdot = grad_qdot(q, qdot)
    for i in 1:n_vars
        t_fwd = t + h; t_bwd = t - h
        q_fwd = [q_func(t_fwd)[j] for j in 1:n_vars]
        q_bwd = [q_func(t_bwd)[j] for j in 1:n_vars]
        qdot_fwd = [qdot_func(t_fwd)[j] for j in 1:n_vars]
        qdot_bwd = [qdot_func(t_bwd)[j] for j in 1:n_vars]
        dqdot_fwd = grad_qdot(q_fwd, qdot_fwd)[i]
        dqdot_bwd = grad_qdot(q_bwd, qdot_bwd)[i]
        d_dt = (dqdot_fwd - dqdot_bwd) / (2h)
        residuals[i] = dL_dq[i] - d_dt
    end
    return residuals
end

"""
场论的 Euler-Lagrange 方程:
S[φ] = ∫ L(φ, ∂_μ φ) d^4x
δS/δφ = ∂L/∂φ - ∂_μ (∂L/∂(∂_μ φ)) = 0
"""
function field_el_residual(L_density::Function, dL_dphi::Function, dL_ddphi::Function,
                            phi::Function, x::Float64, t::Float64; h=1e-6)
    phix = (phi(t, x+h) - phi(t, x-h))/(2h)
    phit = (phi(t+h, x) - phi(t-h, x))/(2h)
    term1 = dL_dphi(phi(t,x), phit, phix)
    dL_dphit = dL_ddphi(phi(t,x), phit, phix)[1]
    dL_dphit_fwd = dL_ddphi(phi(t+h,x), (phi(t+2h,x)-phi(t,x))/(2h), (phi(t+h,x+h)-phi(t+h,x-h))/(2h))[1]
    dL_dphit_bwd = dL_ddphi(phi(t-h,x), (phi(t,x)-phi(t-2h,x))/(2h), (phi(t-h,x+h)-phi(t-h,x-h))/(2h))[1]
    d_dt = (dL_dphit_fwd - dL_dphit_bwd)/(2h)
    dL_dphix = dL_ddphi(phi(t,x), phit, phix)[2]
    dL_dphix_fwd = dL_ddphi(phi(t,x+h), (phi(t+h,x+h)-phi(t-h,x+h))/(2h), (phi(t,x+2h)-phi(t,x))/(2h))[2]
    dL_dphix_bwd = dL_ddphi(phi(t,x-h), (phi(t+h,x-h)-phi(t-h,x-h))/(2h), (phi(t,x)-phi(t,x-2h))/(2h))[2]
    d_dx = (dL_dphix_fwd - dL_dphix_bwd)/(2h)
    return term1 - d_dt - d_dx
end
