# action.jl — 作用量原理与变分积分器
# 参考：Goldstein Ch.2, Marsden & West (2001)

"""Action functional S[q] = integral L dt via trapezoidal rule."""
function action(lagrangian::Function, q_func::Function, qdot_func::Function, t_span::Tuple{Float64,Float64}, n_points::Int=1000)
    t0, tf = t_span; dt = (tf-t0)/(n_points-1); S = 0.0
    for i in 1:n_points
        t = t0 + (i-1)*dt; w = (i==1||i==n_points) ? 0.5 : 1.0
        S += w * lagrangian(q_func(t), qdot_func(t)) * dt
    end
    return S
end

"""Midpoint discrete Lagrangian: L_d(q0,q1,h) = h * L(midpoint, secant_velocity)."""
function midpoint_discrete_lagrangian(L::Function, q0::Vector{Float64}, q1::Vector{Float64}, h::Float64)
    return h * L((q0+q1)/2, (q1-q0)/h)
end

"""
Variational midpoint integrator step.
Solves DEL: D2_Ld(q_prev,q_cur) + D1_Ld(q_cur,q_next) = 0 for q_next.
For L = T - U with constant M, this is equivalent to implicit midpoint.
"""
function variational_midpoint_step(L::Function, grad_L_q::Function, grad_L_qdot::Function,
                                   M_func::Function, q_prev::Vector{Float64}, q_cur::Vector{Float64}, h::Float64;
                                   tol=1e-10, max_iter=20)
    n = length(q_cur)
    M_cur = M_func(q_cur)
    qdot_cur = (q_cur - q_prev) / h
    q_mid_cur = (q_prev + q_cur) / 2
    dL_dq = grad_L_q(q_mid_cur, qdot_cur)
    dL_dqdot = grad_L_qdot(q_mid_cur, qdot_cur)
    # Initial guess from discrete Legendre transform
    rhs = dL_dqdot - h/2 * dL_dq
    q_next = q_cur + h * (M_cur \ rhs)
    for _ in 1:max_iter
        q_mid_next = (q_cur + q_next) / 2
        qdot_next = (q_next - q_cur) / h
        dL_dq_n = grad_L_q(q_mid_next, qdot_next)
        dL_dqdot_n = grad_L_qdot(q_mid_next, qdot_next)
        F = dL_dqdot + dL_dqdot_n - h/2 * (dL_dq_n - dL_dq)
        if norm(F) < tol; return q_next; end
        M_next = M_func(q_next); dq = (h/2) * (M_next \ F)
        q_next -= dq
    end
    return q_next
end

"""Full variational integrator simulation from (q0, v0)."""
function solve_variational(L::Function, grad_L_q::Function, grad_L_qdot::Function, M_func::Function,
                           q0::Vector{Float64}, v0::Vector{Float64}, t_end::Float64, h::Float64)
    M0 = M_func(q0); a0 = M0 \ (grad_L_q(q0, v0) - zeros(length(q0)))
    q_prev = q0 - h*v0 + 0.5*h^2*a0; q_cur = copy(q0)
    results = [(0.0, q0, v0)]; t = 0.0
    while t < t_end
        q_next = variational_midpoint_step(L, grad_L_q, grad_L_qdot, M_func, q_prev, q_cur, h)
        v_next = (q_next - q_prev) / (2*h)
        q_prev, q_cur = q_cur, q_next; t += h
        push!(results, (t, q_next, v_next))
    end
    return results
end

"""Verify least action: compare true path action vs varied path action."""
function verify_least_action(L::Function, q_true::Function, qdot_true::Function, q_varied::Function, qdot_varied::Function, t_span::Tuple{Float64,Float64})
    S_true = action(L, q_true, qdot_true, t_span)
    S_varied = action(L, q_varied, qdot_varied, t_span)
    return S_true, S_varied, S_varied - S_true
end

"""First variation delta S for perturbation eta(t) with eta(t0)=eta(tf)=0."""
function first_variation(L::Function, grad_L_q::Function, grad_L_qdot::Function,
                         q_func::Function, qdot_func::Function, eta::Function, etadot::Function,
                         t_span::Tuple{Float64,Float64}, n_points=1000)
    t0,tf=t_span; dt=(tf-t0)/(n_points-1); dS=0.0
    for i in 1:n_points
        t=t0+(i-1)*dt; w=(i==1||i==n_points)?0.5:1.0
        dS += w*(dot(grad_L_q(q_func(t),qdot_func(t)),eta(t)) + dot(grad_L_qdot(q_func(t),qdot_func(t)),etadot(t)))*dt
    end
    return dS
end

"""Second variation sign: positive definite=minimum, indefinite=saddle."""
function second_variation_sign(L::Function, grad_L_q::Function, grad_L_qdot::Function,
                               q_func::Function, qdot_func::Function, eta::Function, etadot::Function,
                               t_span::Tuple{Float64,Float64}; n_points=500, eps_val=1e-5)
    t0,tf=t_span; dt=(tf-t0)/(n_points-1); d2S=0.0
    for i in 1:n_points
        t=t0+(i-1)*dt; w=(i==1||i==n_points)?0.5:1.0
        q=q_func(t); qd=qdot_func(t); n_q=length(q)
        H_qq=zeros(n_q,n_q); H_qdqd=zeros(n_q,n_q); H_qqd=zeros(n_q,n_q)
        for a in 1:n_q
            qp=copy(q);qp[a]+=eps_val;qm=copy(q);qm[a]-=eps_val
            H_qq[:,a]=(grad_L_q(qp,qd)-grad_L_q(qm,qd))/(2eps_val)
            H_qqd[:,a]=(grad_L_qdot(qp,qd)-grad_L_qdot(qm,qd))/(2eps_val)
            qdp=copy(qd);qdp[a]+=eps_val;qdm=copy(qd);qdm[a]-=eps_val
            H_qdqd[:,a]=(grad_L_qdot(q,qdp)-grad_L_qdot(q,qdm))/(2eps_val)
        end
        e=eta(t);ed=etadot(t)
        d2S+=w*(dot(e,H_qq*e)+2*dot(e,H_qqd*ed)+dot(ed,H_qdqd*ed))*dt
    end
    return d2S>0 ? :minimum : :saddle
end

"""Symplectic midpoint integrator for separable L=T-V (simplified, no implicit solve)."""
function symplectic_midpoint_step(grad_V::Function, M::Matrix{Float64}, q::Vector{Float64}, p::Vector{Float64}, h::Float64)
    p_half = p - 0.5*h*grad_V(q)
    q_new = q + h*(M \ p_half)
    p_new = p_half - 0.5*h*grad_V(q_new)
    return q_new, p_new
end

"""Action for discrete trajectory (list of q values)."""
function discrete_action(L_d::Function, qs::Vector{Vector{Float64}}, h::Float64)
    S = 0.0
    for k in 1:length(qs)-1
        S += L_d(qs[k], qs[k+1], h)
    end
    return S
end

"""Linearized action variation: δ²S for assessing stability."""
function linearized_action_hessian(L::Function, q0::Vector{Float64}, qdot0::Vector{Float64}, h::Float64; eps_val=1e-5)
    n = length(q0)
    H = zeros(n, n)
    L0 = L(q0, qdot0)
    for i in 1:n
        for j in 1:n
            qpp = copy(q0); qpp[i] += eps_val; qpp[j] += eps_val
            qpm = copy(q0); qpm[i] += eps_val; qpm[j] -= eps_val
            qmp = copy(q0); qmp[i] -= eps_val; qmp[j] += eps_val
            qmm = copy(q0); qmm[i] -= eps_val; qmm[j] -= eps_val
            H[i,j] = (L(qpp,qdot0) - L(qpm,qdot0) - L(qmp,qdot0) + L(qmm,qdot0))/(4*eps_val^2)
        end
    end
    return H
end

## ============================================================
## 高阶变分积分器
## ============================================================

"""
四阶辛分块 Runge-Kutta (Gauss-Legendre 配点法).
用于扩展中点法到四阶精度. 基于 collocation at Gauss points c = 1/2 +/- sqrt(3)/6.
"""
function gauss_legendre_4th_order_step(L::Function, dL_dq::Function, dL_dqdot::Function, M::Function, q0::Vector{Float64}, p0::Vector{Float64}, h::Float64; tol=1e-12, max_iter=30)
    n = length(q0)
    a11 = 1/4; a12 = 1/4 - sqrt(3)/6
    a21 = 1/4 + sqrt(3)/6; a22 = 1/4
    b1 = 1/2; b2 = 1/2
    c1 = 1/2 - sqrt(3)/6; c2 = 1/2 + sqrt(3)/6
    Q1 = copy(q0); Q2 = copy(q0)
    P1 = copy(p0); P2 = copy(p0)
    for iter in 1:max_iter
        q1 = q0 + h*(a11*Q1 + a12*Q2)
        q2 = q0 + h*(a21*Q1 + a22*Q2)
        p1 = p0 + h*(a11*P1 + a12*P2)
        p2 = p0 + h*(a21*P1 + a22*P2)
        F1 = P1 - dL_dqdot(q1, Q1)
        F2 = P2 - dL_dqdot(q2, Q2)
        G1 = M(q1)*Q1 - p1
        G2 = M(q2)*Q2 - p2
        if norm([F1;F2;G1;G2]) < tol
            break
        end
        Q1 -= 0.5 * (M(q1) \ F1); Q2 -= 0.5 * (M(q2) \ F2)
        P1 -= 0.5 * G1; P2 -= 0.5 * G2
    end
    q_new = q0 + h*(b1*Q1 + b2*Q2)
    p_new = p0 + h*(b1*P1 + b2*P2)
    return q_new, p_new
end
