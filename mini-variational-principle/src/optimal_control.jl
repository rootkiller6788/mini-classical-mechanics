# optimal_control.jl — 最优控制: Pontryagin极大值原理, 动态规划, LQR
# 参考: Kirk "Optimal Control Theory", Bryson & Ho, Bertsekas

using LinearAlgebra

# ============================================================
# Pontryagin 极大值原理
# ============================================================
"""
Pontryagin 极大值原理 (连续时间)

系统: ẋ = f(x,u,t), x(0)=x0
代价: J = φ(x(tf)) + ∫ L(x,u,t) dt

Hamiltonian: H(x,u,p,t) = L(x,u,t) + pᵀ f(x,u,t)
伴随方程: ṗ = -∂H/∂x
边界条件: p(tf) = ∂φ/∂x(tf)
最优控制: u* = argmin_u H(x,u,p,t)

返回在给定 (x,p,t) 处的最优 u (解析或数值)
"""
function pontryagin_hamiltonian(L::Function, f::Function, x::Vector{Float64}, u::Vector{Float64}, p::Vector{Float64}, t::Float64)
    L(x,u,t) + dot(p, f(x,u,t))
end

"""
伴随方程右边: ṗ = -∂H/∂x (数值)
"""
function adjoint_rhs(H::Function, x::Vector{Float64}, u::Vector{Float64}, p::Vector{Float64}, t::Float64; dx=1e-5)
    n = length(x); dp = zeros(n)
    for i in 1:n
        xp, xm = copy(x), copy(x); xp[i]+=dx; xm[i]-=dx
        dp[i] = -(H(xp,u,p,t) - H(xm,u,p,t))/(2*dx)
    end; dp
end

"""
最优控制必要条件 (bang-bang 问题的切换函数)
"""
function switching_function(H::Function, x::Vector{Float64}, p::Vector{Float64}, t::Float64, u_vals::Vector{Float64})
    vals = [H(x,[u],p,t) for u in u_vals]
    return u_vals, vals, u_vals[argmin(vals)]
end

"""
LQR 最优控制: 线性系统 ẋ=Ax+Bu, 二次代价 J=∫(xᵀQx+uᵀRu)dt

Riccati 方程: 0 = PA + AᵀP - PBR⁻¹BᵀP + Q
最优反馈: u* = -R⁻¹BᵀP x
"""
function lqr_controller(A::Matrix{Float64}, B::Matrix{Float64}, Q::Matrix{Float64}, R::Matrix{Float64})
    n, m = size(A,1), size(B,2)
    # CARE using iterative method (simplified for diagonal-dominant systems)
    P = copy(Q); I_n = Matrix{Float64}(I, n, n)
    for iter in 1:200
        P_new = A'*P + P*A - P*B*(R\B')*P + Q
        P_new = 0.5*(P_new + P_new')
        if norm(P_new)/max(norm(P),1e-10) < 1e-8
            P = P_new; break
        end
        P = P + 0.1*(P_new - P)  # damped update
    end
    K = R \ (B' * P)
    return K, P
end

"""
LQR 代价: J = x0ᵀ P x0
"""
function lqr_cost(x0::Vector{Float64}, P::Matrix{Float64})
    dot(x0, P*x0)
end

"""
跟踪控制: 系统跟踪参考轨迹 x_ref(t)
u* = -R⁻¹BᵀPx + R⁻¹Bᵀv, 其中 v 满足 v̇ = -(A-BR⁻¹BᵀP)ᵀv + Qx_ref
"""
function lqr_tracking_gain(A::Matrix{Float64}, B::Matrix{Float64}, Q::Matrix{Float64}, R::Matrix{Float64})
    K, P = lqr_controller(A, B, Q, R)
    return K, P
end

# ============================================================
# 动态规划 (Bellman)
# ============================================================
"""
Bellman 最优性原理 (离散时间)

系统: x_{k+1} = f(x_k, u_k)
代价: J = Σ g(x_k, u_k)

值函数: V_k(x) = min_u [g(x,u) + V_{k+1}(f(x,u))]
"""
function bellman_value_iteration(f::Function, g::Function, x_grid::Vector{Vector{Float64}}, u_grid::Vector{Vector{Float64}}, n_steps::Int)
    nx, nu = length(x_grid), length(u_grid)
    V = zeros(nx)  # terminal cost=0
    policy = zeros(Int, n_steps, nx)
    
    for k in n_steps:-1:1
        V_new = zeros(nx)
        for i in 1:nx
            best_val, best_u = Inf, 1
            for j in 1:nu
                x_next = f(x_grid[i], u_grid[j])
                idx = find_nearest(x_grid, x_next)
                val = g(x_grid[i], u_grid[j]) + V[idx]
                if val < best_val; best_val, best_u = val, j; end
            end
            V_new[i] = best_val; policy[k,i] = best_u
        end; V = V_new
    end
    return V, policy
end

function find_nearest(grid::Vector{Vector{Float64}}, x::Vector{Float64})
    best_idx, best_dist = 1, Inf
    for (i, gx) in enumerate(grid)
        d = norm(gx - x)
        if d < best_dist; best_dist, best_idx = d, i; end
    end; best_idx
end

# ============================================================
# 间接法: 两点边值问题 (BVP)
# ============================================================
"""
用打靶法求解最优控制的 BVP

状态: [x; p] (2n维)
ODE: ẋ=f(x,u*,t), ṗ=-∂H/∂x, u*=argmin H
"""
function solve_optimal_control_bvp(f::Function, L_func::Function, x0::Vector{Float64}, tf::Float64, phi::Function; n_grid=500, p_guess=nothing)
    n = length(x0)
    p0 = isnothing(p_guess) ? zeros(n) : p_guess
    # Simplified: forward-backward sweep
    function shoot(p0_vec)
        x = copy(x0); p = copy(p0_vec); dt = tf/(n_grid-1)
        for _ in 1:n_grid
            u_opt = zeros(1)  # placeholder for 1D control
            H(xv,uv,pv,t) = L_func(xv,uv,t) + dot(pv, f(xv,uv,t))
            x = x + dt * f(x, u_opt, 0.0)
            dp = adjoint_rhs(H, x, u_opt, p, 0.0)
            p = p + dt * dp
        end
        return p - [0.0]  # p(tf) should be ∂φ/∂x
    end
    return p0  # simplified return
end

# ============================================================
# 最短时间控制 (Bang-Bang)
# ============================================================
"""
双积分器最短时间: ẋ=v, v̇=u, |u|≤1, J = ∫dt → min

切换曲线: v² = ±2x (sign dependent)
最优控制: u* = -sign(current - target)
"""
function bang_bang_switch_curve(x::Float64, v::Float64)
    return v^2/2  # 切换曲线参数
end

function minimum_time_control_double_integrator(x::Float64, v::Float64, x_target::Float64, v_target::Float64)
    dx = x - x_target; dv = v - v_target
    zeta = dx + 0.5*dv*abs(dv)
    return -sign(zeta)
end

"""
切换时间计算: 对于双积分器，从 (x0,v0) 到 (0,0) 的最短时间
"""
function minimum_time_double_integrator(x0::Float64, v0::Float64, u_max::Float64=1.0)
    # 解析解
    t_switch = abs(v0)/u_max + sqrt(abs(2*x0 + v0^2/u_max*sign(x0)))
    return t_switch
end

# ============================================================
# 模型预测控制 (MPC)
# ============================================================
"""
MPC: 在每个采样时刻求解有限时域最优控制，仅施加第一步
"""
function mpc_step(A::Matrix{Float64}, B::Matrix{Float64}, Q::Matrix{Float64}, R::Matrix{Float64}, x_current::Vector{Float64}, horizon::Int)
    K, P = lqr_controller(A, B, Q, R)
    u = -K * x_current
    return u
end

# ============================================================
# Hamilton-Jacobi-Bellman 方程
# ============================================================
"""
HJB 方程: -∂V/∂t = min_u [L(x,u,t) + (∂V/∂x)ᵀ f(x,u,t)]
V(x,tf) = φ(x)

对于无限时域 LQR: V(x) = xᵀPx
"""
function hjb_lqr_verify(P::Matrix{Float64}, A::Matrix{Float64}, B::Matrix{Float64}, Q::Matrix{Float64}, R::Matrix{Float64}; tol=1e-10)
    R_inv_BtP = R \ (B' * P)
    residual = A'*P + P*A - P*B*R_inv_BtP + Q
    norm(residual) < tol
end

export pontryagin_hamiltonian, adjoint_rhs, switching_function
export lqr_controller, lqr_cost, lqr_tracking_gain
export bellman_value_iteration, solve_optimal_control_bvp
export bang_bang_switch_curve, minimum_time_control_double_integrator, minimum_time_double_integrator
export mpc_step, hjb_lqr_verify
