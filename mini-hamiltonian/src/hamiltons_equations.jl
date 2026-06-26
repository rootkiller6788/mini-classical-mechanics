# hamiltons_equations.jl — Hamilton 方程求解器
# dq/dt = ∂H/∂p,  dp/dt = -∂H/∂q

using LinearAlgebra

struct HamiltonianSystem
    n::Int
    H::Function              # H(q, p) -> Float64
    grad_H::Function         # (q, p) -> (∂H/∂q, ∂H/∂p)
end

"""
从解析梯度构建 Hamilton 系统（无数值误差）
"""
function analytic_hamiltonian_system(H::Function, grad_q_H::Function, grad_p_H::Function, n::Int)
    function grad_H(q, p)
        return grad_q_H(q, p), grad_p_H(q, p)
    end
    return HamiltonianSystem(n, H, grad_H)
end

"""
从数值梯度构建 Hamilton 系统（通用，适用于任意 H）
"""
function numerical_hamiltonian_system(H::Function, n::Int)
    function grad_H(q, p)
        dH_dq = zeros(n); dH_dp = zeros(n)
        eps_val = 1e-6
        for i in 1:n
            qp = copy(q); qp[i] += eps_val; qm = copy(q); qm[i] -= eps_val
            dH_dq[i] = (H(qp, p) - H(qm, p)) / (2*eps_val)
            pp = copy(p); pp[i] += eps_val; pm = copy(p); pm[i] -= eps_val
            dH_dp[i] = (H(q, pp) - H(q, pm)) / (2*eps_val)
        end
        return dH_dq, dH_dp
    end
    return HamiltonianSystem(n, H, grad_H)
end

"""
Hamilton 方程 → 一阶 ODE: y = [q; p]
"""
function hamiltons_ode(sys::HamiltonianSystem)
    function f(t, y)
        n = sys.n; q = y[1:n]; p = y[n+1:2n]
        dH_dq, dH_dp = sys.grad_H(q, p)
        return [dH_dp; -dH_dq]
    end
    return f
end

"""
Symplectic Euler (1st order symplectic)
q_{n+1} = q_n + h*∂H/∂p(q_{n+1}, p_n)  ← 隐式
p_{n+1} = p_n - h*∂H/∂q(q_{n+1}, p_n)

简化：显式 symplectic Euler (q先用旧值)
q_{n+1} = q_n + h*∂H/∂p(q_n, p_n)
p_{n+1} = p_n - h*∂H/∂q(q_{n+1}, p_n)   ← 用新q
"""
function symplectic_euler_step(sys::HamiltonianSystem, q::Vector{Float64}, p::Vector{Float64}, dt::Float64)
    _, dH_dp_old = sys.grad_H(q, p)
    q_new = q + dt * dH_dp_old
    dH_dq_new, _ = sys.grad_H(q_new, p)
    p_new = p - dt * dH_dq_new
    return q_new, p_new
end

"""
Stormer-Verlet (2nd order symplectic, from Hamiltonian perspective)
与 mini-newtonian 的 velocity_verlet_step 等价
"""
function stormer_verlet_step(sys::HamiltonianSystem, q::Vector{Float64}, p::Vector{Float64}, dt::Float64)
    # 半步 p
    dH_dq, _ = sys.grad_H(q, p)
    p_half = p - 0.5*dt*dH_dq
    # 全步 q
    _, dH_dp = sys.grad_H(q, p_half)
    q_new = q + dt*dH_dp
    # 完成 p
    dH_dq_new, _ = sys.grad_H(q_new, p_half)
    p_new = p_half - 0.5*dt*dH_dq_new
    return q_new, p_new
end

"""
固定步长求解器
"""
function solve_hamiltonian(sys::HamiltonianSystem, q0::Vector{Float64}, p0::Vector{Float64}, t_end::Float64, dt::Float64; method=:verlet, record_every=1)
    q, p = copy(q0), copy(p0); t = 0.0
    traj = PhaseTrajectory()
    H0 = sys.H(q, p)
    record!(traj, t, q, p, H0)
    step = 0
    while t < t_end
        h = min(dt, t_end - t)
        if method == :verlet
            q, p = stormer_verlet_step(sys, q, p, h)
        elseif method == :symplectic_euler
            q, p = symplectic_euler_step(sys, q, p, h)
        else
            f = hamiltons_ode(sys); y = [q; p]
            y = rk4_step(f, t, y, h); q, p = y[1:sys.n], y[sys.n+1:end]
        end
        t += h; step += 1
        if step % record_every == 0
            record!(traj, t, q, p, sys.H(q, p))
        end
    end
    if traj.ts[end] < t_end
        record!(traj, t, q, p, sys.H(q, p))
    end
    return traj
end

function rk4_step(f, t, y, dt)
    half = dt/2; k1=f(t,y); k2=f(t+half,y+half*k1); k3=f(t+half,y+half*k2); k4=f(t+dt,y+dt*k3)
    return y + (dt/6)*(k1+2k2+2k3+k4)
end
