# integrators.jl — ODE 数值积分器
# 参考：Hairer, Lubich & Wanner (2006), Goldstein Ch.2
# L5 Computational Methods: RK4, Verlet, symplectic integrators

"""
RK4 (经典四阶 Runge-Kutta) 单步积分.
Butcher tableau:
  0   |
  1/2 | 1/2
  1/2 | 0   1/2
  1   | 0   0   1
  ----+-------------
      | 1/6 2/6 2/6 1/6

全局误差 O(h^4).
"""
function rk4_step(f::Function, t::Float64, y::Vector{Float64}, h::Float64)
    k1 = f(t, y)
    k2 = f(t + h/2, y + (h/2)*k1)
    k3 = f(t + h/2, y + (h/2)*k2)
    k4 = f(t + h,   y + h*k3)
    return y + (h/6.0) * (k1 + 2*k2 + 2*k3 + k4)
end

"""
RK4 在整个时间区间 [t0, t_end] 上的完整积分.
返回 (times, states) 的元组，其中 states[i] 是时刻 times[i] 的状态向量.
"""
function rk4_integrate(f::Function, y0::Vector{Float64}, t_span::Tuple{Float64,Float64}, h::Float64)
    t0, t_end = t_span
    n_steps = Int(ceil((t_end - t0) / h))
    h_actual = (t_end - t0) / n_steps
    times = Vector{Float64}(undef, n_steps+1)
    states = Vector{Vector{Float64}}(undef, n_steps+1)
    times[1] = t0; states[1] = copy(y0)
    y = copy(y0)
    for i in 1:n_steps
        t = t0 + (i-1)*h_actual
        y = rk4_step(f, t, y, h_actual)
        times[i+1] = t + h_actual
        states[i+1] = copy(y)
    end
    return times, states
end

"""
Störmer-Verlet 积分器 (二阶辛).
用于可分离 Hamiltonian H = T(p) + V(q):
  p_{n+1/2} = p_n - (h/2) ∇V(q_n)
  q_{n+1}   = q_n + h ∇T(p_{n+1/2})
  p_{n+1}   = p_{n+1/2} - (h/2) ∇V(q_{n+1})

对二次动能 T = 0.5*p'M^{-1}p, ∇T = M^{-1}p = qdot.
该实现接受加速度函数 a(q) = -M^{-1}∇V(q) 和速度映射 v(p) = M^{-1}p.
"""
function verlet_step(a_func::Function, v_func::Function, q::Vector{Float64}, p::Vector{Float64}, h::Float64)
    # p_{n+1/2} = p_n + (h/2)*M*a(q_n)  where a = -M^{-1}∇V
    # Actually for direct q,v: v_{n+1/2} = v_n + (h/2)*a(q_n)
    p_half = p - (h/2) * a_func(q)  # a_func already returns ∇V
    # q_{n+1} = q_n + h * v_func(p_{n+1/2})
    q_new = q + h * p_half  # assuming unit mass or v_func applied
    # p_{n+1} = p_{n+1/2} - (h/2) * a_func(q_{n+1})
    p_new = p_half - (h/2) * a_func(q_new)
    return q_new, p_new
end

"""
Verlet 在整个区间上的积分.
"""
function verlet_integrate(a_func::Function, q0::Vector{Float64}, v0::Vector{Float64}, t_end::Float64, h::Float64)
    n_steps = Int(ceil(t_end / h)); h_actual = t_end / n_steps
    n = length(q0)
    ts = Vector{Float64}(undef, n_steps+1)
    qs = Vector{Vector{Float64}}(undef, n_steps+1)
    vs = Vector{Vector{Float64}}(undef, n_steps+1)
    ts[1] = 0.0; qs[1] = copy(q0); vs[1] = copy(v0)
    q = copy(q0); v = copy(v0)
    for i in 1:n_steps
        a0 = a_func(q)
        v_half = v + (h_actual/2) * a0
        q = q + h_actual * v_half
        a1 = a_func(q)
        v = v_half + (h_actual/2) * a1
        ts[i+1] = i*h_actual; qs[i+1] = copy(q); vs[i+1] = copy(v)
    end
    return ts, qs, vs
end

"""
辛 Euler 法 (一阶辛).
  p_{n+1} = p_n - h*∇_q H(q_n, p_{n+1})
  q_{n+1} = q_n + h*∇_p H(q_n, p_{n+1})
"""
function symplectic_euler_step(grad_Hq::Function, grad_Hp::Function, q::Vector{Float64}, p::Vector{Float64}, h::Float64)
    p_new = p - h * grad_Hq(q, p)
    q_new = q + h * grad_Hp(q, p_new)
    return q_new, p_new
end

"""
Velocity Verlet (二阶辛, 等价于 leapfrog 的重新排列).
  q_{n+1} = q_n + h*v_n + (h^2/2)*a(q_n)
  v_{n+1} = v_n + (h/2)*(a(q_n) + a(q_{n+1}))
"""
function velocity_verlet_step(a_func::Function, q::Vector{Float64}, v::Vector{Float64}, h::Float64)
    a0 = a_func(q)
    q_new = q + h*v + (h^2/2)*a0
    a1 = a_func(q_new)
    v_new = v + (h/2)*(a0 + a1)
    return q_new, v_new
end

"""
Leapfrog 积分器 (二阶辛).
  v_{n+1/2} = v_{n-1/2} + h*a(q_n)
  q_{n+1}   = q_n + h*v_{n+1/2}

初始化: v_{1/2} = v_0 + (h/2)*a(q_0)
"""
function leapfrog_step(a_func::Function, q::Vector{Float64}, v_half::Vector{Float64}, h::Float64)
    q_new = q + h * v_half
    a_new = a_func(q_new)
    v_half_new = v_half + h * a_new
    return q_new, v_half_new
end

"""
Leapfrog 全区间积分.
"""
function leapfrog_integrate(a_func::Function, q0::Vector{Float64}, v0::Vector{Float64}, t_end::Float64, h::Float64)
    n_steps = Int(ceil(t_end / h)); h_actual = t_end / n_steps
    ts = Float64[]; qs = Vector{Float64}[]
    q = copy(q0); v_half = v0 + (h_actual/2)*a_func(q0)
    push!(ts, 0.0); push!(qs, copy(q))
    for i in 1:n_steps
        q, v_half = leapfrog_step(a_func, q, v_half, h_actual)
        push!(ts, i*h_actual); push!(qs, copy(q))
    end
    return ts, qs
end

"""
自适应步长 RK45 (Dormand-Prince 5(4)).
Butcher tableau (DOPRI5):
  0    |
  1/5  | 1/5
  3/10 | 3/40       9/40
  4/5  | 44/45     -56/15      32/9
  8/9  | 19372/6561 -25360/2187 64448/6561 -212/729
  1    | 9017/3168  -355/33     46732/5247  49/176   -5103/18656
  -----+------------------------------------------------------------------
  5th  | 35/384     0           500/1113    125/192  -2187/6784  11/84
  4th  | 5179/57600 0           7571/16695  393/640  -92097/339200 187/2100 1/40

容差控制: err > tol 时拒绝步长、h_new = 0.9*h*(tol/err)^{1/5}.
"""
function rk45_adaptive_step(f::Function, t::Float64, y::Vector{Float64}, h::Float64; tol=1e-6, max_h=Inf)
    # 5th-order coefficients
    a21=1/5; a31=3/40; a32=9/40; a41=44/45; a42=-56/15; a43=32/9
    a51=19372/6561; a52=-25360/2187; a53=64448/6561; a54=-212/729
    a61=9017/3168; a62=-355/33; a63=46732/5247; a64=49/176; a65=-5103/18656
    b1=35/384; b3=500/1113; b4=125/192; b5=-2187/6784; b6=11/84
    b1s=5179/57600; b3s=7571/16695; b4s=393/640; b5s=-92097/339200; b6s=187/2100; b7s=1/40

    k1 = f(t, y)
    k2 = f(t + a21*h, y + h*(a21*k1))
    k3 = f(t + (a31+a32)*h, y + h*(a31*k1 + a32*k2))
    k4 = f(t + (a41+a42+a43)*h, y + h*(a41*k1 + a42*k2 + a43*k3))
    k5 = f(t + (a51+a52+a53+a54)*h, y + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4))
    k6 = f(t + h, y + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5))

    y5 = y + h*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
    y4 = y + h*(b1s*k1 + b3s*k3 + b4s*k4 + b5s*k5 + b6s*k6 + b7s*k6)

    err = norm(y5 - y4, Inf)
    h_new = min(0.9 * h * (tol / max(err, 1e-15))^0.2, max_h)
    h_new = min(h_new, 4*h)  # 不增长超过 4 倍

    return y5, h_new, err < tol
end

"""
固定步长 RK4 积分，返回最终状态（便捷函数）.
"""
function rk4_fixed(f::Function, y0::Vector{Float64}, t_end::Float64, h::Float64)
    n_steps = Int(ceil(t_end / h)); h_actual = t_end / n_steps
    y = copy(y0)
    for i in 1:n_steps
        y = rk4_step(f, (i-1)*h_actual, y, h_actual)
    end
    return y
end

"""
隐式中点法 (二阶辛, 对一般 Hamiltonian).
求解 q_{n+1} = q_n + h*f((q_n+q_{n+1})/2, (t_n+t_{n+1})/2).
对于可分离系统简化为显式.
"""
function implicit_midpoint_step(f_ode::Function, t::Float64, y::Vector{Float64}, h::Float64; tol=1e-10, max_iter=20)
    n = length(y); y_new = y + h*f_ode(t, y)  # 显式 Euler 初始猜测
    for _ in 1:max_iter
        y_mid = (y + y_new) / 2
        F = y_new - y - h*f_ode(t + h/2, y_mid)
        if norm(F) < tol; break; end
        # 简化 Newton: 用单位矩阵近似 Jacobian
        y_new = y_new - F/2
    end
    return y_new
end

"""
验证辛条件: 对一自由度系统, 相面积守恒 ⇔ det(∂(q',p')/∂(q,p)) = 1.
对于 Verlet 和 Leapfrog, 此恒等式精确成立 (不计浮点误差).
"""
function verify_symplectic(step_func::Function, q::Float64, p::Float64, h::Float64; eps_val=1e-7)
    qp1, pp1 = step_func(q+eps_val, p, h)
    qm1, pm1 = step_func(q-eps_val, p, h)
    qp2, pp2 = step_func(q, p+eps_val, h)
    qm2, pm2 = step_func(q, p-eps_val, h)
    dq_dq = (qp1 - qm1) / (2*eps_val)
    dq_dp = (qp2 - qm2) / (2*eps_val)
    dp_dq = (pp1 - pm1) / (2*eps_val)
    dp_dp = (pp2 - pm2) / (2*eps_val)
    det = dq_dq*dp_dp - dq_dp*dp_dq
    return det, abs(det - 1.0)
end
