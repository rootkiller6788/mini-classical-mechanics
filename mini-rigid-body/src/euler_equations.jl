# euler_equations.jl -- Euler 动力学方程与刚体运动数值积分
# 参考: Goldstein Ch.5, Landau Ch.6

using LinearAlgebra

# ============================================================
# Euler 方程（体坐标系，无力矩情况）
# ============================================================

"""
Euler 方程（无力矩情况）

参数:
- I: 惯性张量 (InertiaTensor)
- w: 角速度向量 [ωx, ωy, ωz]（体坐标系分量）

返回: dω/dt = [dωx/dt, dωy/dt, dωz/dt]

公式 (Goldstein 5.44):
  I1 ω̇₁ = (I2 - I3) ω₂ ω₃
  I2 ω̇₂ = (I3 - I1) ω₃ ω₁
  I3 ω̇₃ = (I1 - I2) ω₁ ω₂
"""
function euler_equations_derivative(I::InertiaTensor, w::AbstractVector)
    I1, I2, I3 = principal_moments(I)
    return [
        (I2 - I3) * w[2] * w[3] / I1,
        (I3 - I1) * w[3] * w[1] / I2,
        (I1 - I2) * w[1] * w[2] / I3
    ]
end

"""
Euler 方程（含外力矩）

参数:
- I: 惯性张量
- w: 角速度（体坐标系）
- N: 外力矩（体坐标系分量）

公式:
  I1 ω̇₁ - (I2 - I3) ω₂ ω₃ = N1
  ...
→ ω̇₁ = [(I2 - I3) ω₂ ω₃ + N1] / I1
"""
function euler_with_torque(I::InertiaTensor, w::AbstractVector, N::AbstractVector)
    I1, I2, I3 = principal_moments(I)
    dw = euler_equations_derivative(I, w)
    return [dw[1] + N[1]/I1, dw[2] + N[2]/I2, dw[3] + N[3]/I3]
end

"""
角加速度（含外力矩，解耦近似: 适用于 I 已对角化）
"""
function angular_acceleration(I::InertiaTensor, w::AbstractVector, N::AbstractVector)
    return euler_with_torque(I, w, N)
end

# ============================================================
# 无力矩刚体的运动常数
# ============================================================

"""
无力矩 Euler 方程的两个运动常数 (Goldstein 5.45):
- T = ½ (I1 ω₁² + I2 ω₂² + I3 ω₃²) = const  (动能守恒)
- L² = I1² ω₁² + I2² ω₂² + I3² ω₃² = const  (角动量守恒)

返回: (T, L²)
"""
function motion_constants(I::InertiaTensor, w::AbstractVector)
    I1, I2, I3 = principal_moments(I)
    T = 0.5 * (I1*w[1]^2 + I2*w[2]^2 + I3*w[3]^2)
    L2 = (I1*w[1])^2 + (I2*w[2])^2 + (I3*w[3])^2
    return (T, L2)
end

"""
用运动常数 T 和 L² 描述 Binet 椭球交线
"""
function polhode_check(I::InertiaTensor, w::AbstractVector; tol::Float64=1e-10)
    T, L2 = motion_constants(I, w)
    I1, I2, I3 = principal_moments(I)
    # 验证角速度同时在动能椭球和角动量球面上
    T_check = 0.5 * (I1*w[1]^2 + I2*w[2]^2 + I3*w[3]^2)
    L2_check = (I1*w[1])^2 + (I2*w[2])^2 + (I3*w[3])^2
    return abs(T - T_check) < tol && abs(L2 - L2_check) < tol
end

# ============================================================
# 积分器（体坐标系角速度演化）
# ============================================================

"""
Euler 步（一阶显式）
"""
function euler_step_omega(I::InertiaTensor, w::Vector{Float64}, dt::Float64)
    dw = euler_equations_derivative(I, w)
    return w + dw * dt
end

"""
RK4 步（四阶精度）— 保持运动常数精度
"""
function rk4_step_omega(I::InertiaTensor, w::Vector{Float64}, dt::Float64)
    half_dt = dt / 2.0
    k1 = euler_equations_derivative(I, w)
    k2 = euler_equations_derivative(I, w + half_dt * k1)
    k3 = euler_equations_derivative(I, w + half_dt * k2)
    k4 = euler_equations_derivative(I, w + dt * k3)
    return w + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
end

"""
Velocity Verlet 适用于刚体角速度演化

半步加速 → 全步位置 → 新加速 → 半步加速
"""
function verlet_step_omega(I::InertiaTensor, w::Vector{Float64}, dt::Float64)
    a_cur = euler_equations_derivative(I, w)
    w_half = w + 0.5 * a_cur * dt
    a_new = euler_equations_derivative(I, w_half)
    return w_half + 0.5 * a_new * dt
end

"""
模拟无力矩刚体运动（体坐标系角速度演化）

参数:
- I: 惯性张量
- w0: 初始角速度 [ωx, ωy, ωz]
- t_end: 模拟时长
- dt: 时间步长
- method: :euler, :rk4, :verlet

返回: trajectory::Vector{Vector{Float64}} 每步的 ω
"""
function simulate_free_rigid_body(I::InertiaTensor, w0::AbstractVector, t_end::Float64, dt::Float64;
                                  method::Symbol=:rk4)
    if method == :euler
        step_func = euler_step_omega
    elseif method == :rk4
        step_func = rk4_step_omega
    elseif method == :verlet
        step_func = verlet_step_omega
    else
        error("Unknown method: $method. Use :euler, :rk4, or :verlet")
    end

    w = Vector{Float64}(w0)
    traj = [copy(w)]
    n_steps = Int(ceil(t_end / dt))

    for _ in 1:n_steps
        actual_dt = min(dt, t_end - (length(traj)-1)*dt)
        w = step_func(I, w, actual_dt)
        push!(traj, copy(w))
    end

    return traj
end

"""
模拟含外力矩的刚体运动

参数:
- I: 惯性张量
- w0: 初始角速度
- torque_func: 力矩函数 N(t, w) -> Vector{Float64}
- t_end: 模拟时长
- dt: 时间步长

返回: (times, omegas)
"""
function simulate_rigid_body_with_torque(I::InertiaTensor, w0::AbstractVector,
                                          torque_func::Function,
                                          t_end::Float64, dt::Float64)
    w = Vector{Float64}(w0)
    times = [0.0]
    omegas = [copy(w)]
    n_steps = Int(ceil(t_end / dt))

    for step in 1:n_steps
        t = step * dt
        N = torque_func(t, w)
        # RK4 with torque
        f(w, t) = angular_acceleration(I, w, torque_func(t, w))
        half_dt = dt / 2.0

        k1 = f(w, t - dt)
        k2 = f(w + half_dt * k1, t - dt + half_dt)
        k3 = f(w + half_dt * k2, t - dt + half_dt)
        k4 = f(w + dt * k3, t)

        w = w + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        push!(times, t)
        push!(omegas, copy(w))
    end

    return times, omegas
end

# ============================================================
# Poinsot 几何构造
# ============================================================

"""
角速度在体坐标系中的演化——Polhode 曲线点

Poinsot 构造: 无力矩刚体运动 = 惯性椭球在不变平面上的纯滚动
- Polhode: 角速度在体坐标系惯性椭球上画出的闭合曲线
- Herpolhode: 角速度在空间不变平面上画出的曲线

返回 polhode 点: 角速度在体坐标系的分量轨迹
"""
function polhode_curve(I::InertiaTensor, w0::AbstractVector, t_end::Float64, dt::Float64;
                        method::Symbol=:rk4)
    return simulate_free_rigid_body(I, w0, t_end, dt; method=method)
end

"""
体坐标系 ω(t) → 空间坐标系 ω_s(t) = R(t)·ω(t)

参数:
- traj_w_body: 体坐标系 ω 轨迹
- euler_traj: 每步的欧拉角 (若有)
- 或 I 和 w 直接计算
"""
function body_to_space_omega(w_body::AbstractVector, euler::EulerAngles)
    R = euler_to_rotation(euler.phi, euler.theta, euler.psi)
    return R * w_body
end

# ============================================================
# 运动常数监测
# ============================================================

"""
监测模拟过程中的运动常数漂移

返回: (T_history, L2_history)
"""
function monitor_constants(I::InertiaTensor, traj::Vector{Vector{Float64}})
    T_hist = Float64[]
    L2_hist = Float64[]
    for w in traj
        T, L2 = motion_constants(I, w)
        push!(T_hist, T)
        push!(L2_hist, L2)
    end
    return T_hist, L2_hist
end

"""
报告积分器在运动常数上的漂移
"""
function constant_drift_report(I::InertiaTensor, traj::Vector{Vector{Float64}})
    T_hist, L2_hist = monitor_constants(I, traj)
    T0 = T_hist[1]
    L20 = L2_hist[1]
    T_drift = maximum(abs.(T_hist .- T0)) / max(T0, 1e-15)
    L2_drift = maximum(abs.(L2_hist .- L20)) / max(L20, 1e-15)
    return (T_drift=T_drift, L2_drift=L2_drift)
end

export euler_equations_derivative, euler_with_torque, angular_acceleration
export motion_constants, polhode_check
export euler_step_omega, rk4_step_omega, verlet_step_omega
export simulate_free_rigid_body, simulate_rigid_body_with_torque
export polhode_curve, body_to_space_omega
export monitor_constants, constant_drift_report
