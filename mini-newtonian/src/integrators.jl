# integrators.jl — ODE 积分器（手写，零外部依赖）
# 涵盖：Euler, RK2, RK4, Velocity Verlet, Leapfrog, 自适应 RK45

## ============================================================
## 一阶 ODE 积分器：dy/dt = f(t, y)
## ============================================================

"""
欧拉法（一阶精度）
"""
function euler_step(f::Function, t::Float64, y::Vector{Float64}, dt::Float64)
    return y + dt * f(t, y)
end

"""
中点法 / RK2（二阶精度）
"""
function rk2_step(f::Function, t::Float64, y::Vector{Float64}, dt::Float64)
    k1 = f(t, y)
    k2 = f(t + dt/2, y + (dt/2) * k1)
    return y + dt * k2
end

"""
经典 RK4（四阶精度）— 物理模拟的主力
"""
function rk4_step(f::Function, t::Float64, y::Vector{Float64}, dt::Float64)
    half_dt = dt / 2
    k1 = f(t, y)
    k2 = f(t + half_dt, y + half_dt * k1)
    k3 = f(t + half_dt, y + half_dt * k2)
    k4 = f(t + dt, y + dt * k3)
    return y + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
end

"""
RK4 一步，同时返回 k1..k4 用于后续误差估计
"""
function rk4_step_with_stages(f::Function, t::Float64, y::Vector{Float64}, dt::Float64)
    half_dt = dt / 2
    k1 = f(t, y)
    k2 = f(t + half_dt, y + half_dt * k1)
    k3 = f(t + half_dt, y + half_dt * k2)
    k4 = f(t + dt, y + dt * k3)
    y_new = y + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
    return y_new, (k1, k2, k3, k4)
end

"""
自适应 RK45（Dormand-Prince 5(4) 嵌入对）
容错 tol，返回 (y_new, err_estimate, suggested_dt)
"""
function rk45_step(f::Function, t::Float64, y::Vector{Float64}, dt::Float64, tol::Float64=1e-8)
    # Dormand-Prince 系数: 7 个阶段，5阶解 + 4阶嵌入
    c2, c3, c4, c5, c6, c7 = 1/5, 3/10, 4/5, 8/9, 1.0, 1.0

    k1 = f(t, y)
    k2 = f(t + c2*dt, y + dt*( (1/5) * k1 ))
    k3 = f(t + c3*dt, y + dt*( (3/40)*k1 + (9/40)*k2 ))
    k4 = f(t + c4*dt, y + dt*( (44/45)*k1 + (-56/15)*k2 + (32/9)*k3 ))
    k5 = f(t + c5*dt, y + dt*( (19372/6561)*k1 + (-25360/2187)*k2 + (64448/6561)*k3 + (-212/729)*k4 ))
    k6 = f(t + c6*dt, y + dt*( (9017/3168)*k1 + (-355/33)*k2 + (46732/5247)*k3 + (49/176)*k4 + (-5103/18656)*k5 ))
    k7 = f(t + c7*dt, y + dt*( (35/384)*k1 + (500/1113)*k3 + (125/192)*k4 + (-2187/6784)*k5 + (11/84)*k6 ))

    # 5阶解
    y5 = y + dt*( (35/384)*k1 + (500/1113)*k3 + (125/192)*k4 + (-2187/6784)*k5 + (11/84)*k6 )
    # 4阶嵌入解
    y4 = y + dt*( (5179/57600)*k1 + (7571/16695)*k3 + (393/640)*k4 + (-92097/339200)*k5 + (187/2100)*k6 + (1/40)*k7 )

    # 误差估计（scaled RMS norm）
    len_y = length(y)
    err_sq = 0.0
    for i in 1:len_y
        scale = max(1.0, max(abs(y5[i]), abs(y[i])))
        err_sq += ((y5[i] - y4[i]) / scale)^2
    end
    err = sqrt(err_sq / len_y)

    # 步长控制（PI 控制器）
    safety = 0.9
    dt_new = err > 0 ? dt * safety * (tol / err)^0.2 : dt * 2.0

    return y5, err, min(dt_new, dt * 4.0)
end

## ============================================================
## 二阶 ODE 专用积分器: d^2r/dt^2 = a(r, v, t)
## ============================================================

"""
Euler-Cromer（半隐式欧拉）
对保守系统能量漂移远小于显式 Euler
"""
function euler_cromer_step(a_func::Function, state::State, dt::Float64)
    a_cur = a_func(state.r, state.v, state.t)
    v_new = state.v + a_cur * dt
    r_new = state.r + v_new * dt  # 用新速度！
    return State(state.t + dt, r_new, v_new)
end

"""
速度 Verlet — 辛积分器，守恒系统首选
该算法对称、时间可逆、长时能量守恒优异
"""
function velocity_verlet_step(a_func::Function, state::State, dt::Float64)
    a_cur = a_func(state.r, state.v, state.t)
    v_half = state.v + 0.5 * a_cur * dt
    r_new = state.r + v_half * dt
    a_new = a_func(r_new, v_half, state.t + dt)
    v_new = v_half + 0.5 * a_new * dt
    return State(state.t + dt, r_new, v_new)
end

"""
Leapfrog 蛙跳法初始化
返回 (r0, v_half_0)
"""
function leapfrog_init(a_func::Function, r::Vec3, v::Vec3, t::Float64, dt::Float64)
    a0 = a_func(r, v, t)
    v_half = v + 0.5 * a0 * dt
    return r, v_half
end

function leapfrog_step(a_func::Function, r::Vec3, v_half::Vec3, t::Float64, dt::Float64)
    r_new = r + v_half * dt
    a_new = a_func(r_new, v_half, t + dt)
    v_new_half = v_half + a_new * dt
    return r_new, v_new_half
end

"""
用 RK4 解二阶 ODE（转为两个一阶方程耦合）
比 Verlet 通用（可处理速度依赖力），但无辛性质
"""
function rk4_second_order_step(a_func::Function, state::State, dt::Float64)
    function f(t, y)
        r = Vec3(y[1], y[2], y[3])
        v = Vec3(y[4], y[5], y[6])
        a = a_func(r, v, t)
        return [v.x, v.y, v.z, a.x, a.y, a.z]
    end
    y0 = [state.r.x, state.r.y, state.r.z, state.v.x, state.v.y, state.v.z]
    y1 = rk4_step(f, state.t, y0, dt)
    return State(state.t + dt, Vec3(y1[1],y1[2],y1[3]), Vec3(y1[4],y1[5],y1[6]))
end

## ============================================================
## 通用求解循环
## ============================================================

"""
固定步长积分（粒子问题）
step_func: (a_func, state, dt) -> new_state
"""
function solve_fixed_step(a_func::Function, r0::Vec3, v0::Vec3, t_end::Float64, dt::Float64;
                          method::Symbol=:verlet, record_every::Int=1, callback=nothing)
    if method == :verlet
        step_func = velocity_verlet_step
    elseif method == :euler_cromer
        step_func = euler_cromer_step
    elseif method == :rk4
        step_func = rk4_second_order_step
    else
        error("Unknown method: $method, use :verlet, :euler_cromer, or :rk4")
    end

    state = State(r0, v0)
    traj = Trajectory()
    record!(traj, state)

    step = 0
    while state.t < t_end
        actual_dt = min(dt, t_end - state.t)
        state = step_func(a_func, state, actual_dt)
        step += 1
        if step % record_every == 0
            record!(traj, state)
        end
        if callback !== nothing
            callback(state)
        end
    end

    if traj.ts[end] < t_end
        record!(traj, state)
    end

    return traj
end

"""
自适应步长积分（粒子问题）
使用 RK45 自动调整步长
"""
function solve_adaptive(a_func::Function, r0::Vec3, v0::Vec3, t_end::Float64;
                        dt_init::Float64=0.01, tol::Float64=1e-8,
                        dt_min::Float64=1e-12, dt_max::Float64=0.1,
                        record_every::Int=1)

    function f(t, y)
        r = Vec3(y[1], y[2], y[3])
        v = Vec3(y[4], y[5], y[6])
        a = a_func(r, v, t)
        return [v.x, v.y, v.z, a.x, a.y, a.z]
    end

    y = [r0.x, r0.y, r0.z, v0.x, v0.y, v0.z]
    t = 0.0
    dt = dt_init
    traj = Trajectory()
    record!(traj, t, Vec3(y[1],y[2],y[3]), Vec3(y[4],y[5],y[6]))

    n_accepted = 0
    n_rejected = 0

    while t < t_end
        if t + dt > t_end
            dt = t_end - t
        end

        y_new, err, dt_suggested = rk45_step(f, t, y, dt, tol)

        if err <= tol
            y = y_new
            t = t + dt
            n_accepted += 1
            if n_accepted % record_every == 0
                record!(traj, t, Vec3(y[1],y[2],y[3]), Vec3(y[4],y[5],y[6]))
            end
            dt = clamp(dt_suggested, dt_min, dt_max)
        else
            n_rejected += 1
            dt = clamp(dt_suggested, dt_min, dt)
        end

        if dt < dt_min
            @warn "Step size too small, integration may be unstable" t dt err
            break
        end
    end

    if traj.ts[end] < t_end
        record!(traj, t, Vec3(y[1],y[2],y[3]), Vec3(y[4],y[5],y[6]))
    end

    return traj, n_accepted, n_rejected
end

## ============================================================
## N 体积分
## ============================================================

"""
N 体系统 Velocity Verlet 一步
a_func(st::NBodyState) -> Vector{Vec3} 返回每个粒子的加速度
"""
function nbody_verlet_step(a_func::Function, st::NBodyState, dt::Float64)
    n = n_particles(st)
    acc_cur = a_func(st)

    v_half = [st.velocities[i] + 0.5 * acc_cur[i] * dt for i in 1:n]
    r_new  = [st.positions[i] + v_half[i] * dt for i in 1:n]

    st_new = NBodyState(st.masses, r_new, v_half, t=st.t + dt)
    acc_new = a_func(st_new)

    v_new = [v_half[i] + 0.5 * acc_new[i] * dt for i in 1:n]

    return NBodyState(st.masses, r_new, v_new, t=st.t + dt)
end
