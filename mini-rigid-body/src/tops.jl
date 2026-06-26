# tops.jl -- 陀螺运动: 对称陀螺、Lagrange top、进动与章动
# 参考: Goldstein Ch.5.7, Landau Ch.6

using LinearAlgebra

# ============================================================
# 对称陀螺 — 解析理论
# ============================================================

"""
对称陀螺 (I1 = I2 ≠ I3) 的有效势

参数:
- I1: 横向转动惯量 (I1 = I2)
- I3: 轴向转动惯量
- M: 质量
- g: 重力加速度
- l: 质心到支点的距离
- Lz: 空间 z 轴角动量守恒量 (常数)
- L3: 体轴 z 角动量守恒量 (常数 = I3 ω₃)

有效势 (Goldstein 5.58):
  V_eff(θ) = (Lz - L3 cosθ)²/(2 I1 sin²θ) + M g l cosθ

返回: (θ_grid, V_eff)
"""
function symmetric_top_effective_potential(I1::Float64, I3::Float64,
                                            M::Float64, g::Float64, l::Float64,
                                            Lz::Float64, L3::Float64;
                                            n_pts::Int=200)
    thetas = range(1e-4, π-1e-4, length=n_pts)
    V_eff = Float64[]

    for th in thetas
        V1 = (Lz - L3*cos(th))^2 / (2 * I1 * sin(th)^2)
        V2 = M * g * l * cos(th)
        push!(V_eff, V1 + V2)
    end

    return thetas, V_eff
end

"""
对称陀螺的章动角演化

从能量守恒 E = ½I1(θ̇² + φ̇²sin²θ) + ½I3(ψ̇ + φ̇cosθ)² + Mgl cosθ

得到 θ̇ 的方程 (Goldstein 5.55):
  ½ I1 θ̇² = E - Mgl cosθ - (Lz - L3 cosθ)²/(2I1 sin²θ) - L3²/(2I3)

返回 dθ/dt 的符号函数值 (在给定 θ)
"""
function nutation_rate(I1::Float64, I3::Float64, M::Float64, g::Float64, l::Float64,
                        E::Float64, Lz::Float64, L3::Float64, theta::Float64)
    V_eff = (Lz - L3*cos(theta))^2 / (2 * I1 * sin(theta)^2) + M*g*l*cos(theta) + L3^2/(2*I3)
    arg = (E - V_eff) * 2 / I1

    if arg < 0
        return 0.0  # 禁区
    end
    return sqrt(max(arg, 0.0))
end

"""
对称陀螺的进动速率 φ̇(θ)

φ̇ = (Lz - L3 cosθ) / (I1 sin²θ)  (Goldstein 5.57)
"""
function precession_rate_top(I1::Float64, Lz::Float64, L3::Float64, theta::Float64)
    st = sin(theta)
    if abs(st) < 1e-15
        return 0.0
    end
    return (Lz - L3 * cos(theta)) / (I1 * st * st)
end

"""
对称陀螺的自转速率 ψ̇(θ)

ψ̇ = L3/I3 - φ̇ cosθ = (L3/I3) - (Lz - L3 cosθ)cosθ/(I1 sin²θ)
"""
function spin_rate_top(I3::Float64, I1::Float64, Lz::Float64, L3::Float64, theta::Float64)
    phidot = precession_rate_top(I1, Lz, L3, theta)
    return L3 / I3 - phidot * cos(theta)
end

# ============================================================
# 进动分析
# ============================================================

"""
稳态进动条件

条件: θ = const (θ̇ = 0), E = V_eff(θ₀)

返回稳态进动角 θ₀ 的二分法求解

稳态条件 (Goldstein 5.59):
  dV_eff/dθ = 0
"""
function steady_precession_theta(I1::Float64, I3::Float64, M::Float64, g::Float64,
                                  l::Float64, Lz::Float64, L3::Float64;
                                  bracket::Tuple{Float64,Float64}=(0.1, 1.5))
    function dVeff(th)
        st, ct = sin(th), cos(th)
        # V_eff 的导数
        term1 = -(Lz - L3*ct) * (Lz*ct - L3) / (I1 * st^3)
        term2 = -M * g * l * st
        return term1 + term2
    end

    a, b = bracket[1], bracket[2]
    fa, fb = dVeff(a), dVeff(b)

    if fa * fb > 0
        return bracket[1]  # 无根，返回下界
    end

    for _ in 1:50
        mid = (a + b) / 2.0
        fmid = dVeff(mid)
        if abs(fmid) < 1e-12
            return mid
        end
        fa * fmid < 0 ? (b = mid; fb = fmid) : (a = mid; fa = fmid)
    end

    return (a + b) / 2.0
end

"""
快速陀螺近似 (gyroscopic approximation)

当自转角速度 ω₃ 很大时:
  进动速率 ≈ Mgl / (I3 ω₃)
"""
function fast_top_precession(M::Float64, g::Float64, l::Float64,
                              I3::Float64, omega3::Float64)
    return M * g * l / (I3 * omega3)
end

"""
慢速进动近似

当 Lz ≈ L3 cosθ₀:
  φ̇_slow ≈ Mgl / L3
"""
function slow_precession(M::Float64, g::Float64, l::Float64, L3::Float64)
    return M * g * l / L3
end

# ============================================================
# Sleeping Top 分析
# ============================================================

"""
Sleeping Top 稳定性条件

直立陀螺 (θ=0) 在重力场中的稳定性 (Goldstein 5.63):

  ω₃² > 4 M g l I1 / I3²   →  稳定 (继续"睡觉")
  ω₃² < 4 M g l I1 / I3²   →  不稳定 (开始进动+章动)

返回: (:stable, ω₃_critical) 或 (:unstable, ω₃_critical)
"""
function sleeping_top_stability(M::Float64, g::Float64, l::Float64,
                                 I1::Float64, I3::Float64, omega3::Float64)
    omega_crit_sq = 4 * M * g * l * I1 / (I3 * I3)
    omega_crit = sqrt(omega_crit_sq)

    if omega3 > omega_crit
        return (:stable, omega_crit)
    else
        return (:unstable, omega_crit)
    end
end

"""
直立陀螺的临界自转角速度

ω₃_crit = √(4 M g l I1) / I3
"""
function critical_spin_rate(M::Float64, g::Float64, l::Float64,
                             I1::Float64, I3::Float64)
    return sqrt(4 * M * g * l * I1) / I3
end

# ============================================================
# Lagrange Top 完整运动方程
# ============================================================

"""
Lagrange Top 的完整运动方程（一阶 ODE 系统）

状态变量: (θ, φ, ψ, p_θ, p_φ, p_ψ)
其中 p_φ = Lz (空间 z 轴角动量，守恒)
      p_ψ = L3 (体轴 z 角动量，守恒)

返回状态的导数向量

参考: Goldstein (5.54)-(5.57)
"""
function lagrange_top_equations(state::Vector{Float64}, I1::Float64, I3::Float64,
                                 M::Float64, g::Float64, l::Float64)
    theta, phi, psi, ptheta, pphi, ppsi = state

    st = sin(theta); ct = cos(theta)

    # θ̇ = p_θ / I1
    thetadot = ptheta / I1

    # Lz = p_φ (常数), L3 = p_ψ (常数)
    Lz = pphi
    L3 = ppsi

    # φ̇ = (Lz - L3 cosθ) / (I1 sin²θ)
    phidot = abs(st) > 1e-15 ? (Lz - L3*ct) / (I1 * st * st) : 0.0

    # ψ̇ = L3 / I3 - φ̇ cosθ
    psidot = L3 / I3 - phidot * ct

    # ṗ_θ = ∂L/∂θ
    # = (Lz - L3 cosθ)(Lz cosθ - L3) / (I1 sin³θ) - Mgl sinθ
    pthetadot = 0.0
    if abs(st) > 1e-15
        pthetadot = (Lz - L3*ct) * (Lz*ct - L3) / (I1 * st^3) - M*g*l*st
    end

    # p_φ 和 p_ψ 守恒
    pphidot = 0.0
    ppsidot = 0.0

    return [thetadot, phidot, psidot, pthetadot, pphidot, ppsidot]
end

"""
用 RK4 积分 Lagrange Top 的完整运动方程

参数:
- theta0, phi0, psi0: 初始欧拉角
- theta_dot0: 初始章动角速率
- Lz: 空间 z 轴角动量 (常数)
- L3: 体轴 z 轴角动量 (常数)
- t_end, dt: 积分参数

返回: (times, thetas, phis, psis)
"""
function simulate_lagrange_top(I1::Float64, I3::Float64, M::Float64, g::Float64, l::Float64,
                                theta0::Float64, phi0::Float64, psi0::Float64,
                                theta_dot0::Float64, Lz::Float64, L3::Float64,
                                t_end::Float64, dt::Float64)
    ptheta0 = I1 * theta_dot0
    state = [theta0, phi0, psi0, ptheta0, Lz, L3]

    times = [0.0]
    thetas = [theta0]; phis = [phi0]; psis = [psi0]

    n_steps = Int(ceil(t_end / dt))

    for step in 1:n_steps
        f(s) = lagrange_top_equations(s, I1, I3, M, g, l)

        # RK4
        half_dt = dt / 2.0
        k1 = f(state)
        k2 = f(state + half_dt * k1)
        k3 = f(state + half_dt * k2)
        k4 = f(state + dt * k3)

        state = state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

        push!(times, step * dt)
        push!(thetas, state[1])
        push!(phis, state[2])
        push!(psis, state[3])
    end

    return times, thetas, phis, psis
end

# ============================================================
# 陀螺的章动分析
# ============================================================

"""
章动周期（小振幅近似）

T_nut ≈ 2π √(I1 / (Mgl))   (对于近直立陀螺的小振幅章动)
"""
function nutation_period_small(I1::Float64, M::Float64, g::Float64, l::Float64)
    return 2π * sqrt(I1 / (M * g * l))
end

"""
章动振幅范围

从 V_eff(θ) = E 求解章动角的两个转折点 θ_min 和 θ_max
"""
function nutation_amplitude_range(I1::Float64, I3::Float64, M::Float64, g::Float64,
                                   l::Float64, E::Float64, Lz::Float64, L3::Float64;
                                   n_pts::Int=500)
    thetas = range(1e-4, π - 1e-4, length=n_pts)
    turning_points = Float64[]

    for i in 1:n_pts
        th = thetas[i]
        rate = nutation_rate(I1, I3, M, g, l, E, Lz, L3, th)
        if i > 1
            prev_rate = nutation_rate(I1, I3, M, g, l, E, Lz, L3, thetas[i-1])
            # 零交叉检测
            if prev_rate * rate < 0
                push!(turning_points, th)
            end
        end
        # rate ≈ 0
        if abs(rate) < 1e-10 && length(turning_points) < 2
            if isempty(turning_points) || abs(th - turning_points[end]) > 1e-5
                push!(turning_points, th)
            end
        end
    end

    return unique(turning_points)
end

"""
判断进动类型 (Goldstein Fig.5.7):

- monotonic: 进动单向进行 (φ̇ 不变符号), 出现在 θ 在 Lz>0 区域
- looping: 进动中 θ 周期性地回到同一值 (需要 φ̇ 在某些 θ 改变符号)
- cusped: θ 边界处 φ̇=0 (反转点)
"""
function precession_type(I1::Float64, Lz::Float64, L3::Float64,
                          theta_range::Vector{Float64})
    if length(theta_range) < 2
        return :monotonic
    end

    rates = [precession_rate_top(I1, Lz, L3, th) for th in theta_range]

    # 检查 φ̇ 是否改变符号
    if minimum(rates) < 0 && maximum(rates) > 0
        return :looping
    elseif any(abs.(rates) .< 1e-12)
        return :cusped
    else
        return :monotonic
    end
end

# ============================================================
# 陀螺仪原理
# ============================================================

"""
陀螺仪中的科里奥利力矩

当陀螺仪的支撑被强制旋转时:
  N = ω_forced × L_spin

返回: 陀螺仪感受到的力矩
"""
function gyroscopic_torque(L_spin::AbstractVector, omega_forced::AbstractVector)
    # 陀螺角动量方向试图跟随强制旋转，产生的力矩
    return cross(L_spin, omega_forced)
end

"""
单自由度陀螺仪的进动角速率

在恒定外力矩下的稳态进动:
  Ω_precession = N / L_spin
"""
function gyroscope_precession(N::Float64, L_spin::Float64)
    return abs(L_spin) > 1e-15 ? N / L_spin : 0.0
end

export symmetric_top_effective_potential, nutation_rate
export precession_rate_top, spin_rate_top
export steady_precession_theta, fast_top_precession, slow_precession
export sleeping_top_stability, critical_spin_rate
export lagrange_top_equations, simulate_lagrange_top
export nutation_period_small, nutation_amplitude_range, precession_type
export gyroscopic_torque, gyroscope_precession
