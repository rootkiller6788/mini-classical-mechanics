# central_force.jl — 中心力问题与 Kepler 轨道
# 参考：Goldstein Ch.3, Landau Vol.1 Ch.3-4, MIT 8.012 Lecture 25
# L4 Fundamental Laws: Kepler's laws, conservation of angular momentum
# L6 Canonical Systems: Kepler orbits, scattering

"""
中心力问题: 质点受径向力 F = f(r) r̂, 势能 U = U(r).

守恒量:
  角动量 ℓ = m r² θ̇ (常数)
  能量   E = ½m(ṙ² + r²θ̇²) + U(r) = ½mṙ² + U_eff(r)
  有效势 U_eff(r) = U(r) + ℓ²/(2mr²)
"""
struct CentralForceSystem
    m::Float64               # 质量
    U::Function              # U(r) 势能函数
    dU_dr::Function          # dU/dr
    name::String             # 系统名称
end

"""有效势: U_eff(r) = U(r) + ℓ²/(2mr²)"""
function effective_potential(sys::CentralForceSystem, r::Float64, L::Float64)
    return sys.U(r) + L^2 / (2*sys.m*r^2)
end

"""有效势的导数: dU_eff/dr = dU/dr - ℓ²/(mr³)"""
function effective_potential_derivative(sys::CentralForceSystem, r::Float64, L::Float64)
    return sys.dU_dr(r) - L^2 / (sys.m * r^3)
end

"""
径向运动方程:
  m r̈ = ℓ²/(mr³) - dU/dr
转换为: r̈ = ℓ²/(m²r³) - (1/m) dU/dr
"""
function radial_acceleration(sys::CentralForceSystem, r::Float64, L::Float64)
    return L^2 / (sys.m^2 * r^3) - sys.dU_dr(r) / sys.m
end

"""
Binet 公式: 将轨道方程从 (r, θ) 变换到 (u=1/r, θ).
d²u/dθ² + u = -F(1/u) / (m ℓ² u²)
其中 F = -dU/dr 是力的大小.

对于幂律力 F = -k rⁿ:
  d²u/dθ² + u = k u^{-n-2} / (m ℓ²)
"""
function binet_formula(sys::CentralForceSystem, u::Float64, L::Float64)
    n = 2  # 引力/库仑力情况 k/r² → dU/dr = k/r², F = -k/r²
    r = 1/u
    F_r = -sys.dU_dr(r)  # 力的大小
    return -F_r / (sys.m * L^2 * u^2) - u
end

"""
Kepler 问题: U = -k/r (引力, k > 0).

解析轨道 (圆锥曲线):
  r(θ) = p / (1 + e cos(θ - θ₀))
  半通径 p = ℓ²/(mk)
  偏心率 e = √(1 + 2Eℓ²/(mk²))
  半长轴 a = p/(1-e²) = k/(2|E|)

e=0: 圆, 0<e<1: 椭圆, e=1: 抛物线, e>1: 双曲线
"""
struct KeplerOrbit
    k::Float64     # 力常数 GMm (或 k*q1*q2)
    p::Float64     # 半通径
    e::Float64     # 偏心率
    a::Float64     # 半长轴 (椭圆/双曲线)
    theta0::Float64 # 近心点方向
    T::Float64     # 周期 (for e<1, NaN otherwise)
end

"""
从能量和角动量构造 Kepler 轨道.
"""
function kepler_orbit(m::Float64, k::Float64, E::Float64, L::Float64, theta0::Float64=0.0)
    p = L^2 / (m * k)
    discriminant = 1 + 2*E*L^2/(m*k^2)
    if discriminant < 0
        error("Negative eccentricity squared: invalid (E,L) for Kepler")
    end
    e = sqrt(discriminant)
    if e < 1.0
        a = p / (1 - e^2)
        T = 2π * sqrt(m * a^3 / k)  # Kepler 第三定律
    elseif e == 1.0
        a = Inf; T = Inf
    else
        a = p / (e^2 - 1)
        T = NaN
    end
    return KeplerOrbit(k, p, e, a, theta0, T)
end

"""
计算给定 θ 处的轨道半径 r(θ) = p / (1 + e cos(θ-θ₀)).
"""
function orbit_radius(orb::KeplerOrbit, theta::Float64)
    return orb.p / (1 + orb.e * cos(theta - orb.theta0))
end

"""
计算轨道的笛卡尔坐标.
"""
function orbit_cartesian(orb::KeplerOrbit, thetas::AbstractVector)
    n = length(thetas); xs = zeros(n); ys = zeros(n)
    for i in 1:n
        r = orbit_radius(orb, thetas[i])
        xs[i] = r * cos(thetas[i]); ys[i] = r * sin(thetas[i])
    end
    return xs, ys
end

"""
Kepler 方程: 对于椭圆轨道, M = E_ecc - e sin(E_ecc)
M = 平均近点角 = 2πt/T
E_ecc = 偏心近点角
用 Newton-Raphson 求解.
"""
function solve_kepler_equation(e::Float64, M::Float64; tol=1e-12, max_iter=50)
    E = M  # 初始猜测
    for _ in 1:max_iter
        dE = (E - e*sin(E) - M) / (1 - e*cos(E))
        E -= dE
        if abs(dE) < tol; return E; end
    end
    return E
end

"""
从偏心近点角计算位置: r = a(1 - e cos E); true anomaly θ.
"""
function eccentric_to_true_anomaly(e::Float64, E::Float64)
    x = cos(E) - e; y = sqrt(1 - e^2) * sin(E)
    return atan(y, x)
end

"""
求解中心力运动方程 (数值).
状态: [r, θ, ṙ, θ̇]
返回 (times, states) 用于后续分析.
"""
function solve_central_force(sys::CentralForceSystem, r0::Float64, dr0::Float64, theta0::Float64, dtheta0::Float64, t_end::Float64, h::Float64)
    L0 = sys.m * r0^2 * dtheta0  # 角动量
    function f(t, y)
        r, theta, dr, dtheta = y[1], y[2], y[3], y[4]
        L = sys.m * r^2 * dtheta
        ddr = radial_acceleration(sys, r, L)
        # dtheta = L/(m*r^2), but we use the state variable directly
        return [dr, dtheta, ddr, -2*dr*dtheta/r]  # d²θ/dt² from angular mom conservation
    end
    y0 = [r0, theta0, dr0, dtheta0]
    times, states = rk4_integrate(f, y0, (0.0, t_end), h)
    return times, states
end

"""
验证 Kepler 三定律:
1. 轨道是椭圆 (太阳在一个焦点)
2. 面积速度恒定 (角动量守恒)
3. T² ∝ a³
"""
function verify_kepler_laws(times::Vector{Float64}, states::Vector{Vector{Float64}}, m::Float64, k::Float64)
    n = length(times)
    # 定律 2: 面积速度 dA/dt = ½r²θ̇ = L/(2m) = const
    L_vals = [m * states[i][1]^2 * states[i][4] for i in 1:n]
    L_mean = mean(L_vals); L_std = std(L_vals)
    L_conserved = L_std / abs(L_mean) < 1e-3

    # 定律 1: 拟合椭圆
    xs = [states[i][1]*cos(states[i][2]) for i in 1:n]
    ys = [states[i][1]*sin(states[i][2]) for i in 1:n]
    r_vals = [states[i][1] for i in 1:n]; theta_vals = [states[i][2] for i in 1:n]
    r_min, r_max = minimum(r_vals), maximum(r_vals)
    e_fit = (r_max - r_min) / (r_max + r_min)

    # 定律 3: 从轨道估计周期 T ≈ t_end * n_orbits
    # 检测 r 的周期性
    crosses = 0
    for i in 2:n
        if states[i-1][1] > (r_min+r_max)/2 && states[i][1] <= (r_min+r_max)/2
            crosses += 1
        end
    end
    T_est = crosses > 1 ? times[end] / ((crosses-1)/2) : NaN
    a_est = (r_min + r_max) / 2
    T2_over_a3 = T_est^2 / a_est^3
    kepler3_theory = 4π^2 * m / k

    return (L_conserved=L_conserved, eccentricity=e_fit,
            semi_major_axis=a_est, period_est=T_est,
            T2_over_a3=T2_over_a3, theory=kepler3_theory)
end

"""
Laplace-Runge-Lenz 矢量: A = p × L - m k r̂
对于 Kepler 问题, A 是守恒量, 指向近心点.
|A| = m k e
"""
function laplace_runge_lenz(m::Float64, k::Float64, r::Vector{Float64}, v::Vector{Float64})
    # 2D/3D: r = [x, y] 或 [x, y, z]
    dim = length(r)
    Lz = m * (r[1]*v[2] - r[2]*v[1])  # 角动量 z 分量
    r_norm = norm(r)
    A_x = v[2]*Lz - k*r[1]/r_norm   # p_y*L_z - m*k*x/r
    A_y = -v[1]*Lz - k*r[2]/r_norm   # -p_x*L_z - m*k*y/r
    if dim == 3
        # 全 3D LRL 矢量
        L_vec = m * cross(r, v)
        A = cross(v, L_vec) .- (k/r_norm) .* r
        return A
    else
        return [A_x, A_y]
    end
end

"""
散射角 (Rutherford 散射).

对于排斥库仑势 U = k/r (k > 0):
  b = 瞄准距离, E = 入射能量
  散射角: Θ = π - 2∫_{r_min}^∞ (b/r²) / √(1 - b²/r² - U/E) dr

解析结果 (k>0 排斥, k<0 吸引):
  cot(Θ/2) = 2Eb/k
"""
function scattering_angle(b::Float64, E::Float64, k::Float64, m::Float64)
    v0 = sqrt(2*E/m)
    L = m * v0 * b
    # 最近距离 r_min: E = L²/(2mr²) + k/r
    if k > 0  # 排斥
        r_min = (k + sqrt(k^2 + 2*E*L^2/m)) / (2*E)
    else  # 吸引
        r_min = (abs(k) + sqrt(k^2 + 2*E*L^2/m)) / (2*E)
    end
    # cot(θ/2) = 2Eb/k
    theta = 2 * acot(2*E*b/k)
    return theta, r_min
end

"""
微分散射截面: dσ/dΩ = (b/sinθ) |db/dθ|
对于 Rutherford: dσ/dΩ = (k/(4E))² / sin⁴(θ/2)
"""
function differential_cross_section(theta::Float64, E::Float64, k::Float64)
    return (k/(4*E))^2 / sin(theta/2)^4
end

"""
验证 Bertrand 定理: 只有 k r² (谐振子) 和 k/r (Kepler) 产生闭合轨道.

对一般幂律势 U ∝ r^{n+1}, 检查近圆轨道的径向振荡频率 ω_r.
ω_r² = (3/n) * dU/dr / (mr) + d²U/dr²/m
当 ω_r/ω_θ 为有理数时轨道闭合.
对于 SHO (n=1): ω_r/ω_θ = 2 (闭合)
对于 Kepler (n=-2): ω_r/ω_θ = 1 (闭合)
其他幂律不闭合.
"""
function bertrand_check(U_func::Function, dU_func::Function, d2U_func::Function, r0::Float64, m::Float64)
    L_sq = m * r0^3 * dU_func(r0)  # 圆轨道条件: dU_eff/dr = 0
    omega_theta_sq = L_sq / (m^2 * r0^4)
    omega_r_sq = d2U_func(r0)/m + 3*dU_func(r0)/(m*r0)
    ratio = sqrt(omega_r_sq / omega_theta_sq)
    return omega_r_sq, omega_theta_sq, ratio
end
