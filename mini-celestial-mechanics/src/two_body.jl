# two_body.jl — 二体问题: vis-viva, 轨道传播, 飞行时间
# 参考: Goldstein Ch.3, Murray & Dermott Ch.2

# ============================================================
# Vis-Viva 方程与能量
# ============================================================

"""
Vis-Viva 方程: v² = μ(2/r - 1/a)

参数:
- r: 当前位置半径
- a: 半长径 (双曲线为负)
- mu: 引力参数
"""
vis_viva(r::Float64, a::Float64; mu::Float64=G_SUN) = sqrt(mu * (2.0/r - 1.0/a))

"""
轨道能量: ε = -μ/(2a) (椭圆<0, 抛物线=0, 双曲线>0)
"""
orbital_energy(a::Float64; mu::Float64=G_SUN) = -mu / (2.0 * a)

"""
轨道周期: T = 2π√(a³/μ)
"""
orbital_period(a::Float64; mu::Float64=G_SUN) = 2π * sqrt(a^3 / mu)

"""
平均角速度 (mean motion): n = √(μ/a³) = 2π/T
"""
mean_motion(a::Float64; mu::Float64=G_SUN) = sqrt(mu / a^3)

# ============================================================
# 轨道中的位置/速度函数
# ============================================================

"""
椭圆轨道的径向距离: r(ν) = a(1-e²)/(1+e cos ν)
"""
radial_distance_elliptic(a::Float64, e::Float64, nu::Float64) = a * (1.0 - e^2) / (1.0 + e * cos(nu))

"""
椭圆轨道的位置 (轨道平面内, x 指向近心点)
"""
function position_in_orbital_plane(a::Float64, e::Float64, nu::Float64)
    r_mag = radial_distance_elliptic(a, e, nu)
    return [r_mag * cos(nu), r_mag * sin(nu), 0.0]
end

"""
椭圆轨道的速度 (轨道平面内)
"""
function velocity_in_orbital_plane(a::Float64, e::Float64, nu::Float64; mu::Float64=G_SUN)
    r_mag = radial_distance_elliptic(a, e, nu)
    p = a * (1.0 - e^2)
    h_mag = sqrt(mu * p)
    return [
        -mu / h_mag * sin(nu),
        mu / h_mag * (e + cos(nu)),
        0.0
    ]
end

"""
半通径: p = a(1-e²)，对于抛物线: p = h²/μ
"""
semi_latus_rectum(a::Float64, e::Float64) = a * (1.0 - e^2)

# ============================================================
# 飞行时间
# ============================================================

"""
椭圆弧飞行时间 (Kepler 方程)

Δt = √(a³/μ) · [ (E₂ - e sin E₂) - (E₁ - e sin E₁) ]

参数:
- a, e: 轨道根数
- nu1, nu2: 起点和终点的真近点角
- mu: 引力参数

返回: Δt
"""
function time_of_flight_elliptic(a::Float64, e::Float64, nu1::Float64, nu2::Float64;
                                  mu::Float64=G_SUN)
    E1 = true_to_eccentric_anomaly(nu1, e)
    E2 = true_to_eccentric_anomaly(nu2, e)
    M1 = eccentric_to_mean_anomaly(E1, e)
    M2 = eccentric_to_mean_anomaly(E2, e)
    n = mean_motion(a; mu=mu)
    delta_M = mod(M2 - M1 + 2π, 2π)
    return delta_M / n
end

"""
从平近点角计算经过的时间

t = M / n (从近心点开始)
"""
function time_since_periapsis(a::Float64, e::Float64, nu::Float64; mu::Float64=G_SUN)
    E = true_to_eccentric_anomaly(nu, e)
    M = eccentric_to_mean_anomaly(E, e)
    n = mean_motion(a; mu=mu)
    return M / n
end

"""
给定经过的时间，计算真近点角 (Kepler 传播)

传播初始真近点角 ν₀ 到 Δt 之后
"""
function propagate_anomaly(a::Float64, e::Float64, nu0::Float64, delta_t::Float64;
                            mu::Float64=G_SUN)
    E0 = true_to_eccentric_anomaly(nu0, e)
    M0 = eccentric_to_mean_anomaly(E0, e)
    n = mean_motion(a; mu=mu)
    M = M0 + n * delta_t
    M = mod(M, 2π)
    E = solve_kepler(M, e)
    return eccentric_to_true_anomaly(E, e)
end

# ============================================================
# 近心点/远心点
# ============================================================

"近心距: r_p = a(1-e)"
periapsis_distance(a::Float64, e::Float64) = a * (1.0 - e)

"远心距: r_a = a(1+e) (仅 e<1)"
apoapsis_distance(a::Float64, e::Float64) = a * (1.0 + e)

"近心点速度"
periapsis_velocity(a::Float64, e::Float64; mu::Float64=G_SUN) = sqrt(mu / a * (1.0 + e) / (1.0 - e))

"远心点速度"
apoapsis_velocity(a::Float64, e::Float64; mu::Float64=G_SUN) = sqrt(mu / a * (1.0 - e) / (1.0 + e))

"逃逸速度: v_esc = √(2μ/r)"
escape_velocity(r::Float64; mu::Float64=G_SUN) = sqrt(2.0 * mu / r)

"圆轨道速度: v_circ = √(μ/r)"
circular_velocity(r::Float64; mu::Float64=G_SUN) = sqrt(mu / r)

"""
轨道速度分解 (径向、横向)

v_r = μ/(h) e sin ν    (径向)
v_θ = μ/(h) (1 + e cos ν) = h/r  (横向)
"""
function velocity_components(a::Float64, e::Float64, nu::Float64; mu::Float64=G_SUN)
    p = semi_latus_rectum(a, e)
    h_mag = sqrt(mu * p)
    r_mag = radial_distance_elliptic(a, e, nu)
    v_r = mu / h_mag * e * sin(nu)
    v_theta = h_mag / r_mag
    return (v_r=v_r, v_theta=v_theta)
end

# ============================================================
# 轨道特性
# ============================================================

"""
轨道离心率矢量 Laplace-Runge-Lenz: e = (v×h)/μ - r̂
"""
function eccentricity_vector(r::Vector{Float64}, v::Vector{Float64}; mu::Float64=G_SUN)
    h = cross(r, v)
    return cross(v, h) / mu - r / norm(r)
end

"""
轨道半长径 (从状态向量): a = 1/(2/r - v²/μ)
"""
function semi_major_axis_from_state(r::Vector{Float64}, v::Vector{Float64}; mu::Float64=G_SUN)
    r_mag = norm(r)
    v_mag = norm(v)
    return 1.0 / (2.0/r_mag - v_mag^2/mu)
end

"""
轨道离心率 (从状态向量)
"""
function eccentricity_from_state(r::Vector{Float64}, v::Vector{Float64}; mu::Float64=G_SUN)
    return norm(eccentricity_vector(r, v; mu=mu))
end

export vis_viva, orbital_energy, orbital_period, mean_motion
export radial_distance_elliptic, position_in_orbital_plane, velocity_in_orbital_plane
export semi_latus_rectum
export time_of_flight_elliptic, time_since_periapsis, propagate_anomaly
export periapsis_distance, apoapsis_distance, periapsis_velocity, apoapsis_velocity
export escape_velocity, circular_velocity, velocity_components
export eccentricity_vector, semi_major_axis_from_state, eccentricity_from_state
