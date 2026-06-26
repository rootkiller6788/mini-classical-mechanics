# perturbations.jl — 轨道摄动: J2, 大气阻力, 太阳光压, 第三体
# 参考: Vallado Ch.8-9, Murray & Dermott Ch.10

using LinearAlgebra

# ============================================================
# J2 摄动 (地球扁率)
# ============================================================

"""
J2 摄动加速度 (ECI 坐标系)

参数:
- r: 位置矢量 [km]
- mu: 引力参数
- J2: J2 系数 (地球: 0.00108263)
- R: 中心天体半径 [km]

公式 (Vallado 8.25):
  a_J2 = -3/2 · J2 · μ · R² / r⁵ · [x(1-5z²/r²), y(1-5z²/r²), z(3-5z²/r²)]
"""
function j2_acceleration(r::Vector{Float64}; mu::Float64=398600.5,
                          J2::Float64=0.00108263, R::Float64=6378.1363)
    r_mag = norm(r)
    x, y, z = r
    zr = z / r_mag
    factor = -1.5 * J2 * mu * R^2 / r_mag^5

    return [
        factor * x * (1.0 - 5.0 * zr^2),
        factor * y * (1.0 - 5.0 * zr^2),
        factor * z * (3.0 - 5.0 * zr^2)
    ]
end

"""
J2 引起的长期进动率 (轨道根数变化率)

Ω̇ = -3/2 · n · J2 · (R/p)² · cos i
ω̇ = 3/4 · n · J2 · (R/p)² · (5 cos² i - 1)
Ṁ = n + 3/4 · n · J2 · (R/p)² · √(1-e²) · (3 cos² i - 1)

返回: (Omega_dot, omega_dot, M_dot) [rad/s]

参考: Vallado 8.29
"""
function j2_secular_rates(a::Float64, e::Float64, i::Float64;
                           mu::Float64=398600.5, J2::Float64=0.00108263,
                           R::Float64=6378.1363)
    n = sqrt(mu / a^3)
    p = a * (1.0 - e^2)
    Rp2 = (R / p)^2
    factor = 1.5 * n * J2 * Rp2

    Omega_dot = -factor * cos(i)
    omega_dot = 0.5 * factor * (5.0 * cos(i)^2 - 1.0)
    M_dot = n + 0.5 * factor * sqrt(1.0 - e^2) * (3.0 * cos(i)^2 - 1.0)

    return (Omega_dot=Omega_dot, omega_dot=omega_dot, M_dot=M_dot)
end

"""
太阳同步轨道条件: Ω̇ = 2π/年 ≈ 360°/365.25天

对于给定 a, e 求解 i 使得 Ω̇ = Ω̇_target
"""
function sun_sync_inclination(a::Float64, e::Float64;
                               mu::Float64=398600.5, J2::Float64=0.00108263,
                               R::Float64=6378.1363)
    n = sqrt(mu / a^3)
    p = a * (1.0 - e^2)
    omega_sun = 2π / 365.25 / 86400  # [rad/s]

    # omega_sun = -1.5 * n * J2 * (R/p)^2 * cos(i)
    cos_i = -omega_sun / (1.5 * n * J2 * (R/p)^2)
    cos_i = clamp(cos_i, -1.0, 1.0)
    return acos(cos_i)
end

# ============================================================
# 大气阻力
# ============================================================

"""
大气阻力加速度 (简化指数大气模型)

参数:
- r: 位置矢量 [km]
- v: 速度矢量 [km/s]
- C_D: 阻力系数 (~2.2)
- A_over_m: 面积质量比 [km²/kg]
- rho0: 参考密度 [kg/km³]
- r0: 参考半径 [km]
- H: 大气标高 [km]

公式: a_drag = -½ C_D (A/m) ρ v v_rel
      ρ = ρ0 exp(-(r - r0)/H)
"""
function atmospheric_drag_acceleration(r::Vector{Float64}, v::Vector{Float64};
                                        C_D::Float64=2.2, A_over_m::Float64=0.01,
                                        rho0::Float64=1.225e12, r0::Float64=6378.0,
                                        H::Float64=8.5)
    r_mag = norm(r)
    alt = r_mag - r0
    if alt < 0
        return zeros(3)  # 已再入
    end
    rho = rho0 * exp(-alt / H)
    v_mag = norm(v)
    v_hat = v / v_mag
    return -0.5 * C_D * A_over_m * rho * v_mag^2 * v_hat
end

"""
大气阻力导致的轨道衰降率 (圆轨道近似)

da/dt ≈ -C_D (A/m) ρ n a²

返回: da/dt [km/s]
"""
function drag_decay_rate(a::Float64, C_D::Float64, A_over_m::Float64, rho::Float64;
                          mu::Float64=398600.5)
    n = sqrt(mu / a^3)
    return -C_D * A_over_m * rho * n * a^2
end

# ============================================================
# 太阳光压 (Solar Radiation Pressure)
# ============================================================

"""
太阳光压加速度

参数:
- r_sat: 卫星位置 [km]
- r_sun: 太阳位置 [km]
- C_R: 反射系数 (1.0 完全吸收, 2.0 完全反射)
- A_over_m: 面积质量比 [km²/kg]
- P_sun: 太阳光压 (4.56e-6 N/m² = 4.56e-15 N/km²)

公式: a_SRP = -C_R (A/m) P_sun AU² · (r_sat - r_sun) / |r_sat - r_sun|³
"""
function srp_acceleration(r_sat::Vector{Float64}, r_sun::Vector{Float64};
                           C_R::Float64=1.5, A_over_m::Float64=0.01,
                           P_sun::Float64=4.56e-15)
    delta = r_sat - r_sun
    d_mag = norm(delta)
    if d_mag < 1e-10
        return zeros(3)
    end
    factor = -C_R * A_over_m * P_sun * (149597870.7^2) / d_mag^3
    return factor * delta
end

# ============================================================
# 第三体摄动
# ============================================================

"""
第三体引力摄动加速度

参数:
- r_sat: 卫星位置
- r_body: 摄动体位置
- mu_body: 摄动体 GM

公式: a_3b = -μ{ (r_sat - r_body)/|r_sat - r_body|³ + r_body/|r_body|³ }
"""
function third_body_acceleration(r_sat::Vector{Float64}, r_body::Vector{Float64},
                                  mu_body::Float64)
    delta = r_sat - r_body
    d_mag = norm(delta)
    r_body_mag = norm(r_body)

    if d_mag < 1e-10
        return zeros(3)
    end

    direct = delta / d_mag^3
    indirect = r_body / r_body_mag^3
    return -mu_body * (direct - indirect)
end

# ============================================================
# 长期进动 (合并所有摄动)
# ============================================================

"""
GR (广义相对论) 近心点进动率

ω̇_GR = 3 n · μ / (a · c² · (1-e²))

用于水星进动的经典验证。

参数:
- a: 半长径 [AU]
- e: 离心率
- mu: 引力参数 (太阳单位)
- c: 光速 [AU/yr] (≈ 63239.7 AU/yr)

返回: 进动率 [rad/yr]

参考: Goldstein Ch.7, Weinberg
"""
function gr_precession_rate(a::Float64, e::Float64; mu::Float64=G_SUN,
                             c::Float64=63239.7)
    n = sqrt(mu / a^3)
    return 3.0 * n * mu / (a * c^2 * (1.0 - e^2))
end

"""
经典进动率 (忽略 GR, 仅用于参考)
"""
function precession_rate_classical(a::Float64, m_perturber::Float64, M_central::Float64;
                                    mu::Float64=G_SUN)
    n = sqrt(mu / a^3)
    return 3.0 * n * m_perturber / M_central
end

# ============================================================
# 摄动轨道传播 (简化)
# ============================================================

"""
含 J2 摄动的简化轨道传播 (仅长期项)

使用平均轨道根数的长期变化率更新根数
"""
function propagate_with_j2(el0::OrbitalElements, delta_t::Float64;
                            mu::Float64=398600.5, J2::Float64=0.00108263,
                            R::Float64=6378.1363)
    rates = j2_secular_rates(el0.a, el0.e, el0.i; mu=mu, J2=J2, R=R)

    return OrbitalElements(
        el0.a,                            # a 无长期变化
        el0.e,                            # e 无长期变化
        el0.i,                            # i 无长期变化
        el0.Omega + rates.Omega_dot * delta_t,
        el0.omega + rates.omega_dot * delta_t,
        el0.nu + rates.M_dot * delta_t    # 近似: M_dot ≈ ν̇
    )
end

export j2_acceleration, j2_secular_rates, sun_sync_inclination
export atmospheric_drag_acceleration, drag_decay_rate
export srp_acceleration, third_body_acceleration
export gr_precession_rate, precession_rate_classical
export propagate_with_j2
