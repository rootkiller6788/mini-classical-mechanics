# kepler.jl — Kepler 方程求解与轨道根数转换
# 参考: Goldstein Ch.3, Murray & Dermott Ch.2

using LinearAlgebra

# ============================================================
# Kepler 方程求解
# ============================================================

"""
解 Kepler 方程: M = E - e sin(E) (椭圆)

使用 Newton-Raphson 迭代法。

参数:
- M: 平近点角 [rad]
- e: 离心率
- tol: 收敛精度

返回: E 偏近点角 [rad]

参考: Murray & Dermott (2.68)
"""
function solve_kepler(M::Float64, e::Float64; tol::Float64=1e-12)
    # 好的初始猜测 (Danby 1988)
    E = M + e * sin(M) / (1.0 - sin(M + e) + sin(M))
    E = clamp(E, M - e, M + e)
    for _ in 1:30
        dE = (E - e*sin(E) - M) / (1.0 - e*cos(E))
        E -= dE
        abs(dE) < tol && break
    end
    return E
end

"""
解 Kepler 方程 (双曲线): M_h = e sinh(H) - H

参数:
- M_h: 双曲线平近点角
- e: 离心率 (>1)
- tol: 收敛精度

返回: H
"""
function solve_kepler_hyperbolic(M_h::Float64, e::Float64; tol::Float64=1e-12)
    # 初始猜测
    H = e > 1.5 ? log(2 * M_h / e + 1.8) : asinh(M_h / e)
    for _ in 1:30
        dH = (e * sinh(H) - H - M_h) / (e * cosh(H) - 1.0)
        H -= dH
        abs(dH) < tol && break
    end
    return H
end

"""
解 Kepler 方程 (抛物线): Barker 方程

对于 e=1:
  tan(ν/2) + (1/3)tan^3(ν/2) = √(μ/p³) · t = M_p

使用解析公式直接求解。

参数:
- M_p: 抛物线平近点角
- tol: 收敛精度

返回: 真近点角 ν

参考: Barker (1757)
"""
function solve_kepler_parabolic(M_p::Float64; tol::Float64=1e-12)
    # Barker 方程: 2p tan(ν/2) + (2/3)p tan³(ν/2) = 3t √(μ/p)
    # 简化: A = (3/2)M_p
    A = 1.5 * M_p
    B = (A + sqrt(A^2 + 1.0))^(1/3)
    tan_half_nu = B - 1/B
    return 2.0 * atan(tan_half_nu)
end

"""
通用 Kepler 方程求解 (自动选择椭圆/双曲线/抛物线)
"""
function solve_kepler_universal(x::Float64, e::Float64, t::Float64, mu::Float64=G_SUN)
    if e < 1.0 - 1e-12
        return solve_kepler(x, e)  # E
    elseif e > 1.0 + 1e-12
        return solve_kepler_hyperbolic(x, e)  # H
    else
        return solve_kepler_parabolic(x)  # ν
    end
end

# ============================================================
# 真近点角 ↔ 偏近点角 ↔ 平近点角
# ============================================================

"""
偏近点角 → 真近点角 (椭圆)

tan(ν/2) = √((1+e)/(1-e)) · tan(E/2)
"""
function eccentric_to_true_anomaly(E::Float64, e::Float64)
    return 2.0 * atan(sqrt((1.0 + e) / (1.0 - e)) * tan(E / 2.0))
end

"""
真近点角 → 偏近点角 (椭圆)

tan(E/2) = √((1-e)/(1+e)) · tan(ν/2)
"""
function true_to_eccentric_anomaly(nu::Float64, e::Float64)
    return 2.0 * atan(sqrt((1.0 - e) / (1.0 + e)) * tan(nu / 2.0))
end

"""
偏近点角 → 平近点角: M = E - e sin(E)
"""
eccentric_to_mean_anomaly(E::Float64, e::Float64) = E - e * sin(E)

"""
真近点角 → 平近点角
"""
function true_to_mean_anomaly(nu::Float64, e::Float64)
    E = true_to_eccentric_anomaly(nu, e)
    return eccentric_to_mean_anomaly(E, e)
end

"""
平近点角 → 真近点角
"""
function mean_to_true_anomaly(M::Float64, e::Float64)
    E = solve_kepler(M, e)
    return eccentric_to_true_anomaly(E, e)
end

# ============================================================
# 轨道根数 ↔ 笛卡尔状态向量
# ============================================================

"""
轨道根数 → 笛卡尔状态向量 (r, v)

参数:
- a, e, i, Omega, omega, nu: 经典轨道根数
- mu: 引力参数

算法 (Murray & Dermott 2.122):
  1. 计算轨道平面中的位置/速度
  2. 通过 3-1-3 旋转转到空间坐标系

参考: Murray & Dermott Ch.2.8
"""
function orbital_elements_to_state(a::Float64, e::Float64, i::Float64,
                                    Omega::Float64, omega::Float64,
                                    M_or_nu::Float64;
                                    mu::Float64=G_SUN, use_M::Bool=false)
    # 确定真近点角
    if use_M
        nu = mean_to_true_anomaly(M_or_nu, e)
    else
        nu = M_or_nu
    end

    # 轨道平面中的位置 (x 指向近心点方向)
    r_mag = a * (1.0 - e^2) / (1.0 + e * cos(nu))
    r_orb = [r_mag * cos(nu), r_mag * sin(nu), 0.0]

    # 轨道平面中的速度
    p = a * (1.0 - e^2)
    h_mag = sqrt(mu * p)
    v_orb = [
        -mu / h_mag * sin(nu),
        mu / h_mag * (e + cos(nu)),
        0.0
    ]

    # 旋转矩阵: R_z(Omega) · R_x(i) · R_z(omega)
    Rz_Omega = [cos(Omega) -sin(Omega) 0; sin(Omega) cos(Omega) 0; 0 0 1]
    Rx_i     = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)]
    Rz_omega = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1]

    R = Rz_Omega * Rx_i * Rz_omega

    return R * r_orb, R * v_orb
end

"""
笛卡尔状态向量 → 轨道根数

参数:
- r, v: 位置/速度矢量
- mu: 引力参数

算法:
  1. 计算角动量矢量 h = r × v
  2. 从 h 确定倾角和升交点
  3. 从能量确定半长径
  4. 从 Laplace-Runge-Lenz 矢量确定离心率矢量

参考: Murray & Dermott (2.122逆向)
"""
function state_to_orbital_elements(r::Vector{Float64}, v::Vector{Float64};
                                    mu::Float64=G_SUN)
    # 角动量
    h_vec = cross(r, v)
    h_mag = norm(h_vec)

    # 倾角
    i_val = h_mag > 0 ? acos(h_vec[3] / h_mag) : 0.0

    # 升交点矢量 n = k × h
    n_vec = cross([0.0, 0.0, 1.0], h_vec)
    n_mag = norm(n_vec)

    # 升交点赤经
    if n_mag > 1e-15
        Omega = acos(clamp(n_vec[1] / n_mag, -1.0, 1.0))
        if n_vec[2] < 0; Omega = 2π - Omega; end
    else
        Omega = 0.0
    end

    # 半长径 (从 vis-viva: v² = μ(2/r - 1/a))
    r_mag = norm(r)
    v_mag = norm(v)
    a = 1.0 / (2.0/r_mag - v_mag^2/mu)

    # 离心率矢量 e_vec = (v×h)/μ - r̂
    e_vec = cross(v, h_vec) / mu - r / r_mag
    e = norm(e_vec)

    # 近心点幅角
    if e > 1e-12 && n_mag > 1e-15
        omega = acos(clamp(dot(n_vec, e_vec) / (n_mag * e), -1.0, 1.0))
        if e_vec[3] < 0; omega = 2π - omega; end
    else
        omega = 0.0
    end

    # 真近点角
    if e > 1e-12
        nu = acos(clamp(dot(e_vec, r) / (e * r_mag), -1.0, 1.0))
        if dot(r, v) < 0; nu = 2π - nu; end
    else
        # 圆轨道: 用纬度幅角
        nu = n_mag > 1e-15 ? acos(clamp(dot(n_vec, r) / (n_mag * r_mag), -1.0, 1.0)) : 0.0
        if r[3] < 0; nu = 2π - nu; end
    end

    return OrbitalElements(a, e, i_val, Omega, omega, nu)
end

"""
根数→状态→根数 往返验证

返回: (max_error_in_a, max_error_in_e, ...)
"""
function kepler_roundtrip_error(el::OrbitalElements; mu::Float64=G_SUN)
    r, v = orbital_elements_to_state(el.a, el.e, el.i, el.Omega, el.omega, el.nu; mu=mu)
    el2 = state_to_orbital_elements(r, v; mu=mu)
    return (a_err=abs(el.a - el2.a), e_err=abs(el.e - el2.e),
            i_err=abs(el.i - el2.i), Omega_err=abs(el.Omega - el2.Omega),
            omega_err=abs(el.omega - el2.omega))
end

export solve_kepler, solve_kepler_hyperbolic, solve_kepler_parabolic
export solve_kepler_universal
export eccentric_to_true_anomaly, true_to_eccentric_anomaly
export eccentric_to_mean_anomaly, true_to_mean_anomaly, mean_to_true_anomaly
export orbital_elements_to_state, state_to_orbital_elements, kepler_roundtrip_error
