# three_body.jl — 限制性三体问题: Lagrange 点, Jacobi 常数, 零速度面
# 参考: Murray & Dermott Ch.3, Szebehely "Theory of Orbits"

using LinearAlgebra

# ============================================================
# 限制性三体问题 (CRTBP)
# ============================================================

"""
圆型限制性三体问题 (CRTBP) 运动方程

参数:
- state: [x, y, z, vx, vy, vz] (rotating frame)
- mu: 质量参数 m2/(m1+m2)

参考: Szebehely Ch.1
"""
function crtbp_equations(state::Vector{Float64}, mu::Float64)
    x, y, z, vx, vy, vz = state

    r1 = sqrt((x + mu)^2 + y^2 + z^2)      # 到 m1 的距离
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2)  # 到 m2 的距离

    # 有效势 Ω 的导数
    Ux = x - (1 - mu)*(x + mu)/r1^3 - mu*(x - 1 + mu)/r2^3
    Uy = y - (1 - mu)*y/r1^3 - mu*y/r2^3
    Uz = -(1 - mu)*z/r1^3 - mu*z/r2^3

    # 运动方程: ẍ - 2ẏ = ∂Ω/∂x,  ÿ + 2ẋ = ∂Ω/∂y,  z̈ = ∂Ω/∂z
    ax =  2*vy + Ux
    ay = -2*vx + Uy
    az = Uz

    return [vx, vy, vz, ax, ay, az]
end

"""
CRTBP 有效势 — Jacobi 积分的势能部分

Ω(x, y, z) = (x² + y²)/2 + (1-μ)/r1 + μ/r2 + (1-μ)μ/2
"""
function crtbp_effective_potential(x::Float64, y::Float64, z::Float64, mu::Float64)
    r1 = sqrt((x + mu)^2 + y^2 + z^2)
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2)
    return 0.5*(x^2 + y^2) + (1 - mu)/r1 + mu/r2 + 0.5*(1 - mu)*mu
end

"""
Jacobi 常数 (运动常数)

C_J = 2Ω(x,y,z) - v²

在旋转坐标系中守恒。
"""
function jacobi_constant(state::Vector{Float64}, mu::Float64)
    x, y, z, vx, vy, vz = state
    v2 = vx^2 + vy^2 + vz^2
    Omega = crtbp_effective_potential(x, y, z, mu)
    return 2*Omega - v2
end

"""
零速度面: 给定 Jacobi 常数 C_J，解 v² = 2Ω - C_J = 0

零速度面在 (x,y) 平面上的高度计算
"""
function zero_velocity_surface_z(x::Float64, y::Float64, C_J::Float64, mu::Float64; z::Float64=0.0)
    r1 = sqrt((x + mu)^2 + y^2 + z^2)
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2)
    Omega = 0.5*(x^2 + y^2) + (1 - mu)/r1 + mu/r2 + 0.5*(1 - mu)*mu
    arg = 2*Omega - C_J
    return arg >= 0 ? sqrt(arg) : NaN
end

# ============================================================
# Lagrange 点
# ============================================================

"""
所有 5 个 Lagrange 点 (旋转坐标系, xy 平面)

L1-L3: 共线点 (y=0)
L4-L5: 三角点 (等边三角形)

参数:
- mu: m2/(m1+m2)

返回: LagrangePointData
"""
function lagrange_points(mu::Float64)
    # 共线点
    L1 = find_collinear_lagrange(mu, 0.0, 1.0 - mu)   # m1 和 m2 之间
    L2 = find_collinear_lagrange(mu, 1.0 - mu, 2.0)    # m2 之外
    L3 = find_collinear_lagrange(mu, -2.0, -mu)         # m1 之外

    # 三角点 (等边三角形)
    L4_x = 0.5 - mu
    L4_y = sqrt(3.0) / 2.0

    return LagrangePointData(
        [L1, 0.0, 0.0],
        [L2, 0.0, 0.0],
        [L3, 0.0, 0.0],
        [L4_x, L4_y, 0.0],
        [L4_x, -L4_y, 0.0]
    )
end

"二分法求共线 Lagrange 点"
function find_collinear_lagrange(mu::Float64, x0::Float64, x1::Float64; tol::Float64=1e-12)
    # 共线点方程: x - (1-μ)(x+μ)/|x+μ|³ - μ(x-1+μ)/|x-1+μ|³ = 0
    function f(x)
        return x - (1 - mu)*(x + mu)/abs(x + mu)^3 - mu*(x - 1 + mu)/abs(x - 1 + mu)^3
    end

    a, b = x0, x1
    fa, fb = f(a), f(b)

    for _ in 1:60
        mid = (a + b) / 2.0
        fmid = f(mid)
        if abs(fmid) < tol
            return mid
        end
        fa * fmid < 0 ? (b = mid; fb = fmid) : (a = mid; fa = fmid)
    end
    return (a + b) / 2.0
end

"""
Lagrange 点的 Jacobi 常数
"""
function lagrange_jacobi_constants(mu::Float64)
    pts = lagrange_points(mu)
    C_Js = Float64[]
    for pt in [pts.L1, pts.L2, pts.L3, pts.L4, pts.L5]
        Omega = crtbp_effective_potential(pt[1], pt[2], pt[3], mu)
        push!(C_Js, 2*Omega)  # v=0 at L-points
    end
    return C_Js
end

"""
Hill 区域分析: 给定 C_J, 确定哪些区域的运动在零速度面允许范围内

返回: 连通性描述
"""
function hill_region(c_j::Float64, mu::Float64)
    C_l = lagrange_jacobi_constants(mu)

    if c_j > C_l[1]  # > C_L1
        return "m1 or m2 only (disconnected)"
    elseif c_j > C_l[2]  # C_L2 < C_J < C_L1
        return "m1 and m2 connected, outer region closed"
    elseif c_j > C_l[3]  # C_L3 < C_J < C_L2
        return "inner + outer open, m2 exterior closed"
    else  # < C_L3
        return "all regions connected"
    end
end

# ============================================================
# 影响球
# ============================================================

"""
Laplace 影响球半径

r_SOI = a · (m/M)^(2/5)
"""
sphere_of_influence(a::Float64, m::Float64, M::Float64) = a * (m/M)^(2.0/5.0)

"""
Hill 球半径

r_H = a · (m/(3M))^(1/3)
"""
hill_sphere(a::Float64, m::Float64, M::Float64) = a * (m/(3*M))^(1.0/3.0)

"""
Roche 极限 (刚体卫星)

d = 1.26 · R_primary · (ρ_primary/ρ_satellite)^(1/3)
"""
function roche_limit_rigid(R_primary::Float64, rho_primary::Float64, rho_satellite::Float64)
    return 1.26 * R_primary * (rho_primary/rho_satellite)^(1.0/3.0)
end

"Roche 极限 (流体卫星): d = 2.44 · R · (ρ_M/ρ_m)^(1/3)"
function roche_limit_fluid(R_primary::Float64, rho_primary::Float64, rho_satellite::Float64)
    return 2.44 * R_primary * (rho_primary/rho_satellite)^(1.0/3.0)
end

# ============================================================
# CRTBP 数值积分 (简化)
# ============================================================

"""
在 CRTBP 中积分一小段时间

使用 RK4
"""
function integrate_crtbp(state0::Vector{Float64}, mu::Float64, t_end::Float64, dt::Float64)
    state = copy(state0)
    n_steps = Int(ceil(t_end / dt))
    traj = [copy(state)]

    for _ in 1:n_steps
        f(s) = crtbp_equations(s, mu)
        half_dt = dt / 2.0
        k1 = f(state)
        k2 = f(state + half_dt * k1)
        k3 = f(state + half_dt * k2)
        k4 = f(state + dt * k3)
        state = state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        push!(traj, copy(state))
    end

    return traj
end

export crtbp_equations, crtbp_effective_potential, jacobi_constant
export zero_velocity_surface_z
export lagrange_points, lagrange_jacobi_constants, hill_region
export sphere_of_influence, hill_sphere, roche_limit_rigid, roche_limit_fluid
export integrate_crtbp
