# inertia.jl -- 惯性张量计算、主轴变换、平行轴定理
# 参考: Goldstein Ch.5, MIT 8.012 Lecture 24-25

using LinearAlgebra

# ============================================================
# 惯性张量计算
# ============================================================

"""
从粒子集合计算惯性张量（质心系）

参数:
- masses: 各粒子质量
- positions: 各粒子位置（相对于质心）

返回: InertiaTensor

公式 (Goldstein 5.3):
  Ixx = Σ m(y²+z²),  Iyy = Σ m(x²+z²),  Izz = Σ m(x²+y²)
  Ixy = -Σ mxy,      Ixz = -Σ mxz,      Iyz = -Σ myz
"""
function inertia_tensor(masses::Vector{Float64}, positions::Vector{<:AbstractVector})
    Ixx = Iyy = Izz = Ixy = Ixz = Iyz = 0.0
    for (m, r) in zip(masses, positions)
        x, y, z = r[1], r[2], r[3]
        r2 = x*x + y*y + z*z
        Ixx += m * (r2 - x*x)
        Iyy += m * (r2 - y*y)
        Izz += m * (r2 - z*z)
        Ixy -= m * x * y
        Ixz -= m * x * z
        Iyz -= m * y * z
    end
    return InertiaTensor(Ixx, Iyy, Izz, Ixy, Ixz, Iyz)
end

"""
将质心惯性张量平移到另一点（平行轴定理）

参数:
- I_cm: 质心惯性张量
- M: 总质量
- d: 位移向量 (从质心到新参考点)

公式 (Goldstein 5.18):
  I_new = I_cm + M(d²𝟙 - d⊗d)
"""
function parallel_axis(I_cm::InertiaTensor, M::Float64, d::AbstractVector)
    d2 = d[1]*d[1] + d[2]*d[2] + d[3]*d[3]
    return InertiaTensor(
        I_cm.Ixx + M*(d2 - d[1]*d[1]),
        I_cm.Iyy + M*(d2 - d[2]*d[2]),
        I_cm.Izz + M*(d2 - d[3]*d[3]),
        I_cm.Ixy - M*d[1]*d[2],
        I_cm.Ixz - M*d[1]*d[3],
        I_cm.Iyz - M*d[2]*d[3]
    )
end

# ============================================================
# 主轴分解
# ============================================================

"""
计算主轴: 对角化惯性矩阵

返回: PrincipalAxes(moments, axes)
- moments: [I1, I2, I3] 降序排列 I1 ≥ I2 ≥ I3
- axes: 3x3 矩阵，columns 为对应主轴方向

参考: Goldstein (5.4)-(5.6)
"""
function principal_axes(I::InertiaTensor)
    M = inertia_matrix(I)
    evals, evecs = eigen(M)
    # 降序排列
    idx = sortperm(evals, rev=true)
    return PrincipalAxes(evals[idx], evecs[:, idx])
end

"""
主转动惯量（仅返回值，不返回轴）
"""
function principal_moments(I::InertiaTensor)
    return sort(eigvals(inertia_matrix(I)), rev=true)
end

"""
验证惯性张量是否已对角化（是否在主轴系中）
"""
function is_diagonal(I::InertiaTensor; tol::Float64=1e-12)
    return abs(I.Ixy) < tol && abs(I.Ixz) < tol && abs(I.Iyz) < tol
end

# ============================================================
# 标准几何体惯性张量解析公式（质心系）
# ============================================================

"""
均匀球体（半径 R, 质量 M）
Ixx = Iyy = Izz = (2/5) M R^2
"""
function sphere_inertia(M::Float64, R::Float64)
    I = 0.4 * M * R * R  # 2/5
    return InertiaTensor(I, I, I, 0.0, 0.0, 0.0)
end

"""
均匀球壳（半径 R, 质量 M）
Ixx = Iyy = Izz = (2/3) M R^2
"""
function spherical_shell_inertia(M::Float64, R::Float64)
    I = (2/3) * M * R * R
    return InertiaTensor(I, I, I, 0.0, 0.0, 0.0)
end

"""
均匀实心圆柱（半径 R, 高度 H, 质量 M）
对称轴为 z 轴。
Izz = (1/2) M R^2
Ixx = Iyy = (1/12) M (3R^2 + H^2)
"""
function cylinder_inertia(M::Float64, R::Float64, H::Float64)
    Izz = 0.5 * M * R * R
    Iperp = M * (3*R*R + H*H) / 12.0
    return InertiaTensor(Iperp, Iperp, Izz, 0.0, 0.0, 0.0)
end

"""
均匀长方体（边长 a×b×c, 质量 M）
Ixx = (1/12) M (b^2 + c^2)
Iyy = (1/12) M (a^2 + c^2)
Izz = (1/12) M (a^2 + b^2)
"""
function cuboid_inertia(M::Float64, a::Float64, b::Float64, c::Float64)
    return InertiaTensor(
        M*(b*b + c*c)/12.0,
        M*(a*a + c*c)/12.0,
        M*(a*a + b*b)/12.0,
        0.0, 0.0, 0.0
    )
end

"""
均匀细杆（长度 L, 质量 M, 沿 z 轴放置，质心在原点）
Ixx = Iyy = (1/12) M L^2
Izz = 0
"""
function rod_inertia(M::Float64, L::Float64)
    Iperp = M * L * L / 12.0
    return InertiaTensor(Iperp, Iperp, 0.0, 0.0, 0.0, 0.0)
end

"""
均匀薄圆盘（半径 R, 质量 M, 在 xy 平面）
Ixx = Iyy = (1/4) M R^2
Izz = (1/2) M R^2  (对称轴)
"""
function disk_inertia(M::Float64, R::Float64)
    Iperp = 0.25 * M * R * R
    Izz = 0.5 * M * R * R
    return InertiaTensor(Iperp, Iperp, Izz, 0.0, 0.0, 0.0)
end

# ============================================================
# 惯性椭球
# ============================================================

"""
惯性椭球方程系数 (Goldstein 5.11):
ρ^T · I · ρ = Ixx x^2 + Iyy y^2 + Izz z^2 + 2Ixy xy + 2Ixz xz + 2Iyz yz = 1

返回椭球半轴 (a, b, c):
a = 1/√I1, b = 1/√I2, c = 1/√I3
"""
function inertia_ellipsoid(I::InertiaTensor)
    pm = principal_moments(I)
    a = 1.0 / sqrt(pm[1])  # 最短半轴 (最大惯性)
    b = 1.0 / sqrt(pm[2])
    c = 1.0 / sqrt(pm[3])  # 最长半轴 (最小惯性)
    return (a, b, c)
end

"""
生成惯性椭球表面点云（用于可视化）
"""
function inertia_ellipsoid_points(I::InertiaTensor; n_theta::Int=40, n_phi::Int=20)
    pa = principal_axes(I)
    a, b, c = inertia_ellipsoid(I)
    points = Vector{Float64}[]
    for i in 1:n_theta, j in 1:n_phi
        theta = π * i / (n_theta - 1)
        phi = 2π * j / n_phi
        x_local = a * sin(theta) * cos(phi)
        y_local = b * sin(theta) * sin(phi)
        z_local = c * cos(theta)
        global_pt = pa.axes * [x_local, y_local, z_local]
        push!(points, global_pt)
    end
    return points
end

export inertia_tensor, parallel_axis, principal_axes, principal_moments, is_diagonal
export sphere_inertia, spherical_shell_inertia, cylinder_inertia, cuboid_inertia
export rod_inertia, disk_inertia, inertia_ellipsoid, inertia_ellipsoid_points
