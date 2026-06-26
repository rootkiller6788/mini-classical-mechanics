# types.jl -- 刚体力学核心数据类型
# 参考: Goldstein Ch.4-5, MIT 8.012
# 零外部依赖，手写所有结构体

# ============================================================
# InertiaTensor -- 3x3 惯性张量（对称矩阵）
# ============================================================
"""
惯性张量（3x3 实对称矩阵）

字段:
- Ixx, Iyy, Izz: 对角线分量（转动惯量）
- Ixy, Ixz, Iyz: 非对角线分量（惯性积）

参考: Goldstein (5.3)
"""
struct InertiaTensor
    Ixx::Float64; Iyy::Float64; Izz::Float64
    Ixy::Float64; Ixz::Float64; Iyz::Float64
end

"""
将 InertiaTensor 转为 3x3 Matrix
"""
function inertia_matrix(I::InertiaTensor)
    return [I.Ixx I.Ixy I.Ixz;
            I.Ixy I.Iyy I.Iyz;
            I.Ixz I.Iyz I.Izz]
end

# ============================================================
# EulerAngles -- 欧拉角 (ZXZ convention, 常用于经典力学)
# ============================================================
"""
欧拉角 (ZXZ convention: ϕ=进动, θ=章动, ψ=自转)

参考: Goldstein (4.46)
"""
struct EulerAngles
    phi::Float64      # 进动角 (precession) [0, 2π)
    theta::Float64    # 章动角 (nutation) [0, π]
    psi::Float64      # 自转角 (spin) [0, 2π)
end

EulerAngles() = EulerAngles(0.0, 0.0, 0.0)

# ============================================================
# RigidBodyState -- 刚体运动完整状态
# ============================================================
"""
刚体运动状态（空间坐标系 + 体坐标系）

字段:
- t: 时间
- omega: 角速度向量（体坐标系分量）
- euler: 欧拉角（空间定向）
- L: 角动量（空间坐标系）
- T: 转动动能

参考: Goldstein Ch.5
"""
mutable struct RigidBodyState
    t::Float64
    omega::Vector{Float64}   # 体坐标系角速度 [ωx, ωy, ωz]
    euler::EulerAngles       # 空间定向
    L::Vector{Float64}       # 空间坐标系角动量
    T::Float64               # 转动动能
end

function RigidBodyState(omega::Vector{Float64}, euler::EulerAngles)
    return RigidBodyState(0.0, omega, euler, zeros(3), 0.0)
end

RigidBodyState() = RigidBodyState([0.0,0.0,0.0], EulerAngles())

# ============================================================
# PrincipalAxes -- 主轴系数据
# ============================================================
"""
主轴系:
- moments: 主转动惯量 [I1, I2, I3], 降序排列 I1≥I2≥I3
- axes: 主轴方向 (3x3 正交矩阵，每列为归一化主轴向量)
- transform: 旋转矩阵 R 使 I_principal = R^T * I * R

参考: Goldstein (5.4)-(5.6)
"""
struct PrincipalAxes
    moments::Vector{Float64}   # [I1, I2, I3]
    axes::Matrix{Float64}      # 3x3 正交矩阵
end

# ============================================================
# 几何体描述常量
# ============================================================
"""
标准刚体的惯性张量解析公式（质心系）
"""
struct StandardShape
    name::String
    mass::Float64
    params::Dict{Symbol, Float64}  # 几何参数
end

# ============================================================
# 显示方法
# ============================================================
function Base.show(io::IO, I::InertiaTensor)
    println(io, "InertiaTensor:")
    println(io, "  Ixx=$(I.Ixx), Iyy=$(I.Iyy), Izz=$(I.Izz)")
    print(io,   "  Ixy=$(I.Ixy), Ixz=$(I.Ixz), Iyz=$(I.Iyz)")
end

function Base.show(io::IO, ea::EulerAngles)
    print(io, "EulerAngles(ϕ=$(round(ea.phi,4)), θ=$(round(ea.theta,4)), ψ=$(round(ea.psi,4)))")
end

function Base.show(io::IO, pa::PrincipalAxes)
    println(io, "PrincipalAxes:")
    println(io, "  moments: $(round.(pa.moments, digits=4))")
    print(io,   "  axes: $(round.(pa.axes, digits=4))")
end

export InertiaTensor, EulerAngles, RigidBodyState, PrincipalAxes, StandardShape
export inertia_matrix
