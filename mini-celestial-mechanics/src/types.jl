# types.jl — 天体力学核心数据类型
# 参考: Goldstein Ch.3, Murray & Dermott Ch.2

# ============================================================
# OrbitalElements — 经典轨道根数 (Keplerian elements)
# ============================================================
"""
经典轨道六根数

参数:
- a: 半长径 (>0 椭圆, <0 双曲线, =∞抛物线用 p)
- e: 离心率 (0≤e<1 椭圆, e=1 抛物线, e>1 双曲线)
- i: 倾角 [0, π]
- Omega: 升交点赤经 [0, 2π)
- omega: 近地点幅角 [0, 2π)
- M: 平近点角 (仅椭圆) 或 ν: 真近点角 [0, 2π)

参考: Murray & Dermott (2.122)
"""
struct OrbitalElements
    a::Float64      # 半长径
    e::Float64      # 离心率
    i::Float64      # 倾角 [rad]
    Omega::Float64  # 升交点赤经 [rad]
    omega::Float64  # 近地点幅角 [rad]
    nu::Float64     # 真近点角 [rad]

    function OrbitalElements(a, e, i, Omega, omega, nu)
        @assert e >= 0 "Eccentricity must be non-negative"
        new(a, e, i, Omega, omega, nu)
    end
end

"""
轨道类型判别
"""
function orbit_type(e::Float64)
    if e < 1.0 - 1e-12
        return :elliptic
    elseif abs(e - 1.0) < 1e-12
        return :parabolic
    else
        return :hyperbolic
    end
end

orbit_type(el::OrbitalElements) = orbit_type(el.e)

# ============================================================
# StateVector — 笛卡尔状态向量
# ============================================================
"""
笛卡尔状态向量 (位置 + 速度)
"""
struct StateVector
    r::Vector{Float64}  # 位置 [x, y, z]
    v::Vector{Float64}  # 速度 [vx, vy, vz]
end

StateVector() = StateVector(zeros(3), zeros(3))

# ============================================================
# KeplerOrbit — 完整 Kepler 轨道定义
# ============================================================
"""
完整 Kepler 轨道: 根数 + 引力参数 + 导出量

导出量:
- T: 轨道周期 (椭圆)
- n: 平均角速度 (mean motion)
- p: 半通径
- r_p: 近心距
- r_a: 远心距
- eps: 轨道能量
- h: 角动量矢量
- h_mag: 角动量大小
"""
mutable struct KeplerOrbit
    elements::OrbitalElements
    mu::Float64         # 引力参数 GM
    T::Float64          # 轨道周期
    n::Float64          # 平均角速度
    p::Float64          # 半通径 a(1-e^2)
    r_p::Float64        # 近心距 a(1-e)
    r_a::Float64        # 远心距 a(1+e)
    epsilon::Float64    # 轨道能量 -mu/(2a)
    h::Vector{Float64}  # 角动量矢量
    h_mag::Float64      # 角动量大小
end

function KeplerOrbit(el::OrbitalElements; mu::Float64=4π^2)
    p = el.a * (1 - el.e^2)
    r_p = el.a * (1 - el.e)
    r_a = el.a * (1 + el.e)
    eps = -mu / (2 * el.a)
    T = el.e < 1.0 ? 2π * sqrt(el.a^3 / mu) : Inf
    n = el.e < 1.0 ? 2π / T : 0.0
    h_mag = sqrt(mu * p)
    h = [0.0, 0.0, h_mag]  # 轨道平面法向
    return KeplerOrbit(el, mu, T, n, p, r_p, r_a, eps, h, h_mag)
end

# ============================================================
# Perturbation — 摄动描述
# ============================================================
"""
摄动加速度 (在轨道坐标系中)

用于数值轨道传播
"""
struct PerturbationAccel
    R::Float64  # 径向分量
    S::Float64  # 横向分量 (沿速度方向, 垂直于径向)
    W::Float64  # 法向分量 (角动量方向)
end

# ============================================================
# LagrangePoint — Lagrange 点数据
# ============================================================
struct LagrangePointData
    L1::Vector{Float64}
    L2::Vector{Float64}
    L3::Vector{Float64}
    L4::Vector{Float64}
    L5::Vector{Float64}
end

# ============================================================
# TransferOrbit — 转移轨道
# ============================================================
"""
转移轨道: 包括起始轨道、目标轨道、Δv
"""
struct TransferOrbit
    name::String
    orbit_depart::OrbitalElements
    orbit_arrive::OrbitalElements
    orbit_transfer::OrbitalElements
    delta_v1::Float64       # 第一次脉冲
    delta_v2::Float64       # 第二次脉冲
    delta_v_total::Float64  # 总Δv
    transfer_time::Float64  # 转移时间
end

# ============================================================
# 常量
# ============================================================
const G_SUN = 4π^2           # 引力常数 (太阳单位) [AU³/yr²/M☉]
const AU_VAL = 1.0            # 天文单位 (太阳单位)
const YEAR_VAL = 1.0          # 年 (太阳单位)
const M_SUN_VAL = 1.0         # 太阳质量 (太阳单位)

export OrbitalElements, StateVector, KeplerOrbit, PerturbationAccel
export LagrangePointData, TransferOrbit
export orbit_type
export G_SUN, AU_VAL, YEAR_VAL, M_SUN_VAL
