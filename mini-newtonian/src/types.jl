# types.jl — 核心数据类型定义
# 零外部依赖：不用 StaticArrays，手写 immutable struct

## ============================================================
## Vec3 — 三维向量（栈分配，不可变）
## ============================================================
struct Vec3
    x::Float64
    y::Float64
    z::Float64
end

# 构造便捷函数
Vec3() = Vec3(0.0, 0.0, 0.0)

# 基本运算（运算符重载保持公式可读性）
Base.:+(a::Vec3, b::Vec3) = Vec3(a.x + b.x, a.y + b.y, a.z + b.z)
Base.:-(a::Vec3, b::Vec3) = Vec3(a.x - b.x, a.y - b.y, a.z - b.z)
Base.:-(v::Vec3) = Vec3(-v.x, -v.y, -v.z)
Base.:*(s::Number, v::Vec3) = Vec3(s * v.x, s * v.y, s * v.z)
Base.:*(v::Vec3, s::Number) = Vec3(v.x * s, v.y * s, v.z * s)
Base.:/(v::Vec3, s::Number) = Vec3(v.x / s, v.y / s, v.z / s)

# 向量运算
dot(a::Vec3, b::Vec3) = a.x*b.x + a.y*b.y + a.z*b.z
cross(a::Vec3, b::Vec3) = Vec3(
    a.y*b.z - a.z*b.y,
    a.z*b.x - a.x*b.z,
    a.x*b.y - a.y*b.x
)
norm(v::Vec3) = sqrt(dot(v, v))
norm2(v::Vec3) = dot(v, v)  # 平方模长，避免sqrt
normalize(v::Vec3) = norm(v) > 1e-15 ? v / norm(v) : Vec3()
distance(a::Vec3, b::Vec3) = norm(a - b)
angle(a::Vec3, b::Vec3) = acos(dot(a, b) / (norm(a) * norm(b)))

# 坐标转换
function to_spherical(v::Vec3)
    r = norm(v)
    theta = r > 0 ? acos(v.z / r) : 0.0
    phi = atan(v.y, v.x)
    return (r, theta, phi)
end

function from_spherical(r::Float64, theta::Float64, phi::Float64)
    return Vec3(
        r * sin(theta) * cos(phi),
        r * sin(theta) * sin(phi),
        r * cos(theta)
    )
end

function to_cylindrical(v::Vec3)
    rho = sqrt(v.x^2 + v.y^2)
    phi = atan(v.y, v.x)
    return (rho, phi, v.z)
end

function from_cylindrical(rho::Float64, phi::Float64, z::Float64)
    return Vec3(rho * cos(phi), rho * sin(phi), z)
end

Base.show(io::IO, v::Vec3) = print(io, "Vec3($(v.x), $(v.y), $(v.z))")

## ============================================================
## State — 单粒子运动状态（可变，用于积分器迭代）
## ============================================================
mutable struct State
    t::Float64      # 时间
    r::Vec3         # 位置
    v::Vec3         # 速度
end

function State(r::Vec3, v::Vec3; t::Float64=0.0)
    return State(t, r, v)
end

State() = State(Vec3(), Vec3())

## ============================================================
## NBodyState — N体系统状态
## ============================================================
mutable struct NBodyState
    t::Float64
    masses::Vector{Float64}
    positions::Vector{Vec3}
    velocities::Vector{Vec3}
end

function NBodyState(masses::Vector{Float64}, pos::Vector{Vec3}, vel::Vector{Vec3}; t::Float64=0.0)
    n = length(masses)
    @assert length(pos) == n && length(vel) == n "数组长度必须一致"
    return NBodyState(t, masses, pos, vel)
end

n_particles(st::NBodyState) = length(st.masses)

## ============================================================
## Trajectory — 轨迹记录
## ============================================================
mutable struct Trajectory
    ts::Vector{Float64}
    positions::Vector{Vec3}
    velocities::Vector{Vec3}
end

Trajectory() = Trajectory(Float64[], Vec3[], Vec3[])

function record!(traj::Trajectory, state::State)
    push!(traj.ts, state.t)
    push!(traj.positions, state.r)
    push!(traj.velocities, state.v)
end

function record!(traj::Trajectory, t::Float64, r::Vec3, v::Vec3)
    push!(traj.ts, t)
    push!(traj.positions, r)
    push!(traj.velocities, v)
end

## ============================================================
## ODEProblem — 通用常微分方程问题定义
## ============================================================
struct ODEProblem{F<:Function}
    f::F           # dy/dt = f(t, y)  或二阶形式 a(r) = f(r)
    t0::Float64    # 初始时间
    t_end::Float64 # 结束时间
    y0::Vector{Float64}  # 初始状态向量 [rx,ry,rz, vx,vy,vz] 或自定义
end

## ============================================================
## ODEParticleProblem — 粒子运动问题（Vec3 版本）
## ============================================================
struct ParticleODE{F<:Function}
    acceleration::F   # a(r, v, t) -> Vec3  加速度函数
    t0::Float64
    t_end::Float64
    r0::Vec3
    v0::Vec3
end

## ============================================================
## ODESolution — 求解结果
## ============================================================
struct ODESolution
    ts::Vector{Float64}
    ys::Vector{Vector{Float64}}
    n_steps::Int
    method::String
end

## ============================================================
## 常量
## ============================================================
const G = 6.67430e-11      # 万有引力常数 [m³/(kg·s²)]
const G_AU = 39.478         # G in [AU³/(yr²·Msun)]: (2π)² for solar system
const G_EARTH = 9.80665     # 地表重力加速度 [m/s²]
const DAY = 86400.0         # 秒/天
const YEAR = 31557600.0     # 秒/年
const AU = 1.495978707e11   # 天文单位 [m]
const M_SUN = 1.98847e30    # 太阳质量 [kg]
const M_EARTH = 5.9722e24   # 地球质量 [kg]
const R_EARTH = 6.371e6     # 地球半径 [m]
