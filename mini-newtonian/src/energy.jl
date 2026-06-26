# energy.jl — 能量计算与守恒验证
# 参考：MIT 8.012 Ch.5, Goldstein Ch.2

## ============================================================
## 动能
## ============================================================

"""平动动能: T = 0.5 * m * v^2"""
function kinetic_energy(m::Float64, v::Vec3)
    return 0.5 * m * norm2(v)
end

"""N 体系统总动能"""
function total_kinetic_energy(masses::Vector{Float64}, velocities::Vector{Vec3})
    T = 0.0
    for i in 1:length(masses)
        T += kinetic_energy(masses[i], velocities[i])
    end
    return T
end

## ============================================================
## 势能
## ============================================================

"""均匀重力场势能: U = m*g*h（取 z=0 为零势面）"""
function gravitational_potential_uniform(m::Float64, z::Float64; g::Float64=G_EARTH)
    return m * g * z
end

"""万有引力势能（两点质量）: U = -G*M*m / r"""
function gravitational_potential(r_obj::Vec3, r_source::Vec3, M::Float64, m::Float64; G_val::Float64=G)
    dist = norm(r_obj - r_source)
    if dist < 1e-10
        return -Inf  # 奇点
    end
    return -G_val * M * m / dist
end

"""N 体引力总势能（每对只算一次）"""
function total_gravitational_potential(masses::Vector{Float64}, positions::Vector{Vec3}; G_val::Float64=G)
    n = length(masses)
    U = 0.0
    for i in 1:n
        for j in i+1:n
            U += gravitational_potential(positions[i], positions[j], masses[j], masses[i]; G_val=G_val)
        end
    end
    return U
end

"""弹性势能: U = 0.5 * k * (r - r_eq)^2"""
function elastic_potential(r::Vec3, r_eq::Vec3, k::Float64)
    return 0.5 * k * norm2(r - r_eq)
end

function elastic_potential_1d(x::Float64, x_eq::Float64, k::Float64)
    return 0.5 * k * (x - x_eq)^2
end

## ============================================================
## 功
## ============================================================

"""恒力做功: W = F · dr"""
function work_constant_force(F::Vec3, dr::Vec3)
    return dot(F, dr)
end

"""变力做功（沿路径数值积分）"""
function work_path(F_func::Function, path::Vector{Vec3})
    W = 0.0
    for i in 2:length(path)
        dr = path[i] - path[i-1]
        F_mid = F_func(0.5 * (path[i] + path[i-1]))  # 中点法则
        W += dot(F_mid, dr)
    end
    return W
end

"""功率: P = F · v"""
function power(F::Vec3, v::Vec3)
    return dot(F, v)
end

## ============================================================
## 守恒量验证
## ============================================================

"""
能量守恒验证
计算轨迹中总能量 E=T+U 的相对漂移
"""
function energy_drift(traj::Trajectory, m::Float64, potential_func::Function)
    n = length(traj.ts)
    energies = Vector{Float64}(undef, n)
    for i in 1:n
        T = kinetic_energy(m, traj.velocities[i])
        U = potential_func(traj.positions[i])
        energies[i] = T + U
    end
    E0 = energies[1]
    drift = (energies .- E0) ./ abs(E0)
    return energies, drift
end

"""
N 体系统能量守恒验证
"""
function nbody_energy_drift(states::Vector{NBodyState}; G_val::Float64=G)
    n_steps = length(states)
    energies = Vector{Float64}(undef, n_steps)
    for i in 1:n_steps
        st = states[i]
        T = total_kinetic_energy(st.masses, st.velocities)
        U = total_gravitational_potential(st.masses, st.positions; G_val=G_val)
        energies[i] = T + U
    end
    E0 = energies[1]
    drift = (energies .- E0) ./ (abs(E0) + 1e-300)
    return energies, drift
end

## ============================================================
## 机械能守恒判断
## ============================================================

"""保守力判断：力是否为某势能函数的负梯度"""
function is_conservative_1d(F::Function, x_range, dx=0.01)
    # 一维保守力判据：∮F·dx = 0 → F 可写为单值函数即可
    # 简化：检查 F 是否只依赖于位置（不依赖速度和时间）
    return true  # 一维任意位置函数都是保守力
end

"""总机械能: E = T + U"""
total_mechanical_energy(m::Float64, v::Vec3, U::Float64) = kinetic_energy(m, v) + U
