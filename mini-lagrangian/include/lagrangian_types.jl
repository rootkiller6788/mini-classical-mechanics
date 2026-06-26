# lagrangian_types.jl — 拉格朗日力学核心类型接口定义
# 本文件声明所有核心数据结构和类型，作为模块的 API 合约.
# 对应：Goldstein Ch.1-2, MIT 8.012 Lectures 20-23
# L1 Definitions: core struct/typedef declarations

"""
### 广义坐标类型

GeneralizedCoords: n 个广义坐标 q = (q_1, ..., q_n)
  用于描述系统的位形 (configuration).

GeneralizedVelocities: n 个广义速度 qdot = (qdot_1, ..., qdot_n)

GeneralizedState: 系统在时刻 t 的完整力学状态
  包含时间 t, 广义坐标 q, 广义速度 qdot.
"""
struct GeneralizedCoords
    n::Int
    values::Vector{Float64}
end

struct GeneralizedVelocities
    n::Int
    values::Vector{Float64}
end

struct GeneralizedState
    t::Float64
    q::GeneralizedCoords
    qdot::GeneralizedVelocities
end

"""
### Euler-Lagrange 系统类型

EulerLagrangeSystem: 封装一个拉格朗日系统的完整数学描述.
  字段: n_dof (自由度), lagrangian L(q,qdot),
        grad_L_q (partial L / partial q), grad_L_qdot (partial L / partial qdot),
        mass_matrix M(q).

StandardLagrangian: 标准力学系统 L = T - U 的简化描述.
  字段: n_dof, M_func(q), grad_U(q).
"""
struct EulerLagrangeSystem
    n_dof::Int
    lagrangian::Function
    grad_L_q::Function
    grad_L_qdot::Function
    mass_matrix::Function
end

struct StandardLagrangian
    n_dof::Int
    M_func::Function
    grad_M_func::Function
    grad_U::Function
end
