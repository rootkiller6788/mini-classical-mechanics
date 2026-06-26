# types.jl — 混沌动力学类型定义
# 参考: Strogatz, Goldstein Ch.11

"动力学系统: dx/dt = f(x, params...)"
struct DynamicalSystem{F<:Function}
    f::F
    dim::Int
    name::String
    params::Dict{Symbol, Float64}
end

"分岔数据: 参数值和对应的稳态值"
mutable struct BifurcationData
    param_name::String
    param_vals::Vector{Float64}
    steady_vals::Vector{Vector{Float64}}
    periods::Vector{Int}
end

"Lyapunov 谱结果"
struct LyapunovSpectrum
    exponents::Vector{Float64}
    sum_pos::Float64
    kaplan_yorke_dim::Float64
end

export DynamicalSystem, BifurcationData, LyapunovSpectrum
