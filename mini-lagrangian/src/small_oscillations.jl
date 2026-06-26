# small_oscillations.jl — 小振动理论
# 参考：Goldstein Ch.6, MIT 8.012 Lecture 24

using LinearAlgebra

"""
给定势能 U(q) 和质量矩阵 M(q)，找到平衡点并分析小振动

平衡条件: ∂U/∂q = 0
小振动: q = q_eq + η, L ≈ 1/2*ηdot'*M*ηdot - 1/2*η'*K*η
其中 K_{ij} = ∂²U/∂q_i∂q_j|_{q_eq}（Hessian）
简正模: det(K - ω²M) = 0 的解
"""

struct SmallOscillationSystem
    n_dof::Int
    M::Matrix{Float64}           # 在平衡点的质量矩阵
    K::Matrix{Float64}           # Hessian ∂²U/∂q_i∂q_j
    q_eq::Vector{Float64}        # 平衡位置
    frequencies::Vector{Float64} # 简正频率 ω_α
    modes::Matrix{Float64}       # 简正模向量（列），归一化
    effective_masses::Vector{Float64}
end

"""
用 Newton-Raphson 找平衡点: ∇U(q) = 0
"""
function find_equilibrium(grad_U::Function, hess_U::Function, q0::Vector{Float64}; max_iter=100, tol=1e-10)
    q = copy(q0)
    for _ in 1:max_iter
        g = grad_U(q)
        if norm(g) < tol
            return q
        end
        H = hess_U(q)
        dq = H \ g
        q = q - dq
        if norm(dq) < tol
            return q
        end
    end
    error("Newton-Raphson did not converge")
end

"""
对小振动系统求解简正模

M: 质量矩阵（在平衡点求值）
K: Hessian 矩阵（在平衡点求值）

求解广义特征值问题: K*a = ω²*M*a
"""
function solve_normal_modes(M::Matrix{Float64}, K::Matrix{Float64}, q_eq::Vector{Float64})
    # 广义特征值: K*a = λ*M*a
    eigen_decomp = eigen(K, Symmetric(M))
    omega_sq = eigen_decomp.values
    modes_raw = eigen_decomp.vectors

    n = length(omega_sq)
    # 过滤非正频率（不稳定平衡点）
    positive_mask = omega_sq .> 1e-12
    n_pos = sum(positive_mask)

    if n_pos == 0
        @warn "All modes unstable: no positive eigenvalues" omega_sq
        return SmallOscillationSystem(n, M, K, q_eq, Float64[], Matrix{Float64}(undef,n,0), Float64[])
    end

    omegas = sqrt.(omega_sq[positive_mask])
    modes = modes_raw[:, positive_mask]

    # 归一化: a_α'*M*a_β = δ_{αβ}
    for alpha in 1:size(modes, 2)
        norm_factor = sqrt(dot(modes[:,alpha], M * modes[:,alpha]))
        modes[:,alpha] /= norm_factor
    end

    # 有效质量
    eff_masses = [dot(modes[:,a], M * modes[:,a]) for a in 1:size(modes, 2)]

    return SmallOscillationSystem(n, M, K, q_eq, omegas, modes, eff_masses)
end

"""
数值 Hessian: ∂²U/∂q_i∂q_j
"""
function numerical_hessian(U::Function, q::Vector{Float64}; h=1e-5)
    n = length(q)
    H = zeros(n, n)
    U0 = U(q)
    for i in 1:n
        qp = copy(q); qp[i] += h; Up = U(qp)
        qm = copy(q); qm[i] -= h; Um = U(qm)
        H[i,i] = (Up - 2*U0 + Um) / h^2
        for j in i+1:n
            qpp = copy(q); qpp[i] += h; qpp[j] += h; Upp = U(qpp)
            qpm = copy(q); qpm[i] += h; qpm[j] -= h; Upm = U(qpm)
            qmp = copy(q); qmp[i] -= h; qmp[j] += h; Ump = U(qmp)
            qmm = copy(q); qmm[i] -= h; qmm[j] -= h; Umm = U(qmm)
            H[i,j] = (Upp - Upm - Ump + Umm) / (4*h^2)
            H[j,i] = H[i,j]
        end
    end
    return H
end

"""
数值梯度: ∂U/∂q_i
"""
function numerical_gradient(U::Function, q::Vector{Float64}; h=1e-6)
    n = length(q)
    grad = zeros(n)
    for i in 1:n
        qp = copy(q); qp[i] += h
        qm = copy(q); qm[i] -= h
        grad[i] = (U(qp) - U(qm)) / (2*h)
    end
    return grad
end

"""
从简正坐标重建物理坐标
q(t) = q_eq + Σ_α c_α*cos(ω_α*t + φ_α)*a_α
"""
function normal_mode_solution(sys::SmallOscillationSystem, amplitudes::Vector{Float64}, phases::Vector{Float64}, t::Float64)
    n = sys.n_dof
    q = copy(sys.q_eq)
    n_active = length(sys.frequencies)
    for alpha in 1:n_active
        q += amplitudes[alpha] * cos(sys.frequencies[alpha]*t + phases[alpha]) * sys.modes[:,alpha]
    end
    return q
end

"""
多自由度耦合谐振子的完整分析流水线
"""
function analyze_small_oscillations(M_func::Function, U_func::Function, q0::Vector{Float64})
    # 找平衡点
    grad_U(q) = numerical_gradient(U_func, q)
    hess_U(q) = numerical_hessian(U_func, q)
    q_eq = find_equilibrium(grad_U, hess_U, q0)
    println("Equilibrium found at: $q_eq")

    # Hessian + 质量矩阵
    K = hess_U(q_eq)
    M = M_func(q_eq)
    println("M = $M")
    println("K = $K")

    # 简正模
    sys = solve_normal_modes(M, K, q_eq)
    println("\nNormal modes:")
    for alpha in 1:length(sys.frequencies)
        println("  ω_$alpha = $(round(sys.frequencies[alpha], digits=4)) rad/s")
        println("  mode_$alpha = $(round.(sys.modes[:,alpha], digits=4))")
        println("  effective mass = $(round(sys.effective_masses[alpha], digits=4))")
    end

    return sys
end

"""
受迫小振动: M η̈ + C η̇ + K η = F(t)

求解稳态频率响应 (假设 F(t)=F0 cos(ωt)).
返回振幅向量 |η(ω)|.
"""
function forced_response(M::Matrix{Float64}, K::Matrix{Float64}, C::Matrix{Float64}, F0::Vector{Float64}, omega::Float64)
    n = size(M,1)
    Z = -omega^2*M + 1im*omega*C + K
    eta_complex = Z \ F0
    return abs.(eta_complex), angle.(eta_complex)
end

"""
扫频: 计算多自由度系统的频率响应函数 (FRF).
"""
function frequency_sweep(M::Matrix{Float64}, K::Matrix{Float64}, C::Matrix{Float64}, F0::Vector{Float64}, omega_range::AbstractRange)
    n = size(M,1); m = length(omega_range)
    amplitudes = zeros(n, m)
    for (j, omega) in enumerate(omega_range)
        amp, _ = forced_response(M, K, C, F0, omega)
        amplitudes[:, j] = amp
    end
    return collect(omega_range), amplitudes
end

"""
阻尼比: ζ = c/(2√(mk)) = c/(2mω_n).
"""
damping_ratio(c::Float64, m::Float64, omega_n::Float64) = c/(2*m*omega_n)

"""
品质因子: Q = 1/(2ζ) = ω_n m / c.
"""
quality_factor(c::Float64, m::Float64, omega_n::Float64) = omega_n*m/c

"""
共振频率 (有阻尼): ω_r = ω_n √(1-2ζ²).
"""
resonant_frequency(omega_n::Float64, zeta::Float64) = omega_n*sqrt(1-2*zeta^2)

"""
Rayleigh 阻尼: C = α M + β K.
"""
function rayleigh_damping(M::Matrix{Float64}, K::Matrix{Float64}, alpha::Float64, beta::Float64)
    return alpha*M + beta*K
end
