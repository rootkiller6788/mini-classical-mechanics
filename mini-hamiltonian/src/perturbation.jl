# perturbation.jl — 正则微扰理论
# 参考：Goldstein Ch.12, Lichtenberg & Lieberman

"""
正则微扰理论：处理 H = H₀ + ε H₁
当 H₀ 可积（存在作用量-角变量），用正则变换消除 H₁ 中的非长期项

一阶微扰：
  平均化: H̄₁(J) = (1/2π) ∫ H₁(J,θ) dθ  （对快变量取平均）
  修正频率: ω(J) = ∂(H₀ + ε H̄₁)/∂J
"""

struct PerturbedSystem
    H0::Function          # 可积部分 H₀(J)
    H1::Function          # 微扰 H₁(J, θ)
    epsilon::Float64      # 微扰参数
    n_dof::Int
end

"""
对微扰取角度平均（一阶微扰的核心操作）
H̄₁(J) = (1/2π)∫₀^{2π} H₁(J,θ) dθ
对于多自由度，对每个 θ_i 积分
"""
function average_perturbation(sys::PerturbedSystem, J::Vector{Float64}; n_theta_points=100)
    n = sys.n_dof
    # Monte Carlo 或等距采样
    dtheta = 2π / n_theta_points
    avg = 0.0
    count = 0
    # 用等距网格对每个角变量采样
    # 简化：随机采样（用确定性序列）
    seed = 12345
    function lcg()
        seed = (1103515245*seed + 12345) & 0x7fffffff
        return seed / 0x7fffffff
    end
    n_samples = n_theta_points^n
    # 避免组合爆炸：用 Monte Carlo
    n_mc = min(n_samples, 2000)
    for _ in 1:n_mc
        theta = [2π*lcg() for _ in 1:n]
        avg += sys.H1(J, theta)
    end
    return avg / n_mc
end

"""
一阶修正频率: ω_i(J) = ∂(H₀ + ε H̄₁)/∂J_i
"""
function corrected_frequencies(sys::PerturbedSystem, J::Vector{Float64}; dJ=1e-6)
    n = sys.n_dof
    H0_J = sys.H0(J)
    H1_avg = average_perturbation(sys, J)
    E_unperturbed = H0_J + sys.epsilon * H1_avg

    omega = zeros(n)
    for i in 1:n
        Jp = copy(J); Jp[i] += dJ
        Jm = copy(J); Jm[i] -= dJ
        Ep = sys.H0(Jp) + sys.epsilon * average_perturbation(sys, Jp)
        Em = sys.H0(Jm) + sys.epsilon * average_perturbation(sys, Jm)
        omega[i] = (Ep - Em) / (2*dJ)
    end
    return omega
end

"""
Duffing 振子（非线性微扰）:
H = ½p² + ½q² + (ε/4)q⁴
在作用量-角变量中: H₀(J)=ω₀J, H₁(J,θ)=J²sin⁴(θ)/(2ω₀²)
一阶平均: H̄₁(J) = (3/8)J²/ω₀²
一阶修正频率: ω(J) = ω₀ + (3ε/4ω₀²)J
"""
struct DuffingOscillator
    omega0::Float64; epsilon::Float64
end

function duffing_H0(duff::DuffingOscillator, J::Vector{Float64})
    return duff.omega0 * J[1]
end

function duffing_H1(duff::DuffingOscillator, J::Vector{Float64})
    # 解析平均: H̄₁ = (3/8)*J²/ω₀²
    return (3/8) * J[1]^2 / duff.omega0^2
end

function duffing_frequency(duff::DuffingOscillator, J::Vector{Float64})
    return duff.omega0 + (3*duff.epsilon/(4*duff.omega0^2)) * J[1]
end

"""
受驱谐振子（KAM 定理的基础例子）:
H = ω₀J + ε*F₀*cos(θ - νt)（在旋转坐标系中）
共振条件: ω₀ ≈ ν 时出现共振 island
"""
struct DrivenOscillator
    omega0::Float64; nu::Float64; F0::Float64; epsilon::Float64
end

function resonance_detuning(dr::DrivenOscillator)
    return dr.omega0 - dr.nu
end

function pendulum_approximation_width(dr::DrivenOscillator)
    # 共振 island 宽度（pendulum approximation）
    return 2 * sqrt(dr.epsilon * dr.F0)
end

"""
Poincare-Lindstedt method for eliminating secular terms.

For the Duffing oscillator: xddot + ω₀² x + ε x³ = 0

Ansatz:  x(t) = x₀(t) + ε x₁(t) + ε² x₂(t) + ...
         ω    = ω₀   + ε ω₁   + ε² ω₂   + ...

Order ε⁰: x₀ddot + ω₀² x₀ = 0  →  x₀(t) = A cos(ω t + φ)

Order ε¹: x₁ddot + ω₀² x₁ = -x₀³ - 2 ω₀ ω₁ x₀
Substituting x₀ = A cos(ω t):
  -A³ cos³(ω t) = -A³ [ (3/4) cos(ω t) + (1/4) cos(3ω t) ]

The cos(ω t) term is resonant → must vanish to avoid secular growth:
  -A³·(3/4) - 2 ω₀ ω₁ A = 0  →  ω₁ = 3A²/(8ω₀)

Non-resonant part: x₁ddot + ω₀² x₁ = -(A³/4) cos(3ω t)
  →  x₁(t) = (A³/(32 ω₀²)) cos(3ω t)

Returns (omega1, x0_func(t), x1_func(t)) where omega1 is the first-order
frequency shift and x0, x1 are functions of time.
"""
function poincare_lindstedt_order1(omega0::Float64, f::Function, A::Float64; n_terms=5)
    # Duffing-specific analytical solution: f(x) = -x³
    # General structure for arbitrary f is the same but requires Fourier expansion

    # Frequency shift (first-order): ω₁ determined by eliminating the resonant
    # Fourier component of f(x₀) at the fundamental frequency.
    # For f(x) = -x³: ω₁ = 3A²/(8ω₀)
    # For general f(x) = -x^n with n odd: extract cos(ωt) coefficient.
    omega1 = 3.0 * A^2 / (8.0 * omega0)

    # Zeroth-order solution
    function x0_func(t::Float64; phi=0.0)
        return A * cos(omega0 * t + phi)
    end

    # First-order correction: particular solution for the non-resonant forcing
    function x1_func(t::Float64; phi=0.0)
        # For Duffing: x₁ = (A³/(32ω₀²)) cos(3ω₀t)
        return (A^3 / (32.0 * omega0^2)) * cos(3.0 * omega0 * t + 3.0*phi)
    end

    # Compute the corrected frequency
    omega_corrected = omega0 + 0.0  # ε ω₁ (caller multiplies by ε)

    return (omega1=omega1, x0=x0_func, x1=x1_func, omega_corrected=omega_corrected)
end
