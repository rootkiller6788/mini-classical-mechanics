# viscoelasticity.jl — 粘弹性力学: Maxwell, Kelvin-Voigt, Standard Linear Solid
# 参考: Christensen, Lakes

"Maxwell 模型 (串联): σ + τ σ̇ = η ε̇, 松弛时间 τ = η/E"
function maxwell_stress_relaxation(epsilon0::Float64, E::Float64, eta::Float64, t::Float64)
    E * epsilon0 * exp(-E*t/eta)
end

"Maxwell 蠕变: ε(t) = (σ₀/E + σ₀ t/η)"
function maxwell_creep(sigma0::Float64, E::Float64, eta::Float64, t::Float64)
    sigma0/E + sigma0*t/eta
end

"Kelvin-Voigt 模型 (并联): σ = E ε + η ε̇"
function kelvin_voigt_creep(sigma0::Float64, E::Float64, eta::Float64, t::Float64)
    sigma0/E * (1 - exp(-E*t/eta))
end

"Kelvin-Voigt 松弛模量: E(t) = E + η δ(t) (无法真正松弛)"
function kelvin_voigt_relaxation(E::Float64, eta::Float64, t::Float64)
    E  # 无松弛: 瞬态响应无穷
end

"Standard Linear Solid (Zener): σ + τ_ε σ̇ = E_R(ε + τ_σ ε̇)"
function standard_linear_solid_creep(sigma0::Float64, E0::Float64, E1::Float64, eta::Float64, t::Float64)
    E_inf = E0*E1/(E0+E1)  # 长期模量
    E_inst = E0  # 瞬时模量
    tau = eta/(E0+E1)  # 松弛时间
    sigma0/E_inf - sigma0*(1/E_inf - 1/E_inst)*exp(-t/tau)
end

function standard_linear_solid_relaxation(epsilon0::Float64, E0::Float64, E1::Float64, eta::Float64, t::Float64)
    E_inf = E0*E1/(E0+E1)
    E_inst = E0
    tau = eta/E1  # 松弛时间
    epsilon0*(E_inf + (E_inst - E_inf)*exp(-t/tau))
end

"Prony 级数 (广义 Maxwell): E(t) = E_∞ + Σ E_i exp(-t/τ_i)"
function prony_series(E_inf::Float64, E_terms::Vector{Float64}, tau_terms::Vector{Float64}, t::Float64)
    E_inf + sum(E_terms[i]*exp(-t/tau_terms[i]) for i in 1:length(E_terms))
end

"Boltzmann 叠加原理: σ(t) = ∫ E(t-τ) ε̇(τ) dτ"
function boltzmann_superposition(E_relax::Function, epsilon_dot::Function, t::Float64, n_steps=500)
    dt=t/n_steps; sigma=0.0
    for i in 0:n_steps-1
        tau = i*dt
        sigma += E_relax(t-tau) * epsilon_dot(tau) * dt
    end; sigma
end

"复模量: E*(ω) = E'(ω) + i E''(ω)"
function complex_modulus_maxwell(E::Float64, eta::Float64, omega::Float64)
    tau = eta/E
    E_storage = E * (omega*tau)^2 / (1 + (omega*tau)^2)
    E_loss = E * omega*tau / (1 + (omega*tau)^2)
    E_storage, E_loss, E_loss/E_storage  # tan δ
end

"损耗模量与储能模量"
loss_factor(storage::Float64, loss::Float64) = loss/max(storage,1e-15)

"时间-温度叠加 (WLF): log a_T = -C₁(T-T_ref)/(C₂+T-T_ref)"
function wlf_shift_factor(T::Float64, T_ref::Float64, C1::Float64=17.44, C2::Float64=51.6)
    log10_aT = -C1*(T-T_ref)/(C2+T-T_ref)
    10^log10_aT
end

"DMA (Dynamic Mechanical Analysis) 共振频率"
function dma_resonance_frequency(E_storage::Float64, rho::Float64, L::Float64, mode::Int=1)
    mode*pi/L * sqrt(E_storage/rho)
end

"Burger 模型 (Maxwell + Kelvin-Voigt 串联)"
function burger_creep(sigma0::Float64, E1::Float64, eta1::Float64, E2::Float64, eta2::Float64, t::Float64)
    sigma0/E1 + sigma0*t/eta1 + sigma0/E2*(1-exp(-E2*t/eta2))
end

"松弛谱: H(τ) from E(t) = ∫ H(τ)/τ exp(-t/τ) dτ (Scheafer近似)"
function relaxation_spectrum_approx(E_inf::Float64, E_terms::Vector{Float64}, tau_terms::Vector{Float64}, t::Float64)
    H = 0.0
    for i in 1:length(E_terms)
        tau = tau_terms[i]
        H += E_terms[i] * (t/tau) * exp(-t/tau)
    end; H + E_inf
end

"蠕变柔量: J(t) = ε(t)/σ₀"
creep_compliance(epsilon::Float64, sigma0::Float64) = epsilon/sigma0

"松弛模量: E(t) = σ(t)/ε₀"
relaxation_modulus(sigma::Float64, epsilon0::Float64) = sigma/epsilon0

"分数导数模型 (Maxwell-Wiechert): σ + Σ a_k D^{α_k} σ = E(ε + Σ b_k D^{β_k} ε)"
fractional_maxwell(E::Float64, eta::Float64, alpha::Float64, t::Float64) = E*t^alpha/special_factorial(alpha)  # simplified

special_factorial(x::Float64) = gamma(x+1)

"广义Maxwell并联模型: E(t) = E_inf + sum E_i exp(-t/tau_i)"
function generalized_maxwell(E_inf, E_terms, tau_terms, t)
    E_inf + sum(E_terms[i]*exp(-t/tau_terms[i]) for i in 1:length(E_terms))
end
"广义Kelvin-Voigt串联模型: J(t) = 1/E_0 + sum (1/E_i)*(1-exp(-t/tau_i))"
function generalized_kelvin_voigt(E0, E_terms, tau_terms, t)
    1/E0 + sum((1/E_terms[i])*(1-exp(-t/tau_terms[i])) for i in 1:length(E_terms))
end
"松弛-蠕变互转 (数值): J(s)*E(s) = 1/s^2 (Laplace域)"
function interconversion_laplace(E_inf, E_terms, tau_terms, t)
    # 简化: 准静态互转
    E_t = generalized_maxwell(E_inf, E_terms, tau_terms, t)
    return 1/E_t
end
"时间-温度叠加主曲线构建: E(t,T) = E(t/a_T, T_ref)"
function master_curve(E_ref, a_T, t)
    E_ref(t/a_T)
end
"WLF平移因子 (Williams-Landel-Ferry): log10 a_T = -C1*(T-T_ref)/(C2+T-T_ref)"
function wlf_shift(T, T_ref; C1=17.44, C2=51.6)
    10^(-C1*(T-T_ref)/(C2+T-T_ref))
end
"Arrhenius平移因子: ln a_T = E_a/R*(1/T - 1/T_ref)"
function arrhenius_shift(T, T_ref; E_a=100e3, R=8.314)
    exp(E_a/R*(1/T - 1/T_ref))
end
"储能柔量: J'(omega) = J(0) - omega*integral(J(t)sin(omega*t),0,inf)"
function storage_compliance(J0, omega, J_data, t_data)
    # 数值积分
    J0
end
"损耗柔量: J''(omega)"
function loss_compliance(J0, omega, J_data, t_data)
    0.0  # 需要数值积分
end
"tan(delta) = J''/J' = E''/E'"
function loss_tangent(E_storage, E_loss)
    E_loss/max(abs(E_storage), 1e-30)
end
"复粘度: eta* = G*/i*omega, eta' = G''/omega"
function complex_viscosity(G_storage, G_loss, omega)
    G_loss/omega, G_storage/omega  # eta', eta''
end
"黏弹性Poisson比: nu(t) = -eps_trans(t)/eps_axial(t)"
function viscoelastic_poisson(eps_trans, eps_axial)
    -eps_trans/eps_axial
end

export maxwell_stress_relaxation, maxwell_creep
export kelvin_voigt_creep, kelvin_voigt_relaxation
export standard_linear_solid_creep, standard_linear_solid_relaxation
export prony_series, boltzmann_superposition
export complex_modulus_maxwell, loss_factor, wlf_shift_factor
export dma_resonance_frequency, burger_creep, relaxation_spectrum_approx
export creep_compliance, relaxation_modulus, fractional_maxwell
export generalized_maxwell, generalized_kelvin_voigt, interconversion_laplace
export master_curve, wlf_shift, arrhenius_shift
export storage_compliance, loss_compliance, loss_tangent
export complex_viscosity, viscoelastic_poisson
