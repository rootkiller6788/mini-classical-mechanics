# failure.jl — 强度理论: von Mises, Tresca, Mohr-Coulomb, Drucker-Prager, 疲劳
# 参考: Juvinall "Engineering Considerations of Stress", Dowling "Mechanical Behavior"

using LinearAlgebra

"von Mises 等效应力"
function von_mises(sigma::StressTensor)
    s=principal_stresses(sigma); sqrt(0.5*((s[1]-s[2])^2+(s[2]-s[3])^2+(s[3]-s[1])^2))
end

"von Mises from invariants: σ_vm = √(3J2)"
von_mises_J2(J2::Float64) = sqrt(3*J2)

"Tresca 等效应力"
function tresca(sigma::StressTensor)
    s=principal_stresses(sigma); max(abs(s[1]-s[2]),abs(s[2]-s[3]),abs(s[3]-s[1]))
end

"最大主应力准则 (Rankine)"
max_principal(sigma::StressTensor) = maximum(principal_stresses(sigma))

"最大主应变准则 (Saint-Venant)"
function max_principal_strain(sigma::StressTensor, E::Float64, nu::Float64)
    s=principal_stresses(sigma)
    max(s[1]-nu*(s[2]+s[3]), s[2]-nu*(s[1]+s[3]), s[3]-nu*(s[1]+s[2]))/E
end

"Mohr-Coulomb 剪切强度: τ=c+σtanφ"
mohr_coulomb_shear(c::Float64, sigma_n::Float64, phi::Float64) = c+sigma_n*tan(phi)

"Mohr-Coulomb 主应力形式: σ₁/σ_ut + σ₃/σ_uc = 1"
function mohr_coulomb_envelope(sigma::StressTensor, sigma_ut::Float64, sigma_uc::Float64)
    s=principal_stresses(sigma); s[1]/sigma_ut - s[3]/sigma_uc
end

"Drucker-Prager: √J2 + αI1 = k"
function drucker_prager(sigma::StressTensor, alpha::Float64)
    inv=stress_invariants(sigma); sqrt(inv.J2)+alpha*inv.I1
end

"安全系数 (静载): n = σ_yield / σ_vm"
safety_factor_static(yield_stress::Float64, sigma_vm::Float64) = yield_stress/sigma_vm

"安全系数 (Mohr-Coulomb): n = (c+σ_n tanφ)/τ"
function safety_factor_mc(c::Float64, sigma_n::Float64, phi::Float64, tau::Float64)
    (c+sigma_n*tan(phi))/max(tau,1e-15)
end

# --- 疲劳 ---
"Basquin S-N: σ_a = σ_f'(2Nf)^b"
basquin_sn_curve(Nf::Float64, sigma_f_prime::Float64, b::Float64) = sigma_f_prime*(2*Nf)^b

"疲劳寿命: Nf = ½(σ_a/σ_f')^(1/b)"
fatigue_life(sigma_a::Float64, sigma_f_prime::Float64, b::Float64) = 0.5*(sigma_a/sigma_f_prime)^(1/b)

"Coffin-Manson (应变-寿命): ε_a = (σ_f'/E)(2Nf)^b + ε_f'(2Nf)^c"
function strain_life(Nf::Float64, E::Float64, sigma_f_prime::Float64, b::Float64, eps_f_prime::Float64, c::Float64)
    sigma_f_prime/E*(2*Nf)^b + eps_f_prime*(2*Nf)^c
end

"Goodman 平均应力修正: σ_a/σ_ar + σ_m/σ_uts = 1"
goodman(sigma_a::Float64, sigma_m::Float64, sigma_uts::Float64) = sigma_a/(1-sigma_m/sigma_uts)

"Gerber 抛物线修正: σ_a/σ_ar + (σ_m/σ_uts)² = 1"
gerber(sigma_a::Float64, sigma_m::Float64, sigma_uts::Float64) = sigma_a/(1-(sigma_m/sigma_uts)^2)

"Soderberg 修正: σ_a/σ_ar + σ_m/σ_y = 1"
soderberg(sigma_a::Float64, sigma_m::Float64, sigma_y::Float64) = sigma_a/(1-sigma_m/sigma_y)

"Morrow 修正: σ_ar = (σ_f'-σ_m)(2Nf)^b"
function morrow_equivalent(sigma_a::Float64, sigma_m::Float64, sigma_f_prime::Float64, b::Float64)
    (sigma_f_prime-sigma_m)*(sigma_a/sigma_f_prime)^(1/b)
end

"Palmgren-Miner 累积损伤: D = Σ n_i/N_fi"
function miner_damage(cycles::Vector{Float64}, lives::Vector{Float64})
    sum(cycles./lives)
end

"雨流计数 (简化: 仅计数极值)"
function simple_cycle_count(signal::Vector{Float64})
    cycles=Float64[]; means=Float64[]
    for i in 2:length(signal)
        push!(cycles,abs(signal[i]-signal[i-1]))
        push!(means,(signal[i]+signal[i-1])/2)
    end
    cycles, means
end

"Neuber 缺口修正: K_t² = K_σ K_ε"
neuber_notch(Kt::Float64) = Kt^2, Kt, Kt  # K_t^2 = K_σ * K_ε

export von_mises, von_mises_J2, tresca, max_principal, max_principal_strain
export mohr_coulomb_shear, mohr_coulomb_envelope, drucker_prager
export safety_factor_static, safety_factor_mc
export basquin_sn_curve, fatigue_life, strain_life
export goodman, gerber, soderberg, morrow_equivalent
export miner_damage, simple_cycle_count, neuber_notch
