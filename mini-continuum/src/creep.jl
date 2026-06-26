# creep.jl — 蠕变力学: 幂律蠕变、扩散蠕变、Larson-Miller、松弛
# 参考: Kassner "Fundamentals of Creep", Nabarro & de Villiers

"幂律蠕变 (Norton): d(eps_c)/dt = A sigma^n exp(-Q/RT)"
function norton_creep_rate(A, sigma, n, Q, R, T)
    A * sigma^n * exp(-Q/(R*T))
end
"二次蠕变应变: eps = A sigma^n t exp(-Q/RT)"
function norton_creep_strain(A, sigma, n, Q, R, T, t)
    norton_creep_rate(A, sigma, n, Q, R, T) * t
end
"扩散蠕变 (Nabarro-Herring, 体扩散): eps_dot = B sigma Omega D_v/(kT d^2)"
function nabarro_herring_rate(B, sigma, Omega, Dv, k, T, d)
    B * sigma * Omega * Dv / (k * T * d^2)
end
"Coble蠕变 (晶界扩散): eps_dot = B sigma Omega delta D_gb/(kT d^3)"
function coble_creep_rate(B, sigma, Omega, delta, Dgb, k, T, d)
    B * sigma * Omega * delta * Dgb / (k * T * d^3)
end
"扩散系数 (Arrhenius): D = D0 exp(-Q/RT)"
function arrhenius_diffusion(D0, Q, R, T)
    D0 * exp(-Q/(R*T))
end
"Larson-Miller参数: LMP = T*(C + log10(t_r)), t_r rupture time [h]"
function larson_miller_parameter(T, C, t_r_hours)
    T * (C + log10(t_r_hours))
end
"Larson-Miller 寿命反算: t_r = 10^(LMP/T - C)"
function larson_miller_life(LMP, T, C)
    10^(LMP/T - C)
end
"Monkman-Grant关系: d(eps)/dt * t_r^m = C_MG"
function monkman_grant(eps_dot_ss, t_r, m, C_MG)
    eps_dot_ss * t_r^m - C_MG  # residual
end
"Omega蠕变模型 (Prasad): eps = eps_0 exp(Omega*eps)"
function omega_creep(eps0, Omega)
    eps0 * exp(Omega * eps0)
end
"应力松弛 (Maxwell): sigma(t) = sigma_0 exp(-t*E/eta)"
function stress_relaxation_maxwell(sigma0, E, eta, t)
    sigma0 * exp(-E*t/eta)
end
"应力松弛 幂律: d(sigma)/dt = -k sigma^n"
function stress_relaxation_powerlaw(sigma0, k, n, t)
    if n == 1.0; sigma0 * exp(-k*t)
    else; (sigma0^(1-n) - k*(1-n)*t)^(1/(1-n)); end
end
"蠕变柔量: J(t) = eps(t)/sigma_0"
creep_compliance_c(eps_t, sigma0) = eps_t/sigma0
"蠕变回复 (Kelvin): eps(t) = eps_0 exp(-t/tau)"
function creep_recovery(eps0, tau, t)
    eps0 * exp(-t/tau)
end
"稳态蠕变速率 (Sherby-Dorn): d(eps)/dt = A D (sigma/E)^n"
function sherby_dorn_rate(A, D, sigma, E, n)
    A * D * (sigma/E)^n
end
"蠕变断裂时间 (Kachanov损伤): t_r = 1/(A*(nu+1)*sigma^nu)"
function kachanov_rupture_time(A, nu, sigma)
    1.0 / (A * (nu + 1) * sigma^nu)
end
"Kachanov损伤演变: d(omega)/dt = A sigma^nu / (1-omega)^m"
function kachanov_damage_rate(A, sigma, nu, omega, m)
    A * sigma^nu / ((1 - omega)^m)
end
"蠕变疲劳交互 (线性累积): D_total = D_creep + D_fatigue"
function creep_fatigue_damage(D_creep, D_fatigue)
    D_creep + D_fatigue
end
"温度加速因子: AF = exp(Q/R * (1/T_use - 1/T_test))"
function thermal_acceleration_factor(Q, R, T_use, T_test)
    exp(Q/R * (1/T_use - 1/T_test))
end
"晶粒尺寸效应 (Hall-Petch for creep): sigma ~ d^(-1/2)"
function grain_size_strength(sigma0, k_hp, d)
    sigma0 + k_hp / sqrt(d)
end
"Ashby变形机制图: 归一化应力 vs 归一化温度"
function ashby_normalized_stress(sigma, G)
    sigma / G
end

export norton_creep_rate, norton_creep_strain
export nabarro_herring_rate, coble_creep_rate, arrhenius_diffusion
export larson_miller_parameter, larson_miller_life, monkman_grant
export omega_creep, stress_relaxation_maxwell, stress_relaxation_powerlaw
export creep_compliance_c, creep_recovery, sherby_dorn_rate
export kachanov_rupture_time, kachanov_damage_rate, creep_fatigue_damage
export thermal_acceleration_factor, grain_size_strength, ashby_normalized_stress
