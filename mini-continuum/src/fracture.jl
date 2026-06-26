# fracture.jl — 断裂力学: 应力强度因子, 能量释放率, J积分
# 参考: Anderson "Fracture Mechanics", Broberg

"应力强度因子 (mode I, 无限板中心裂纹): K_I = σ √(πa)"
function stress_intensity_factor_mode1(sigma::Float64, a::Float64, beta::Float64=1.0)
    beta * sigma * sqrt(pi*a)
end

"mode II: K_II = τ √(πa)"
function stress_intensity_factor_mode2(tau::Float64, a::Float64, beta::Float64=1.0)
    beta * tau * sqrt(pi*a)
end

"能量释放率 (平面应力): G = K²/E"
function energy_release_rate(K::Float64, E::Float64, nu::Float64=0.0, plane_stress::Bool=true)
    if plane_stress; K^2/E; else; K^2*(1-nu^2)/E; end
end

"断裂韧性判据: K_I ≥ K_IC → 失稳扩展"
function fracture_criterion(K_I::Float64, K_IC::Float64)
    K_I >= K_IC
end

"Paris 疲劳裂纹扩展: da/dN = C (ΔK)^m"
function paris_law(delta_K::Float64, C::Float64, m::Float64)
    C * delta_K^m
end

"塑性区尺寸 (Irwin): r_p = (1/2π)(K/σ_y)² (平面应力)"
function plastic_zone_size(K::Float64, sigma_y::Float64)
    0.5/pi * (K/sigma_y)^2
end

"J 积分 (与路径无关, 弹塑性断裂)"
function j_integral_estimate(K::Float64, E::Float64, nu::Float64=0.0)
    K^2 / E  # 弹性: J=G
end

"CTOD (Crack Tip Opening Displacement): δ = 4K²/(π E σ_y)"
function ctod_estimate(K::Float64, E::Float64, sigma_y::Float64)
    4*K^2/(pi*E*sigma_y)
end

"Weibull 统计断裂: P_f = 1 - exp(-(σ/σ₀)^m)"
function weibull_failure_probability(sigma::Float64, sigma0::Float64, m::Float64)
    1.0 - exp(-(sigma/sigma0)^m)
end

"Weibull 尺度参数 (最弱链模型)"
function weibull_scale_parameter(sigma0::Float64, V::Float64, V0::Float64, m::Float64)
    sigma0 * (V0/V)^(1/m)
end

"Griffith 断裂准则: σ_cr = √(2Eγ/(πa))"
function griffith_critical_stress(E::Float64, gamma_s::Float64, a::Float64)
    sqrt(2*E*gamma_s/(pi*a))
end

"混合模式断裂 (mode I+II): 最大切向应力准则"
function mixed_mode_fracture_angle(K1::Float64, K2::Float64)
    # θ_c = 2 arctan[(1 - √(1+8(K2/K1)²))/(4 K2/K1)]
    r = K2 / max(abs(K1), 1e-15)
    theta = 2*atan((1 - sqrt(1+8*r^2))/(4*r))
    theta
end

export stress_intensity_factor_mode1, stress_intensity_factor_mode2
export energy_release_rate, fracture_criterion, paris_law
export plastic_zone_size, j_integral_estimate, ctod_estimate
export weibull_failure_probability, weibull_scale_parameter
export griffith_critical_stress, mixed_mode_fracture_angle
