# contact.jl — 接触力学: Hertz接触, 粘着, 摩擦
# 参考: Johnson "Contact Mechanics", Popov

"Hertz 球-球接触 (弹性): a = (3PR*/4E*)^(1/3), 接触半径"
function hertz_contact_radius(P::Float64, R1::Float64, R2::Float64, E1::Float64, E2::Float64, nu1::Float64, nu2::Float64)
    R_star = 1/(1/R1 + 1/R2)
    E_star = 1/((1-nu1^2)/E1 + (1-nu2^2)/E2)
    (3*P*R_star/(4*E_star))^(1/3)
end

"Hertz 最大接触压力: p_max = 3P/(2πa²)"
function hertz_max_pressure(P::Float64, a::Float64)
    3*P/(2*pi*a^2)
end

"Hertz 压入深度: δ = a²/R*"
function hertz_indentation(a::Float64, R_star::Float64)
    a^2/R_star
end

"JKR 粘着理论: 包含表面能的接触半径"
function jkr_contact_radius(P::Float64, R_star::Float64, E_star::Float64, delta_gamma::Float64)
    # 由 P + 3πΔγR + √(6πΔγRP + (3πΔγR)²) 决定
    P_eff = P + 3*pi*delta_gamma*R_star + sqrt(6*pi*delta_gamma*R_star*P + (3*pi*delta_gamma*R_star)^2)
    (3*P_eff*R_star/(4*E_star))^(1/3)
end

"JKR pull-off 力: P_c = -3πΔγR/2"
function jkr_pull_off_force(R_star::Float64, delta_gamma::Float64)
    -1.5*pi*delta_gamma*R_star
end

"DMT 理论 pull-off 力: P_c = -2πΔγR"
function dmt_pull_off_force(R_star::Float64, delta_gamma::Float64)
    -2*pi*delta_gamma*R_star
end

"Coulomb 摩擦: F_f = μ |F_n| sign(v)"
coulomb_friction(mu::Float64, Fn::Float64, v::Float64) = mu*abs(Fn)*sign(v)

"Stribeck 曲线 (摩擦 vs 速度): μ = μ_c + (μ_s-μ_c)exp(-|v|/v_s)"
function stribeck_friction(mu_s::Float64, mu_c::Float64, v_s::Float64, v::Float64)
    mu_c + (mu_s-mu_c)*exp(-abs(v)/v_s)
end

"Mindlin 切向接触 (微滑移)"
function mindlin_tangential_stiffness(G_star::Float64, a::Float64, nu::Float64)
    8*G_star*a/(2-nu)
end

"等效模量 (接触问题)"
function contact_effective_modulus(E1::Float64, E2::Float64, nu1::Float64, nu2::Float64)
    1.0/((1-nu1^2)/E1 + (1-nu2^2)/E2)
end

"接触压力的椭圆形分布: p(r) = p_max √(1-r²/a²)"
function hertz_pressure_distribution(p_max::Float64, r::Float64, a::Float64)
    abs(r) <= a ? p_max*sqrt(1-(r/a)^2) : 0.0
end

"线接触 (圆柱): a = √(4PR*/(πE*L))"
function line_contact_width(P::Float64, R_star::Float64, E_star::Float64, L::Float64)
    sqrt(4*P*R_star/(pi*E_star*L))
end

"表面能 (Dupre粘着功): W_ad = gamma_1 + gamma_2 - gamma_12"
function dupre_adhesion_work(gamma1, gamma2, gamma12)
    gamma1 + gamma2 - gamma12
end
"Tabor参数: mu = (R*Delta_gamma^2/(E*^2*z0^3))^(1/3)"
function tabor_parameter(R_star, Delta_gamma, E_star, z0)
    (R_star * Delta_gamma^2 / (E_star^2 * z0^3))^(1/3)
end
"Maugis-Dugdale过渡: JKR(mu>>1) → DMT(mu<<1)"
function maugis_dugdale_transition(mu)
    if mu > 5.0; :JKR; elseif mu < 0.1; :DMT; else; :transition; end
end
"Cattaneo-Mindlin微滑移: 切向力-位移关系"
function cattaneo_mindlin_traction(mu, Fn, delta_t, a, G_star)
    F_max = mu * Fn
    F_t = F_max * (1 - (1 - 16*G_star*a*delta_t/(3*mu*Fn))^(3/2))
    min(F_t, F_max)
end
"滚动接触 (Carter): 蠕滑率-力关系"
function carter_creep_force(xi, a, b, G)
    0.0  # 需要数值积分
end
"磨损 (Archard): V = K*F_n*s/H"
function archard_wear(K, Fn, s, H)
    K * Fn * s / H
end
"粗糙表面接触 (Greenwood-Williamson): A_r = pi*eta*A_nom*R*integral(z-d)^*phi(z)dz"
function greenwood_williamson_area(eta, A_nom, R, sigma, d)
    # 简化: 高斯高度分布
    A_nom * eta * pi * R * sigma * exp(-d^2/(2*sigma^2))
end

export hertz_contact_radius, hertz_max_pressure, hertz_indentation
export jkr_contact_radius, jkr_pull_off_force, dmt_pull_off_force
export coulomb_friction, stribeck_friction, mindlin_tangential_stiffness
export contact_effective_modulus, hertz_pressure_distribution, line_contact_width
export dupre_adhesion_work, tabor_parameter, maugis_dugdale_transition
export cattaneo_mindlin_traction, carter_creep_force, archard_wear
export greenwood_williamson_area
