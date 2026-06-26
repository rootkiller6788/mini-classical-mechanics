# composites.jl — 复合材料力学: 混合律, 层合板, 微观力学
# 参考: Jones "Mechanics of Composite Materials", Halpin-Tsai

"混合律 (Rule of Mixtures): E₁ = V_f E_f + V_m E_m (纵向)"
function rule_of_mixtures_E1(Vf::Float64, Ef::Float64, Em::Float64)
    Vf*Ef + (1-Vf)*Em
end

"横向模量 (Reuss): 1/E₂ = V_f/E_f + V_m/E_m"
function rule_of_mixtures_E2(Vf::Float64, Ef::Float64, Em::Float64)
    1/(Vf/Ef + (1-Vf)/Em)
end

"Halpin-Tsai 方程 (半经验): M/M_m = (1+ξηV_f)/(1-ηV_f)"
function halpin_tsai(Mf::Float64, Mm::Float64, Vf::Float64, xi::Float64)
    eta = (Mf/Mm - 1)/(Mf/Mm + xi)
    Mm * (1 + xi*eta*Vf) / (1 - eta*Vf)
end

"层合板 ABD 矩阵: N = Aε⁰ + Bκ, M = Bε⁰ + Dκ"
function laminate_abd_matrix(plies::Vector{Tuple{Float64,Float64,Matrix{Float64}}})
    # (theta, thickness, Q) for each ply
    n=length(plies); A=zeros(3,3); B=zeros(3,3); D=zeros(3,3)
    z0=0.0; z_prev=0.0
    for (theta,t,Q) in plies
        z1=z_prev - t; z_mid=(z0+z1)/2
        # Transform Q to laminate coordinates
        c=cos(theta); s=sin(theta)
        T=[c^2 s^2 2c*s; s^2 c^2 -2c*s; -c*s c*s c^2-s^2]
        Qbar=T'*Q*T
        A+=Qbar*(z0-z1); B+=Qbar*(z0^2-z1^2)/2; D+=Qbar*(z0^3-z1^3)/3
        z_prev=z1; z0=z1
    end; A, B, D
end

"经典层合板理论: ε(z) = ε⁰ + zκ"
function laminate_strain(epsilon0::Vector{Float64}, kappa::Vector{Float64}, z::Float64)
    epsilon0 + z*kappa
end

"单层板强度 (Tsai-Hill 准则): σ₁²/X² - σ₁σ₂/X² + σ₂²/Y² + τ₁₂²/S² = 1"
function tsai_hill_criterion(sigma1::Float64, sigma2::Float64, tau12::Float64, X::Float64, Y::Float64, S::Float64)
    sigma1^2/X^2 - sigma1*sigma2/X^2 + sigma2^2/Y^2 + tau12^2/S^2
end

"Tsai-Wu 准则: F_i σ_i + F_ij σ_i σ_j = 1"
function tsai_wu_criterion(sigma1::Float64, sigma2::Float64, tau12::Float64, Xt::Float64, Xc::Float64, Yt::Float64, Yc::Float64, S::Float64)
    F1=1/Xt-1/Xc; F2=1/Yt-1/Yc; F11=1/(Xt*Xc); F22=1/(Yt*Yc); F66=1/S^2; F12=-0.5*sqrt(F11*F22)
    F1*sigma1+F2*sigma2+F11*sigma1^2+F22*sigma2^2+F66*tau12^2+2*F12*sigma1*sigma2
end

"Hashin 失效准则 (纤维压缩/拉伸, 基体压缩/拉伸)"
function hashin_fiber_tension(sigma1::Float64, tau12::Float64, Xt::Float64, S12::Float64)
    sigma1 > 0 ? (sigma1/Xt)^2 + (tau12/S12)^2 : 0.0
end

"混合比: V_f = ρ_c w_f / ρ_f"
volume_fraction_from_weight(wf::Float64, rhoc::Float64, rhof::Float64) = rhoc*wf/rhof

"比模量: E/ρ"
specific_modulus(E::Float64, rho::Float64) = E/rho

export rule_of_mixtures_E1, rule_of_mixtures_E2, halpin_tsai
export laminate_abd_matrix, laminate_strain
export tsai_hill_criterion, tsai_wu_criterion, hashin_fiber_tension
export volume_fraction_from_weight, specific_modulus
"混合律扩展 — 短纤维复合材料 (Cox 剪切滞后模型)"
function cox_shear_lag(Ef::Float64, Em::Float64, Vf::Float64, L::Float64, d::Float64, Gm::Float64)
    s = L/d  # 长径比
    beta = sqrt(2*Gm/(Ef*log(pi/(4*Vf))))
    eta_l = 1 - tanh(beta*s/2)/(beta*s/2)
    Em*Vf*eta_l*Ef*Vf + Em*(1-Vf)  # 简化
end

"Eshelby 夹杂理论 (椭圆体夹杂)"
function eshelby_tensor(aspect_ratio::Float64, nu::Float64)
    # 对于球体: S1111 = (7-5ν)/(15(1-ν))
    S1111 = (7-5*nu)/(15*(1-nu))
    S1122 = (5*nu-1)/(15*(1-nu))
    S1212 = (4-5*nu)/(15*(1-nu))
    S1111, S1122, S1212
end

"Mori-Tanaka 方法有效模量"
function mori_tanaka_modulus(Km::Float64, Gm::Float64, Ki::Float64, Gi::Float64, Vf::Float64)
    K_eff = Km + Vf*(Ki-Km)*Km/(Km + (1-Vf)*(Ki-Km)*(1+Km/(Km+4*Gm/3))/3)
    K_eff
end

"纤维临界长度: L_c = σ_f d/(2τ)"
critical_fiber_length(sigma_f::Float64, d::Float64, tau::Float64) = sigma_f*d/(2*tau)

"层间应力 (自由边界效应)"
function free_edge_stress(Qbar::Matrix{Float64}, delta_T::Float64, alpha::Vector{Float64})
    sigma = -Qbar * (alpha * delta_T)
    sigma
end

export cox_shear_lag, eshelby_tensor, mori_tanaka_modulus
export critical_fiber_length, free_edge_stress
