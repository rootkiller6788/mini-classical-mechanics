# plasticity.jl — 塑性力学: 屈服面, 流动法则, 硬化
# 参考: Hill, Lubliner, Simo & Hughes

"von Mises 屈服函数: f = √(3J₂) - σ_y"
function von_mises_yield_function(J2::Float64, sigma_y::Float64)
    sqrt(3*J2) - sigma_y
end

"Tresca 屈服函数"
function tresca_yield_function(s1::Float64, s3::Float64, sigma_y::Float64)
    0.5*abs(s1-s3) - sigma_y
end

"关联流动法则 (associated flow): dεᵖ = dλ ∂f/∂σ"
function associated_flow_direction(sigma::StressTensor)
    sdev=deviatoric_stress(sigma)
    StressTensor(sdev.xx,sdev.yy,sdev.zz,sdev.xy,sdev.xz,sdev.yz)
end

"等向硬化 (isotropic hardening): σ_y = σ_y0 + H ε̄ᵖ"
function isotropic_hardening(sigma_y0::Float64, H::Float64, eq_plastic_strain::Float64)
    sigma_y0 + H * eq_plastic_strain
end

"随动硬化 (kinematic hardening): f = √(3J₂(s-α)) - σ_y, α̇ = C ε̇ᵖ"
function kinematic_hardening_backstress(C::Float64, plastic_strain::Float64)
    C * plastic_strain
end

"塑性乘子 (plastic multiplier) 一致性条件"
function plastic_multiplier(H::Float64, strain_increment::Float64, C_elastic::Float64)
    # dλ = (C_elastic * dε) / (C_elastic + H) for 1D
    C_elastic * strain_increment / (C_elastic + H)
end

"J2 塑性 (径向返回算法)"
function radial_return(sigma_trial::StressTensor, sigma_y::Float64, G::Float64)
    s_trial = deviatoric_stress(sigma_trial)
    phi_trial = sqrt(3*stress_invariants(sigma_trial).J2) - sigma_y
    if phi_trial <= 0.0; return sigma_trial, 0.0; end
    # 径向返回
    dlambda = phi_trial / (3*G)  # 塑性乘子
    factor = (sigma_y)/(sqrt(3*stress_invariants(sigma_trial).J2))
    s_new = StressTensor(s_trial.xx*factor,s_trial.yy*factor,s_trial.zz*factor,
                          s_trial.xy*factor,s_trial.xz*factor,s_trial.yz*factor)
    p = hydrostatic_pressure(sigma_trial)
    sigma_new = StressTensor(s_new.xx+p,s_new.yy+p,s_new.zz+p,s_new.xy,s_new.xz,s_new.yz)
    return sigma_new, dlambda
end

"Drucker-Prager 塑性 (地质材料)"
function drucker_prager_yield(sigma::StressTensor, alpha::Float64, k::Float64)
    inv=stress_invariants(sigma); alpha*inv.I1 + sqrt(inv.J2) - k
end

"Mohr-Coulomb 塑性 (六边形屈服面)"
function mohr_coulomb_yield(s1::Float64, s3::Float64, phi::Float64, c::Float64)
    (s1 - s3) - (s1 + s3)*sin(phi) - 2*c*cos(phi)
end

"Cam-Clay 模型 (临界状态土力学)"
function cam_clay_yield(p::Float64, q::Float64, M::Float64, p0::Float64)
    q^2/(M^2) + p*(p - p0)
end

"塑性功: Wᵖ = ∫ σ:dεᵖ"
plastic_work(sigma::StressTensor, d_eps_p::StrainTensor) = sigma.xx*d_eps_p.xx + sigma.yy*d_eps_p.yy + sigma.zz*d_eps_p.zz + 2*(sigma.xy*d_eps_p.xy + sigma.xz*d_eps_p.xz + sigma.yz*d_eps_p.yz)

"应变硬化率 (tangent modulus): E_tan = dσ/dε"
tangent_modulus(E::Float64, H::Float64) = E*H/(E+H)

export von_mises_yield_function, tresca_yield_function, associated_flow_direction
export isotropic_hardening, kinematic_hardening_backstress, plastic_multiplier
export radial_return, drucker_prager_yield, mohr_coulomb_yield, cam_clay_yield
export plastic_work, tangent_modulus
