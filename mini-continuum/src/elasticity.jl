# elasticity.jl — 弹性本构关系: Hooke定律, 弹性常数互转, 各向异性
# 参考: Landau Vol.7, Ting "Anisotropic Elasticity"

using LinearAlgebra

# ============================================================
# 各向同性 Hooke 定律
# ============================================================
"σ = λ tr(ε)I + 2μ ε"
function hookes_law_isotropic(strain::StrainTensor, lam::Float64, mu::Float64)
    eps=strain_matrix(strain); tr=eps[1,1]+eps[2,2]+eps[3,3]
    StressTensor(lam*tr*Matrix{Float64}(I,3,3) + 2mu*eps)
end

"ε = 1/(2μ)σ - λ/(2μ(3λ+2μ))tr(σ)I"
function hookes_law_inverse(sigma::StressTensor, lam::Float64, mu::Float64)
    s=stress_matrix(sigma); tr=s[1,1]+s[2,2]+s[3,3]
    StrainTensor( (1/(2mu))*s - (lam/(2mu*(3lam+2mu)))*tr*Matrix{Float64}(I,3,3) )
end

"Lame 常数: λ=Eν/((1+ν)(1-2ν)), μ=E/(2(1+ν))"
lame_constants(E::Float64, nu::Float64) = E*nu/((1+nu)*(1-2nu)), E/(2*(1+nu))

"工程常数 ← Lame: E=μ(3λ+2μ)/(λ+μ), ν=λ/(2(λ+μ))"
function engineering_constants(lam::Float64, mu::Float64)
    mu*(3lam+2mu)/(lam+mu), lam/(2*(lam+mu))
end

# ============================================================
# 弹性常数互转
# ============================================================
"完整弹性常数表: E,ν,G,K,λ,μ"
function elastic_constants_table(E::Float64, nu::Float64)
    lam,mu=lame_constants(E,nu); G=mu; K=E/(3*(1-2nu))
    Dict(:E=>E,:nu=>nu,:G=>G,:K=>K,:lam=>lam,:mu=>mu)
end

"从 K,G 反推 E,ν"
function from_bulk_shear(K::Float64, G::Float64)
    E=9K*G/(3K+G); nu=(3K-2G)/(2*(3K+G)); E,nu
end

"从 λ,μ 计算所有"
function from_lame(lam::Float64, mu::Float64)
    E,nu=engineering_constants(lam,mu); K=lam+2mu/3
    Dict(:E=>E,:nu=>nu,:G=>mu,:K=>K,:lam=>lam,:mu=>mu)
end

# ============================================================
# 平面问题
# ============================================================
"平面应力刚度: Q = E/(1-ν²)[1 ν 0; ν 1 0; 0 0 (1-ν)/2]"
function plane_stress_stiffness(E::Float64, nu::Float64)
    f=E/(1-nu^2); [f f*nu 0; f*nu f 0; 0 0 f*(1-nu)/2]
end

"平面应变刚度: D = E/[(1+ν)(1-2ν)][1-ν ν 0; ν 1-ν 0; 0 0 (1-2ν)/2]"
function plane_strain_stiffness(E::Float64, nu::Float64)
    f=E/((1+nu)*(1-2nu)); [f*(1-nu) f*nu 0; f*nu f*(1-nu) 0; 0 0 f*(1-2nu)/2]
end

"平面应力→平面应变 (E,ν)"
function plane_stress_to_strain(E::Float64, nu::Float64)
    E/(1-nu^2), nu/(1-nu)
end

# ============================================================
# 各向异性弹性
# ============================================================
"正交各向异性 (9常数): 柔度矩阵 [S]"
function orthotropic_compliance(E1::Float64, E2::Float64, E3::Float64,
                                nu12::Float64, nu13::Float64, nu23::Float64,
                                G12::Float64, G13::Float64, G23::Float64)
    nu21=nu12*E2/E1; nu31=nu13*E3/E1; nu32=nu23*E3/E2
    [1/E1 -nu21/E2 -nu31/E3 0 0 0;
     -nu12/E1 1/E2 -nu32/E3 0 0 0;
     -nu13/E1 -nu23/E2 1/E3 0 0 0;
     0 0 0 1/G12 0 0;
     0 0 0 0 1/G13 0;
     0 0 0 0 0 1/G23]
end

"横观各向同性 (5常数)"
function transversely_isotropic_stiffness(E1::Float64, E3::Float64, nu12::Float64, nu31::Float64, G31::Float64)
    nu13=nu31*E1/E3; Delta=(1-nu12-2*nu13*nu31)
    C11=E1*(1-nu13*nu31)/Delta; C12=E1*(nu12+nu13*nu31)/Delta
    C13=E1*nu31/Delta; C33=E3*(1-nu12)/Delta; C44=G31; C66=E1/(2*(1+nu12))
    [C11 C12 C13 0 0 0; C12 C11 C13 0 0 0; C13 C13 C33 0 0 0;
     0 0 0 C44 0 0; 0 0 0 0 C44 0; 0 0 0 0 0 C66]
end

# ============================================================
# 能量与热弹性
# ============================================================
"应变能密度: W = ½λ(trε)² + μ tr(ε²)"
function strain_energy_density(strain::StrainTensor, lam::Float64, mu::Float64)
    eps=strain_matrix(strain); tr=eps[1,1]+eps[2,2]+eps[3,3]
    0.5*lam*tr^2 + mu*sum(eps.^2)
end

"互补能密度: W* = ½σᵢⱼSᵢⱼₖₗσₖₗ"
function complementary_energy(sigma::StressTensor, E::Float64, nu::Float64)
    s=stress_matrix(sigma); tr=s[1,1]+s[2,2]+s[3,3]
    (1/(2E))*sum(s.^2) - (nu/(2E))*tr^2
end

"热弹性 Hooke (含热应变): σ = λ tr(ε)I + 2με - (3λ+2μ)αΔT I"
function thermoelastic_stress(strain::StrainTensor, lam::Float64, mu::Float64, alpha::Float64, delta_T::Float64)
    sigma=hookes_law_isotropic(strain,lam,mu)
    thermal = (3lam+2mu)*alpha*delta_T
    StressTensor(sigma.xx-thermal, sigma.yy-thermal, sigma.zz-thermal, sigma.xy, sigma.xz, sigma.yz)
end

# ============================================================
# 3D各向同性刚度/柔度矩阵 (Voigt 6x6)
# ============================================================
"刚度矩阵 C_6x6 (Voigt): sigma_i = C_ij eps_j"
function stiffness_isotropic_6x6(lam, mu)
    c11 = lam + 2mu; c12 = lam; c44 = mu
    [c11 c12 c12 0 0 0; c12 c11 c12 0 0 0; c12 c12 c11 0 0 0;
     0 0 0 c44 0 0; 0 0 0 0 c44 0; 0 0 0 0 0 c44]
end
"柔度矩阵 S_6x6 (Voigt)"
function compliance_isotropic_6x6(E, nu)
    s11 = 1/E; s12 = -nu/E; s44 = 2*(1+nu)/E
    [s11 s12 s12 0 0 0; s12 s11 s12 0 0 0; s12 s12 s11 0 0 0;
     0 0 0 s44 0 0; 0 0 0 0 s44 0; 0 0 0 0 0 s44]
end
"Bond变换: 旋转刚度矩阵 theta about z-axis"
function bond_transformation_6x6(C, theta)
    c=cos(theta); s=sin(theta)
    T = [c^2 s^2 0 0 0 2c*s; s^2 c^2 0 0 0 -2c*s; 0 0 1 0 0 0;
         0 0 0 c -s 0; 0 0 0 s c 0; -c*s c*s 0 0 0 c^2-s^2]
    T' * C * T
end
"柔度→刚度: C = inv(S)"
function compliance_to_stiffness(S)
    inv(S)
end

# ============================================================
# Saint-Venant 相容方程 (2D)
# d^2(eps_xx)/dy^2 + d^2(eps_yy)/dx^2 = 2 d^2(eps_xy)/dxdy
# 保证位移场单值性
# ============================================================
"Airy应力函数偏导数 → 应力分量"
function airy_stress(d2phi_dx2, d2phi_dy2, d2phi_dxdy)
    # sigma_xx = d^2(phi)/dy^2
    # sigma_yy = d^2(phi)/dx^2
    # sigma_xy = -d^2(phi)/dxdy
    (sig_xx=d2phi_dy2, sig_yy=d2phi_dx2, sig_xy=-d2phi_dxdy)
end
"重调和方程检查 (无体力的相容性): nabla^4 phi = 0"
function biharmonic_check(d4phi_dx4, d4phi_dx2dy2, d4phi_dy4)
    d4phi_dx4 + 2*d4phi_dx2dy2 + d4phi_dy4  # 应为 0
end
"多项式Airy函数: phi = sum A_mn x^m y^n"
function airy_polynomial(x, y, coeffs)
    result = 0.0
    for ((m,n), c) in coeffs
        result += c * x^m * y^n
    end; result
end

# ============================================================
# 各向异性弹性 (Cubic对称性, 3常数)
# ============================================================
"Cubic刚度: C11, C12, C44"
function cubic_stiffness_6x6(C11, C12, C44)
    [C11 C12 C12 0 0 0; C12 C11 C12 0 0 0; C12 C12 C11 0 0 0;
     0 0 0 C44 0 0; 0 0 0 0 C44 0; 0 0 0 0 0 C44]
end
"Cubic Zener各向异性比: A = 2*C44/(C11-C12), A=1→各向同性"
function zener_ratio(C11, C12, C44)
    2*C44/(C11-C12)
end

# ============================================================
# 非线性弹性 (Neo-Hookean, Mooney-Rivlin)
# ============================================================
"Neo-Hookean应变能 (不可压缩): W = C10*(I1-3)"
function neo_hookean_energy(C10, F)
    C = right_cauchy_green(F)
    I1 = tr(C)
    C10 * (I1 - 3)
end
"Mooney-Rivlin: W = C10*(I1-3) + C01*(I2-3)"
function mooney_rivlin_energy(C10, C01, F)
    C = right_cauchy_green(F)
    I1 = tr(C)
    I2 = 0.5*(I1^2 - tr(C*C))
    C10*(I1-3) + C01*(I2-3)
end
"Neo-Hookean第二PK应力: S = 2 dW/dC"
function neo_hookean_pk2(C10, F)
    C = right_cauchy_green(F); Cinv = inv(C)
    I1 = tr(C)
    2*C10*(Matrix{Float64}(I,3,3) - Cinv*I1/3)  # 不可压缩修正
end

export hookes_law_isotropic, hookes_law_inverse, lame_constants, engineering_constants
export elastic_constants_table, from_bulk_shear, from_lame
export plane_stress_stiffness, plane_strain_stiffness, plane_stress_to_strain
export orthotropic_compliance, transversely_isotropic_stiffness
export strain_energy_density, complementary_energy, thermoelastic_stress
export stiffness_isotropic_6x6, compliance_isotropic_6x6, bond_transformation_6x6
export compliance_to_stiffness
export airy_stress, biharmonic_check, airy_polynomial
export cubic_stiffness_6x6, zener_ratio
export neo_hookean_energy, mooney_rivlin_energy, neo_hookean_pk2
