# micromechanics.jl — 微观力学: Eshelby夹杂, 均质化, 位错, 晶体塑性
# 参考: Mura "Micromechanics of Defects", Nemat-Nasser & Hori
using LinearAlgebra

"Eshelby张量 (球体, 各向同性基体): S_1111 = (7-5nu)/(15(1-nu))"
function eshelby_sphere_S(nu)
    S1111 = (7-5*nu)/(15*(1-nu))
    S1122 = (5*nu-1)/(15*(1-nu))
    S1212 = (4-5*nu)/(15*(1-nu))
    (S1111=S1111, S1122=S1122, S1212=S1212)
end
"Eshelby张量 (长纤维, a1>>a2=a3): S_2222=S_3333=(5-4nu)/(8(1-nu))"
function eshelby_fiber_S(nu)
    S2222 = (5-4*nu)/(8*(1-nu))
    S2233 = (4*nu-1)/(8*(1-nu))
    S2323 = (3-4*nu)/(8*(1-nu))
    S2211 = nu/(2*(1-nu))
    (S2222=S2222, S2233=S2233, S2323=S2323, S2211=S2211)
end
"Mori-Tanaka 有效体积模量: K = K0 + c*(K1-K0)*K0/(K0+alpha*(1-c)*(K1-K0))"
function mori_tanaka_bulk(K0, G0, K1, G1, c)
    alpha0 = (1+4*G0/(3*K0))/3  # 球体
    K0 + c*(K1-K0)*K0/(K0 + alpha0*(1-c)*(K1-K0))
end
"Mori-Tanaka 有效剪切模量"
function mori_tanaka_shear(K0, G0, G1, c)
    beta0 = 2*(4-5*G0/(K0+4*G0/3))/15  # 球体
    G0 + c*(G1-G0)*G0/(G0 + beta0*(1-c)*(G1-G0))
end
"自洽法 (Self-Consistent): K_sc = K0 + c*(K1-K0)*K_sc/(K_sc+...)"
function self_consistent_bulk(K0, G0, K1, G1, c; tol=1e-6, maxiter=100)
    K = K0; G = G0
    for _ in 1:maxiter
        alpha = (1+4*G/(3*K))/3
        K_new = K0 + c*(K1-K0)*K/(K+alpha*(1-c)*(K1-K0))
        if abs(K_new-K) < tol*abs(K); K = K_new; break; end
        K = K_new
    end; K
end
"Voigt上界: K_V = (1-c)*K0 + c*K1"
function voigt_bound(K0, K1, c)
    (1-c)*K0 + c*K1
end
"Reuss下界: 1/K_R = (1-c)/K0 + c/K1"
function reuss_bound(K0, K1, c)
    1/((1-c)/K0 + c/K1)
end
"Hashin-Shtrikman上下界"
function hashin_shtrikman_bounds(K0, G0, K1, G1, c)
    if K1 > K0 && G1 > G0  # 刚性夹杂
        K_lower = K0 + c/(1/(K1-K0) + 3*(1-c)/(3*K0+4*G0))
        K_upper = K1 + (1-c)/(1/(K0-K1) + 3*c/(3*K1+4*G1))
    else
        K_lower = K1 + (1-c)/(1/(K0-K1) + 3*c/(3*K1+4*G1))
        K_upper = K0 + c/(1/(K1-K0) + 3*(1-c)/(3*K0+4*G0))
    end
    (K_lower=K_lower, K_upper=K_upper)
end
"浓度张量 A: eps_inclusion = A : eps_applied (Mori-Tanaka)"
function strain_concentration_tensor(K0, G0, K1, G1)
    alpha = 1/3*(1+4*G0/(3*K0))
    beta = 2*(4-5*G0/(K0+4*G0/3))/15
    (alpha=alpha, beta=beta)
end
"特征应变 (eigenstrain): eps* = (C1-C0):eps + C0:(eps*-S:eps*)"
function eigenstrain_eshelby(C0, C1, S, eps_applied)
    deltaC = C1 - C0
    # eps* = -((C1-C0):S + C0)^{-1}:(C1-C0):eps_applied
    -((deltaC*S + C0) \ (deltaC * eps_applied))
end
"位错应力场 (刃型): sigma_xx = -D*y(3x^2+y^2)/(x^2+y^2)^2"
function dislocation_edge_stress(D, G, b, nu, x, y)
    D = G*b/(2*pi*(1-nu))
    sxx = -D*y*(3*x^2+y^2)/(x^2+y^2)^2
    syy = D*y*(x^2-y^2)/(x^2+y^2)^2
    sxy = D*x*(x^2-y^2)/(x^2+y^2)^2
    (sxx=sxx, syy=syy, sxy=sxy)
end
"位错能量 (单位长度): E = G*b^2/(4*pi*K)*ln(R/r0)"
function dislocation_energy(G, b, nu, R, r0)
    K = 1 - nu  # 刃型
    G * b^2 / (4*pi*K) * log(R/r0)
end
"Peierls应力: tau_P = 2*G/(1-nu)*exp(-2*pi*a/(b*(1-nu)))"
function peierls_stress(G, nu, a, b)
    2*G/(1-nu) * exp(-2*pi*a/(b*(1-nu)))
end
"Schmid定律: tau_RSS = sigma * cos(phi) * cos(lambda)"
function schmid_resolved_shear(sigma, phi, lambda)
    sigma * cos(phi) * cos(lambda)
end
"Schmid因子: m = cos(phi)*cos(lambda)"
function schmid_factor(phi, lambda)
    cos(phi) * cos(lambda)
end
"Taylor硬化: tau = alpha*G*b*sqrt(rho_disl)"
function taylor_hardening(alpha, G, b, rho_disl)
    alpha * G * b * sqrt(rho_disl)
end
"Orowan应力: tau_Or = G*b/L"
function orowan_stress(G, b, L)
    G * b / L
end
"Hall-Petch: sigma_y = sigma_i + k_y/sqrt(d)"
function hall_petch(sigma_i, k_y, d)
    sigma_i + k_y / sqrt(d)
end
"晶体取向 (Bunge Euler角): g = R_z(phi1)*R_x(Phi)*R_z(phi2)"
function crystal_orientation_matrix(phi1, Phi, phi2)
    Rz1 = [cos(phi1) sin(phi1) 0; -sin(phi1) cos(phi1) 0; 0 0 1]
    Rx = [1 0 0; 0 cos(Phi) sin(Phi); 0 -sin(Phi) cos(Phi)]
    Rz2 = [cos(phi2) sin(phi2) 0; -sin(phi2) cos(phi2) 0; 0 0 1]
    Rz2 * Rx * Rz1
end
"代表性体积单元 (RVE) 尺寸检查: L_RVE >= L_char*N_inclusions^(1/3)"
function rve_size_check(L_RVE, L_char, N_inclusions)
    L_min = L_char * N_inclusions^(1/3)
    L_RVE >= L_min
end

export eshelby_sphere_S, eshelby_fiber_S
export mori_tanaka_bulk, mori_tanaka_shear, self_consistent_bulk
export voigt_bound, reuss_bound, hashin_shtrikman_bounds
export strain_concentration_tensor, eigenstrain_eshelby
export dislocation_edge_stress, dislocation_energy, peierls_stress
export schmid_resolved_shear, schmid_factor, taylor_hardening
export orowan_stress, hall_petch, crystal_orientation_matrix, rve_size_check
