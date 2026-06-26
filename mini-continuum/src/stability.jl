# stability.jl — 结构稳定性: 屈曲模态, 后屈曲, 动力稳定性, 分岔
# 参考: Timoshenko & Gere, Bazant & Cedolin, Thompson & Hunt
using LinearAlgebra

"简支柱屈曲载荷: P_cr(n) = n^2 pi^2 EI/L^2"
function column_buckling_load(EI, L, n)
    (n*pi/L)^2 * EI
end
"固支-自由柱: P_cr = pi^2 EI/(4L^2)"
function cantilever_buckling_load(EI, L)
    pi^2 * EI / (4*L^2)
end
"固支-固支柱: P_cr = 4 pi^2 EI/L^2"
function fixed_fixed_buckling_load(EI, L)
    4*pi^2 * EI / L^2
end
"固支-铰支柱: P_cr ~ 20.19 EI/L^2"
function fixed_pinned_buckling_load(EI, L)
    20.19 * EI / L^2
end
"Southwell图: delta/P = (1/P_cr)*delta + delta_0/P_cr (初始缺陷法)"
function southwell_plot(delta, P; delta_0=0.0)
    delta ./ P
end
"Perry-Robertson公式 (初始弯曲柱): sigma_max = sigma_cr/2*(1+eta+sqrt((1+eta)^2-4*sigma_cr/sigma_y))"
function perry_robertson(sigma_cr, sigma_y, eta)
    0.5*sigma_cr*(1+eta+sqrt((1+eta)^2 - 4*sigma_cr/sigma_y))
end
"初始缺陷参数: eta = delta_0*c/r^2"
function imperfection_parameter(delta_0, c, r)
    delta_0 * c / r^2
end
"板屈曲载荷 (简支矩形): N_cr = k pi^2 D/b^2, k=(m*b/a+a/(m*b))^2"
function plate_buckling_load(D, a, b, m)
    k = (m*b/a + a/(m*b))^2
    k * pi^2 * D / b^2
end
"板屈曲系数k (各种边界): ssss, cccc, scsc等"
function plate_buckling_coefficient(boundary, a_over_b)
    coeffs = Dict(:ssss=>4.0, :cccc=>6.97, :scsc=>5.34, :ccss=>5.0)
    k0 = get(coeffs, boundary, 4.0)
    k0  # 简化; 实际k是a/b的函数
end
"壳屈曲 (轴压圆柱): sigma_cr = E*h/(R*sqrt(3*(1-nu^2))) (经典)"
function shell_buckling_classical(E, nu, h, R)
    E*h/(R*sqrt(3*(1-nu^2)))
end
"壳屈曲 (考虑缺陷): sigma_cr = alpha * sigma_classical, alpha=0.2~0.5"
function shell_buckling_knockdown(alpha, sigma_classical)
    alpha * sigma_classical
end
"Koiter后屈曲理论: P/P_cr = 1 + a*delta + b*delta^2 + ..."
function koiter_postbuckling(delta, P_cr, a, b)
    P_cr * (1 + a*delta + b*delta^2)
end
"后屈曲稳定性: b>0稳定(板), b<0不稳定(壳)"
function postbuckling_stability(b)
    b > 0 ? :stable : :unstable
end
"分岔类型: 超临界分岔 (stable), 亚临界分岔 (unstable)"
function bifurcation_type(a, b)
    if a == 0.0 && b > 0.0; :supercritical_pitchfork
    elseif a == 0.0 && b < 0.0; :subcritical_pitchfork
    elseif a != 0.0; :transcritical
    else; :unknown; end
end
"动力屈曲 (Budiansky-Roth): 阶跃载荷 P=P0*H(t)"
function dynamic_buckling_load(P_static_cr, omega_n, t_duration)
    # 简化: P_dyn_cr ~ P_static_cr * f(omega_n*t)
    tau = omega_n * t_duration
    factor = tau < pi ? 2.0 : 1.0  # 脉冲载荷近似
    P_static_cr * factor
end
"稳定性行列式法: det(K - lambda*Kg) = 0"
function stability_determinant(K, Kg, lambda)
    det(K - lambda * Kg)
end
"几何刚度矩阵Kg (杆单元, 2x2): Kg = P/L * [1 -1; -1 1]"
function geometric_stiffness_bar(P, L)
    P/L * [1 -1; -1 1]
end
"切线刚度: K_tan = K_elastic + K_geometric"
function tangent_stiffness(K_elastic, K_geometric)
    K_elastic + K_geometric
end
"切线刚度正定性检查: min(eigvals(K_tan)) > 0 → 稳定"
function is_stable(K_tan)
    minimum(eigvals(K_tan)) > 0
end
"弧长法增量 (Riks): sqrt(delta_u^2 + delta_lambda^2*psi) = ds"
function arc_length_constraint(delta_u, delta_lambda, psi, ds)
    sqrt(delta_u'*delta_u + delta_lambda^2*psi) - ds
end

export column_buckling_load, cantilever_buckling_load
export fixed_fixed_buckling_load, fixed_pinned_buckling_load
export southwell_plot, perry_robertson, imperfection_parameter
export plate_buckling_load, plate_buckling_coefficient
export shell_buckling_classical, shell_buckling_knockdown
export koiter_postbuckling, postbuckling_stability, bifurcation_type
export dynamic_buckling_load, stability_determinant
export geometric_stiffness_bar, tangent_stiffness, is_stable, arc_length_constraint
