# conservation.jl — 守恒定律与平衡方程
# 参考: Holzapfel, Gurtin, Truesdell & Noll
using LinearAlgebra

"连续性方程残差: drho/dt + rho*div(v)"
function continuity_residual(rho, drho_dt, div_v)
    drho_dt + rho * div_v
end
"物质密度: rho_0 = J * rho"
material_density(rho, J) = J * rho
"质量守恒误差 (1D积分)"
function mass_conservation_error(rho_field, dx)
    sum(rho_field) * dx
end
"Cauchy平衡残差 1D: d(sigma)/dx + rho*b"
function equilibrium_residual_1d(sigma, dx, rho, b)
    n = length(sigma); r = zeros(n)
    for i in 2:n-1; r[i] = (sigma[i+1]-sigma[i-1])/(2*dx) + rho*b; end
    return r
end
"3D平衡检查: div(sigma) + rho*b = 0"
function equilibrium_3d(sigma, body_force, grad_sigma_x, grad_sigma_y, grad_sigma_z, rho)
    bx, by, bz = body_force
    rx = grad_sigma_x + rho*bx
    ry = grad_sigma_y + rho*by
    rz = grad_sigma_z + rho*bz
    norm([rx, ry, rz])
end
"动力学动量残差: rho*a - div(sigma) - rho*b"
function momentum_residual(rho, a, div_sigma, b)
    [rho*a[i] - div_sigma[i] - rho*b[i] for i in 1:3]
end
"角动量守恒 → Cauchy应力对称性检查"
function angular_momentum_symmetry_check(sigma)
    sigma.xy, sigma.xz, sigma.yz
end
"内能变化率: rho*de/dt = sigma:D - div(q) + rho*r"
function energy_balance_rate(sigma, D, rho, r, div_q)
    sd = sigma.xx*D[1,1]+sigma.yy*D[2,2]+sigma.zz*D[3,3] + 2*(sigma.xy*D[1,2]+sigma.xz*D[1,3]+sigma.yz*D[2,3])
    sd + rho*r - div_q
end
"应变能变化率: dW/dt = sigma:D (应力功)"
function strain_energy_rate(sigma, D)
    sigma.xx*D[1,1]+sigma.yy*D[2,2]+sigma.zz*D[3,3] + 2*(sigma.xy*D[1,2]+sigma.xz*D[1,3]+sigma.yz*D[2,3])
end
"非弹性耗散: D_mech = sigma:D - dPsi/dt >= 0"
function mechanical_dissipation(sigma, D, dPsi_dt)
    strain_energy_rate(sigma, D) - dPsi_dt
end
"Clausius-Duhem不等式: sigma:D - rho*(dPsi/dt+s*dT/dt) - q·grad(T)/T >= 0"
function clausius_duhem(sigma, D, rho, dPsi_dt, s, dT_dt, q_vec, grad_T, T)
    mech = strain_energy_rate(sigma, D) - rho*dPsi_dt
    thermal = -rho*s*dT_dt - dot(q_vec, grad_T)/T
    mech + thermal
end
"Fourier不等式: q·grad(T) <= 0"
function fourier_check(q_vec, grad_T)
    dot(q_vec, grad_T)
end
"物质时间导数: Df/Dt = partial f/partial t + v·grad(f)"
function material_time_derivative(df_dt, v, grad_f)
    df_dt + dot(v, grad_f)
end
"物质加速度: a = dv/dt + grad(v)·v"
function material_acceleration(dv_dt, grad_v, v)
    dv_dt + grad_v * v
end
"第一Piola-Kirchhoff应力: P = J sigma F^{-T} (两点张量)"
function first_pk_stress(sigma, F)
    J = defgrad_det(F); Fm = defgrad_to_matrix(F)
    sm = stress_matrix(sigma); J * sm * inv(Fm)'
end
"第二Piola-Kirchhoff应力: S = J F^{-1} sigma F^{-T}"
function second_pk_stress(sigma, F)
    J = defgrad_det(F); Fm = defgrad_to_matrix(F)
    sm = stress_matrix(sigma); J * inv(Fm) * sm * inv(Fm)'
end
"Kirchhoff应力: tau = J sigma"
function kirchhoff_stress(sigma, J)
    J * stress_matrix(sigma)
end
"功共轭验证: P:dF = sigma:D"
function power_conjugacy(sigma, D, P, F_dot)
    pow1 = strain_energy_rate(sigma, D)
    pow2 = sum(P .* F_dot)
    abs(pow1 - pow2)
end
"Eshelby能量-动量张量: Sigma = Psi*I - F^T P"
function eshelby_stress(Psi, F, P)
    Fm = defgrad_to_matrix(F); Psi*Matrix{Float64}(I,3,3) - Fm'*P
end
"J积分路径贡献: J_contour = int (W*dy - T·du/dx ds)"
function j_integral_path(W, T_vec, du_dx, dy, ds)
    W*dy - dot(T_vec, du_dx)*ds
end

export continuity_residual, material_density, mass_conservation_error
export equilibrium_residual_1d, equilibrium_3d, momentum_residual
export angular_momentum_symmetry_check
export energy_balance_rate, strain_energy_rate, mechanical_dissipation
export clausius_duhem, fourier_check
export material_time_derivative, material_acceleration
export first_pk_stress, second_pk_stress, kirchhoff_stress
export power_conjugacy, eshelby_stress, j_integral_path
