# plates.jl — 板壳理论: Kirchhoff板, Mindlin板, 膜, 简正模
# 参考: Timoshenko "Theory of Plates and Shells", Leissa "Vibration of Plates"

using LinearAlgebra

# ============================================================
# 板弯曲基础
# ============================================================
"Kirchhoff 板弯曲刚度: D = E h³/(12(1-ν²))"
plate_flexural_rigidity(E::Float64, nu::Float64, h::Float64) = E*h^3/(12*(1-nu^2))

"Mindlin 板剪切刚度: S = k G h (k≈5/6 for rectangle)"
mindlin_shear_stiffness(G::Float64, h::Float64, k=5/6) = k*G*h

"板弯矩-曲率关系 (Kirchhoff): Mx = -D(∂²w/∂x² + ν∂²w/∂y²)"
function plate_moment(D::Float64, nu::Float64, d2w_dx2::Float64, d2w_dy2::Float64)
    Mx = -D*(d2w_dx2 + nu*d2w_dy2)
    My = -D*(d2w_dy2 + nu*d2w_dx2)
    Mxy = -D*(1-nu)*0.0  # 需要 ∂²w/∂x∂y
    return (Mx=Mx, My=My, Mxy=Mxy)
end

# ============================================================
# 圆板
# ============================================================
"圆板轴对称挠度 (简支, 均布载荷): w(r)= q(a²-r²)[(5+ν)a²/(1+ν)-r²]/(64D)"
function circular_plate_deflection(D::Float64, a::Float64, q::Float64, r_vals, nu::Float64)
    [q*(a^2-r^2)*((5+nu)*a^2/(1+nu)-r^2)/(64D) for r in r_vals]
end

"圆板最大挠度 (简支, 中心)"
circular_plate_max_deflection(D::Float64, a::Float64, q::Float64, nu::Float64) = q*a^4*(5+nu)/(64*D*(1+nu))

"圆板最大挠度 (固支): w_max = qa⁴/(64D)"
clamped_circular_plate_max(D::Float64, a::Float64, q::Float64) = q*a^4/(64D)

"圆板中心集中力: w(r) = P/(16πD)[(3+ν)a²/(1+ν) - r² + 2r²ln(r/a)]"
function circular_plate_point_load(D::Float64, a::Float64, P::Float64, r_vals, nu::Float64)
    w_center = P*a^2*(3+nu)/(16pi*D*(1+nu))
    [w_center - P*r^2/(16pi*D)*(1 + 2*log(r/a)) for r in r_vals]
end

"圆板临界屈曲载荷 (简支): P_cr = 14.68 D/a²"
circular_plate_buckling(D::Float64, a::Float64) = 14.68*D/a^2

# ============================================================
# 矩形板振动
# ============================================================
"矩形板简正模频率 (简支, 边长 a×b)"
function plate_natural_frequencies(D::Float64, rho::Float64, h::Float64, a::Float64, b::Float64, m_max=3, n_max=3)
    freqs = Float64[]
    for m in 1:m_max, n in 1:n_max
        omega = ((m*pi/a)^2 + (n*pi/b)^2) * sqrt(D/(rho*h))
        push!(freqs, omega)
    end
    sort(freqs)
end

"矩形板模态形状: w_mn(x,y) = sin(mπx/a)sin(nπy/b)"
function plate_mode_shape(m::Int, n::Int, x::Float64, y::Float64, a::Float64, b::Float64)
    sin(m*pi*x/a) * sin(n*pi*y/b)
end

"矩形板简正模频率 — 其他边界条件系数"
function plate_frequency_coefficient(D::Float64, rho::Float64, h::Float64, a::Float64, b::Float64, lambda_sq::Float64)
    lambda_sq/a^2 * sqrt(D/(rho*h))
end

"常见边界 λ² 系数 (简支: π², 固支: 36.0, 自由: 13.5)"
const BOUNDARY_LAMBDA = Dict(
    :ssss => pi^2,    :cccc => 36.0,   :scsc => 28.95,
    :ccss => 23.65,   :cfcf => 22.37,  :ssff => 9.63
)

# ============================================================
# 膜振动
# ============================================================
"薄膜振动频率: ω_mn = πc√((m/a)²+(n/b)²), c = √(T/ρ)"
function membrane_frequencies(c::Float64, a::Float64, b::Float64, m_max=3, n_max=3)
    freqs = Float64[]
    for m in 1:m_max, n in 1:n_max
        push!(freqs, pi*c*sqrt((m/a)^2 + (n/b)^2))
    end
    sort(freqs)
end

"膜波速: c = √(T/ρ_s), T为张力, ρ_s为面密度"
membrane_wave_speed(T::Float64, rho_s::Float64) = sqrt(T/rho_s)

"矩形膜模态形状"
function membrane_mode_shape(m::Int, n::Int, x::Float64, y::Float64, a::Float64, b::Float64)
    sin(m*pi*x/a) * sin(n*pi*y/b)
end

"圆膜频率: ω_mn = α_mn c/R, α_mn为Bessel零点"
function circular_membrane_frequencies(c::Float64, R::Float64, m_max=3, n_max=3)
    # α_mn: Bessel J_m 的第 n 个零点 (近似)
    alpha = [2.4048 5.5201 8.6537; 3.8317 7.0156 10.1735; 5.1356 8.4172 11.6198]
    freqs = Float64[]
    for m in 1:min(m_max,3), n in 1:min(n_max,3)
        push!(freqs, alpha[m,n]*c/R)
    end
    sort(freqs)
end

# ============================================================
# 壳理论基础
# ============================================================
"圆柱壳弯曲刚度: D_shell = E h³/(12(1-ν²)) 同板"
shell_flexural_rigidity = plate_flexural_rigidity

"圆柱壳膜刚度: K_membrane = E h/(1-ν²)"
shell_membrane_stiffness(E::Float64, nu::Float64, h::Float64) = E*h/(1-nu^2)

"圆柱壳简正模频率 (Donnell-Mushtari 理论)"
function cylindrical_shell_frequency(E::Float64, nu::Float64, rho::Float64, R::Float64, h::Float64, L::Float64, m::Int, n::Int)
    D = E*h^3/(12*(1-nu^2))
    K = E*h/(1-nu^2)
    km = m*pi/L; kn = n/R
    omega_ring = sqrt(E/(rho*R^2))  # 环频率
    omega_bend = (km^2 + kn^2)^2*D/(rho*h) + E*km^4/(rho*R^2*(km^2+kn^2)^2)
    sqrt(max(omega_bend, 0.0))
end

export plate_flexural_rigidity, mindlin_shear_stiffness, plate_moment
export circular_plate_deflection, circular_plate_max_deflection, clamped_circular_plate_max
export circular_plate_point_load, circular_plate_buckling
export plate_natural_frequencies, plate_mode_shape, plate_frequency_coefficient, BOUNDARY_LAMBDA
export membrane_frequencies, membrane_wave_speed, membrane_mode_shape, circular_membrane_frequencies
export shell_flexural_rigidity, shell_membrane_stiffness, cylindrical_shell_frequency
