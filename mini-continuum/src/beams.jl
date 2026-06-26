# beams.jl — 梁理论: Euler-Bernoulli, Timoshenko, 屈曲, 稳定性
# 参考: Timoshenko, Gere & Goodno "Mechanics of Materials"

using LinearAlgebra

# ============================================================
# Euler-Bernoulli 梁
# ============================================================
"Euler-Bernoulli 梁 FDM 求解: d²/dx²(EI d²w/dx²)=q(x) (简支)"
function euler_bernoulli_deflection(EI::Float64, L::Float64, q_load::Function; nx=200)
    xs=range(0,L,length=nx); dx=L/(nx-1); A=zeros(nx,nx); b=zeros(nx)
    for i in 3:nx-2; A[i,i-2:i+2]=[1,-4,6,-4,1]/dx^4; b[i]=q_load(xs[i])/EI; end
    A[1,1]=1; b[1]=0; A[2,2]=1; b[2]=0; A[nx,nx]=1; b[nx]=0; A[nx-1,nx-1]=1; b[nx-1]=0
    w = A\b; return xs, w
end

"简支梁 均布载荷: w_max = 5qL⁴/(384EI)"
beam_max_deflection_uniform(EI::Float64, L::Float64, q::Float64) = 5*q*L^4/(384*EI)

"简支梁 集中力 中点: w_max = PL³/(48EI)"
beam_max_deflection_point(EI::Float64, L::Float64, P::Float64) = P*L^3/(48*EI)

"悬臂梁 挠度: w(x)=Px²(3L-x)/(6EI), w_max=PL³/(3EI)"
function cantilever_deflection(EI::Float64, L::Float64, P::Float64, x::Float64)
    P*x^2*(3*L-x)/(6*EI)
end

"悬臂梁 均布载荷: w(x)=qx²(6L²-4Lx+x²)/(24EI)"
function cantilever_uniform(EI::Float64, L::Float64, q::Float64, x::Float64)
    q*x^2*(6*L^2 - 4*L*x + x^2)/(24*EI)
end

"梁弯矩: M = -EI d²w/dx²"
beam_moment(EI::Float64, curvature::Float64) = -EI*curvature

"梁剪应力 (矩形截面): τ_max = 3V/(2A)"
beam_shear_stress(V::Float64, A::Float64) = 3*V/(2*A)

# ============================================================
# Timoshenko 梁 (计入剪切变形)
# ============================================================
"Timoshenko 梁挠度 (简支 中点集中力): w = PL³/(48EI) + PL/(4kGA)"
function timoshenko_deflection_point(EI::Float64, L::Float64, P::Float64, GAKs::Float64)
    w_bend = P*L^3/(48*EI)
    w_shear = P*L/(4*GAKs)
    return w_bend + w_shear, w_bend, w_shear
end

"剪切变形重要性: α = (EI/GAKs) · π²/L²"
function shear_deformation_ratio(EI::Float64, GAKs::Float64, L::Float64)
    EI/GAKs * pi^2/L^2
end

# ============================================================
# 屈曲 (Buckling)
# ============================================================
"Euler 临界屈曲载荷: P_cr = n²π²EI/(KL)²"
function euler_buckling_load(EI::Float64, L::Float64, n::Int=1, K::Float64=1.0)
    (n*pi)^2 * EI / (K*L)^2
end

"有效长度因子 K 的标准值"
const BUCKLING_K = Dict(
    :pinned_pinned => 1.0,  :fixed_free => 2.0,
    :fixed_fixed => 0.5,    :fixed_pinned => 0.7
)

"临界应力: σ_cr = P_cr/A = π²E/(KL/r)², r=√(I/A)"
function critical_stress(E::Float64, slenderness::Float64)
    pi^2 * E / slenderness^2
end

"回转半径: r = √(I/A)"
radius_of_gyration(I::Float64, A::Float64) = sqrt(I/A)

"长细比: λ = KL/r"
slenderness_ratio(K::Float64, L::Float64, I::Float64, A::Float64) = K*L / sqrt(I/A)

"屈曲模态形状"
buckling_mode(x::Float64, L::Float64, n=1, A=1.0) = A*sin(n*pi*x/L)

# ============================================================
# 梁振动
# ============================================================
"简支梁简正频率: ω_n = (nπ/L)²√(EI/ρA)"
function beam_natural_frequencies(EI::Float64, rhoA::Float64, L::Float64, nmax=5)
    [(n*pi/L)^2 * sqrt(EI/rhoA) for n in 1:nmax]
end

"梁模态形状 (简支)"
beam_mode_shape(n::Int, x::Float64, L::Float64) = sin(n*pi*x/L)

"悬臂梁基频: ω₁ = 3.516/L² √(EI/ρA)"
cantilever_fundamental_frequency(EI::Float64, rhoA::Float64, L::Float64) = 3.516/L^2*sqrt(EI/rhoA)

# ============================================================
# 不稳定性
# ============================================================
"横向屈曲 (lateral-torsional buckling): M_cr = π/L √(EIy GJ)"
function lateral_torsional_buckling(EIy::Float64, GJ::Float64, L::Float64, Cb=1.0)
    Cb * pi/L * sqrt(EIy*GJ)
end

"圆柱壳屈曲 (轴压): σ_cr = E h/(R√(3(1-ν²)))"
cylinder_axial_buckling_stress(E::Float64, nu::Float64, h::Float64, R::Float64) = E*h/(R*sqrt(3*(1-nu^2)))

"柱 (column) 偏心加载: σ_max = P/A + Pec/I (secant formula)"
function secant_formula(P::Float64, A::Float64, e::Float64, c::Float64, I::Float64, L::Float64, E::Float64)
    r = sqrt(I/A)
    arg = L/(2r) * sqrt(P/(E*A))
    P/A * (1 + e*c/r^2 * 1/cos(arg))
end

# ============================================================
# 连续梁 (Three-Moment Equation / Clapeyron)
# ============================================================
"连续梁 三弯矩方程: M_i*L_i/I_i + 2M_{i+1}*(L_i/I_i+L_{i+1}/I_{i+1}) + M_{i+2}*L_{i+1}/I_{i+1} = -6*load_term"
function three_moment_coefficients(L::Vector{Float64}, I::Vector{Float64})
    n = length(L) - 1
    A = zeros(n, n); b = zeros(n)
    for i in 1:n
        if i > 1; A[i,i-1] = L[i]/I[i]; end
        A[i,i] = 2*(L[i]/I[i] + L[i+1]/I[i+1])
        if i < n; A[i,i+1] = L[i+1]/I[i+1]; end
    end
    A, b
end

"Winkler弹性地基梁特征长度: l_c = (4EI/k)^{1/4}"
function winkler_characteristic_length(EI::Float64, k::Float64)
    (4*EI/k)^(1/4)
end

"Winkler地基梁挠度 (无限长, 集中力): w(x) = P*lambda/(2k)*exp(-lambda*|x|)*(cos+sin)"
function winkler_deflection_point(EI::Float64, k::Float64, P::Float64, x::Float64)
    lambda = 1/winkler_characteristic_length(EI, k)
    P*lambda/(2*k) * exp(-lambda*abs(x)) * (cos(lambda*abs(x)) + sin(lambda*abs(x)))
end

"曲线梁应力 (Winkler-Bach公式): sigma = M*y/(A*e*(r_n-y))"
function curved_beam_stress(M::Float64, A::Float64, e::Float64, r_n::Float64, y::Float64)
    M*y/(A*e*(r_n - y))
end

"剪心位置 (槽钢): e = b^2*h^2*t_f/(4*I_x)"
function shear_center_channel(b::Float64, h::Float64, t_f::Float64, t_w::Float64)
    I_x = b*h^2*t_f/2 + h^3*t_w/12  # 近似
    b^2*h^2*t_f/(4*I_x)
end

"薄壁截面 Saint-Venant 扭转常数: J = sum(b_i*t_i^3/3)"
function torsion_constant_thin_wall(segments::Vector{Tuple{Float64,Float64}})
    J = 0.0
    for (b, t) in segments; J += b*t^3/3; end
    J
end

"翘曲常数 (工字钢): C_w = I_y*h^2/4"
function warping_constant_I_beam(I_y::Float64, h::Float64)
    I_y * h^2 / 4
end

"弯扭屈曲临界载荷: P_cr = (P_y*P_phi + ...)"
function flexural_torsional_buckling(P_y::Float64, P_phi::Float64, e_0::Float64, r_0::Float64)
    # 对称截面: P_cr = min(P_y, P_phi)
    min(P_y, P_phi)
end

"Southwell图法: delta/P vs delta 线性拟合提取P_cr"
function southwell_fit(delta::Vector{Float64}, P::Vector{Float64})
    y = delta ./ P
    # 线性拟合: y = a*delta + b, P_cr = 1/a
    n = length(delta)
    sx = sum(delta); sy = sum(y)
    sxx = sum(delta.^2); sxy = sum(delta .* y)
    a = (n*sxy - sx*sy)/(n*sxx - sx^2)
    b = (sy - a*sx)/n
    P_cr = 1/a
    P_cr, a, b
end

"柱的Perry-Robertson公式 (初始缺陷): sigma_max = 0.5*sigma_cr*(1+eta+sqrt((1+eta)^2-4*sigma_cr/sigma_y))"
function perry_robertson(sigma_cr::Float64, sigma_y::Float64, eta::Float64)
    0.5*sigma_cr*(1+eta+sqrt((1+eta)^2 - 4*sigma_cr/sigma_y))
end

export euler_bernoulli_deflection, beam_max_deflection_uniform, beam_max_deflection_point
export cantilever_deflection, cantilever_uniform, beam_moment, beam_shear_stress
export timoshenko_deflection_point, shear_deformation_ratio
export euler_buckling_load, BUCKLING_K, critical_stress
export radius_of_gyration, slenderness_ratio, buckling_mode
export beam_natural_frequencies, beam_mode_shape, cantilever_fundamental_frequency
export lateral_torsional_buckling, cylinder_axial_buckling_stress, secant_formula
export three_moment_coefficients, winkler_characteristic_length, winkler_deflection_point
export curved_beam_stress, shear_center_channel, torsion_constant_thin_wall
export warping_constant_I_beam, flexural_torsional_buckling, southwell_fit, perry_robertson
