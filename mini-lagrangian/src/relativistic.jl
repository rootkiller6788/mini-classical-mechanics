# relativistic.jl — 相对论性拉格朗日量
# 参考：Goldstein Ch.7, Landau Vol.2 Ch.2, MIT 8.033
# L4 Fundamental Laws: relativistic mechanics, E=mc^2
# L8 Advanced Topics: covariant formulation

"""
相对论性自由粒子拉格朗日量.

Manifestly covariant (参数化形式):
  S = -mc ∫ ds = -mc ∫ √(dx_μ dx^μ)
  = -mc ∫ √(c² - v²) dt  (取坐标时 t 为参数)

L(v) = -mc² √(1 - v²/c²)
小速度展开: L ≈ -mc² + ½mv² + ... → 非相对论极限给出正确的动能.
"""
function relativistic_particle_L(m::Float64, v_sq::Float64; c=299792458.0)
    gamma_inv = sqrt(1 - v_sq/c^2)
    return -m*c^2 * gamma_inv
end

"""
相对论性能量:
E = ∂L/∂v · v - L
  = γ m c² = m c² / √(1 - v²/c²)
静止能量: E₀ = m c²
"""
function relativistic_energy(m::Float64, v_sq::Float64; c=299792458.0)
    gamma = 1 / sqrt(1 - v_sq/c^2)
    return gamma * m * c^2
end

"""
相对论性动量:
p = ∂L/∂v = γ m v
量值: |p| = γ m v
"""
function relativistic_momentum(m::Float64, v::Vector{Float64}; c=299792458.0)
    v_sq = dot(v, v)
    gamma = 1 / sqrt(1 - v_sq/c^2)
    return gamma * m .* v
end

"""
相对论性能量-动量关系:
E² = p²c² + m²c⁴
"""
function energy_momentum_relation(p_sq::Float64, m::Float64; c=299792458.0)
    return sqrt(p_sq*c^2 + m^2*c^4)
end

"""
固有时参数化的拉格朗日量 (manifestly covariant):
作用量: S = -mc ∫ dτ (dτ = 固有时)
L_proper = -mc

四速度: u^μ = dx^μ/dτ = (γc, γv)
约束: u_μ u^μ = c² (时间-like)
"""
function proper_time_lagrangian(m::Float64; c=299792458.0)
    return -m*c
end

"""
四速度计算: u^μ = dx^μ/dτ
输入: 三维速度 v = (vx, vy, vz)
输出: (u^0, u^1, u^2, u^3)
"""
function four_velocity(v::Vector{Float64}; c=299792458.0)
    v_sq = dot(v, v)
    gamma = 1 / sqrt(1 - v_sq/c^2)
    return [gamma*c, gamma*v[1], gamma*v[2], gamma*v[3]]
end

"""
四动量: p^μ = m u^μ = (E/c, p)
"""
function four_momentum(m::Float64, v::Vector{Float64}; c=299792458.0)
    v_sq = dot(v, v)
    gamma = 1 / sqrt(1 - v_sq/c^2)
    return [gamma*m*c, gamma*m*v[1], gamma*m*v[2], gamma*m*v[3]]
end

"""
相对论性带电粒子在电磁场中的拉格朗日量:
L = -mc² √(1 - v²/c²) - q φ + q v · A

运动方程 (Lorentz 力):
d/dt (γ m v) = q (E + v × B)
"""
function relativistic_charged_particle_L(m::Float64, q::Float64, v::Vector{Float64}, phi::Float64, A::Vector{Float64}; c=299792458.0)
    v_sq = dot(v, v)
    gamma_inv = sqrt(1 - v_sq/c^2)
    L0 = -m*c^2 * gamma_inv
    L_int = -q*phi + q*dot(v, A)
    return L0 + L_int
end

"""
Lorentz 力 (从相对论性 L 推导):
dp/dt = q(E + v × B)
"""
function lorentz_force(q::Float64, E::Vector{Float64}, B::Vector{Float64}, v::Vector{Float64})
    return q .* (E + cross(v, B))
end

"""
相对论性速度加法 (Lorentz 变换下):
若 S' 相对 S 以速度 V 运动, 粒子在 S' 中的速度为 v',
则在 S 中的速度为:
u_parallel = (v'_parallel + V)/(1 + v'_parallel V/c²)
u_perp = v'_perp / (γ_V (1 + v'_parallel V/c²))
"""
function relativistic_velocity_addition(v_prime::Vector{Float64}, V::Float64; c=299792458.0)
    v_par = v_prime[1]
    v_perp = norm(v_prime[2:3])
    gamma_V = 1/sqrt(1 - V^2/c^2)
    u_par = (v_par + V) / (1 + v_par*V/c^2)
    u_perp = v_perp / (gamma_V * (1 + v_par*V/c^2))
    return u_par, u_perp
end

"""
相对论性多普勒效应:
观测频率 f' = f √((1-β)/(1+β))  (源远离)
其中 β = v/c.
"""
function relativistic_doppler(f0::Float64, v::Float64; c=299792458.0, approaching=true)
    beta = v/c
    if approaching
        return f0 * sqrt((1+beta)/(1-beta))
    else
        return f0 * sqrt((1-beta)/(1+beta))
    end
end

"""
相对论性拉格朗日量的非相对论展开 (小 v/c):
L = -mc² √(1 - v²/c²)
  = -mc² [1 - v²/(2c²) - v⁴/(8c⁴) - O(v⁶/c⁶)]
  = -mc² + ½mv² + mv⁴/(8c²) + ...

第一修正项: mv⁴/(8c²) (与牛顿力学的偏差).
"""
function relativistic_correction_1st(m::Float64, v::Float64; c=299792458.0)
    beta = v/c
    L_newton = 0.5*m*v^2
    L_rel = -m*c^2*sqrt(1-beta^2) + m*c^2  # subtract rest energy
    correction = L_rel - L_newton
    theory_correction = m*v^4/(8*c^2)
    return L_newton, L_rel, correction, theory_correction
end

"""
双生子佯谬: 计算旅行者和地球上的固有时.

旅行者以速度 v 旅行距离 D, 然后掉头返回.
地球固有时: T_earth = 2D/v
旅行者固有时: T_traveler = 2(D/v) √(1 - v²/c²)
"""
function twin_paradox(D::Float64, v::Float64; c=299792458.0)
    T_earth = 2*D/v
    T_traveler = T_earth * sqrt(1 - v^2/c^2)
    return T_earth, T_traveler, T_earth - T_traveler
end

"""
相对论性 Lagrange 形式下的守恒量 (闵氏时空对称性):
- 时空平移 → 四动量守恒
- Lorentz 变换 → 四角动量守恒

闵氏度规: η_{μν} = diag(-1, 1, 1, 1) 或 diag(1, -1, -1, -1).
此处用 (-,+,+,+) 约定.
"""
function minkowski_dot(a::Vector{Float64}, b::Vector{Float64})
    return -a[1]*b[1] + a[2]*b[2] + a[3]*b[3] + a[4]*b[4]
end

"""
Lorentz boost 变换矩阵 (沿 x 轴):
Λ^μ_ν = [ γ    -γβ  0  0
          -γβ   γ   0  0
           0    0   1  0
           0    0   0  1 ]
其中 β = v/c, γ = 1/√(1-β²).
"""
function lorentz_boost_matrix(v::Float64; c=299792458.0)
    beta = v/c
    gamma = 1/sqrt(1-beta^2)
    return [gamma -gamma*beta 0 0; -gamma*beta gamma 0 0; 0 0 1 0; 0 0 0 1]
end

"""
闵氏时空间隔的不变性:
ds² = -c² dt² + dx² + dy² + dz²
对于光信号: ds² = 0 (类光)
对于实物粒子: ds² < 0 (类时)
对于超光速: ds² > 0 (类空)
"""
function spacetime_interval(cdt::Float64, dx::Float64, dy::Float64, dz::Float64)
    return -cdt^2 + dx^2 + dy^2 + dz^2
end
