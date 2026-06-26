# integrable_systems.jl — 可积性与 Arnold-Liouville 定理
# 参考：Arnold Ch.10, Goldstein Ch.10

"""
可积系统：n 自由度系统，存在 n 个相互对合的运动积分
F_1=H, F_2, ..., F_n，满足 {F_i, F_j} = 0

Arnold-Liouville 定理：紧致连通能级集同胚于 n 维环面 T^n，
运动在环面上是准周期的。
"""

"""
检测给定的一组函数是否两两对合（相互 Poisson 可交换）
"""
function are_in_involution(fs::Vector{Function}, q::Vector{Float64}, p::Vector{Float64}; tol=1e-8)
    m = length(fs)
    for i in 1:m, j in i+1:m
        pb = poisson_bracket(fs[i], fs[j], q, p)
        if abs(pb) > tol
            return false
        end
    end
    return true
end

"""
Kepler 问题的运动积分:
  H = p²/(2m) - k/r (能量)
  L = r×p (角动量矢量，3个分量)
  A = p×L/m - k*r/r (Runge-Lenz 矢量)

其中只有 5 个独立积分（6 个分量 - 2 个约束：L·A=0, A²=2mEL²+k²）
"""

struct KeplerProblem
    m::Float64; k::Float64  # k = GMm
end

function kepler_energy(kp::KeplerProblem, r::Vector{Float64}, p::Vector{Float64})
    r_mag = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
    return sum(p.^2)/(2*kp.m) - kp.k/r_mag
end

function kepler_angular_momentum(r::Vector{Float64}, p::Vector{Float64})
    return [
        r[2]*p[3] - r[3]*p[2],
        r[3]*p[1] - r[1]*p[3],
        r[1]*p[2] - r[2]*p[1]
    ]
end

function runge_lenz_vector(kp::KeplerProblem, r::Vector{Float64}, p::Vector{Float64})
    L = kepler_angular_momentum(r, p)
    r_mag = sqrt(sum(r.^2))
    # A = p×L/m - k*r̂
    p_cross_L = [
        p[2]*L[3] - p[3]*L[2],
        p[3]*L[1] - p[1]*L[3],
        p[1]*L[2] - p[2]*L[1]
    ]
    A = p_cross_L / kp.m - kp.k * r / r_mag
    return A
end

"""
Kepler 轨道参数（从守恒量推导）:
  离心率: e = |A|/k
  半长轴: a = k/(2|E|) (椭圆 E<0)
  半通径: p_semi = L²/(mk)
"""
function kepler_orbit_params(kp::KeplerProblem, r::Vector{Float64}, p::Vector{Float64})
    L_vec = kepler_angular_momentum(r, p)
    L2 = sum(L_vec.^2)
    E = kepler_energy(kp, r, p)
    A_vec = runge_lenz_vector(kp, r, p)
    A_mag = sqrt(sum(A_vec.^2))
    e = A_mag / kp.k
    # 半长轴（仅对椭圆 E<0 有意义）
    a = E < 0 ? kp.k / (2*abs(E)) : Inf
    p_semi = L2 / (kp.m * kp.k)
    return (eccentricity=e, semi_major=a, semi_latus_rectum=p_semi, L2=L2, E=E, A=A_vec)
end

"""
检测 Kepler 问题中 n=3 自由度的 5 个独立积分是否对合

注意：L 的三个分量不对合（{L_x, L_y}=L_z），但 L² 与所有分量对合。
可积性来自: H, L², L_z 三个对合积分 + Kepler 特有的超可积性(A的存在)。
"""
function verify_kepler_integrals(kp::KeplerProblem, r::Vector{Float64}, p::Vector{Float64})
    H_func(qq,pp) = kepler_energy(kp, qq, pp)
    L2_func(qq,pp) = sum(kepler_angular_momentum(qq, pp).^2)
    Lz_func(qq,pp) = kepler_angular_momentum(qq, pp)[3]
    Az_func(qq,pp) = runge_lenz_vector(kp, qq, pp)[3]

    f_list = [H_func, L2_func, Lz_func]
    involutive = are_in_involution(f_list, r, p)
    
    println("Kepler integrals at r=$r, p=$p:")
    println("  H  = $(round(H_func(r,p),digits=4))")
    println("  L² = $(round(L2_func(r,p),digits=4))")
    println("  Lz = $(round(Lz_func(r,p),digits=4))")
    println("  Az = $(round(Az_func(r,p),digits=4))")
    println("  {H,L²}={H,Lz}={L²,Lz}=0: $involutive")
    
    return involutive
end
