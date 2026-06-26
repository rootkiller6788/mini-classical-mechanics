# momentum.jl — 动量、冲量、碰撞
# 参考：MIT 8.012 Ch.6, Goldstein Ch.2

## ============================================================
## 单粒子动量
## ============================================================

"""线动量: p = m*v"""
linear_momentum(m::Float64, v::Vec3) = m * v

"""角动量 (绕原点): L = r × p = m(r × v)"""
function angular_momentum(r::Vec3, v::Vec3, m::Float64)
    return m * cross(r, v)
end

"""角动量（绕任意点 origin）"""
function angular_momentum_about(r::Vec3, v::Vec3, m::Float64, origin::Vec3)
    return m * cross(r - origin, v)
end

## ============================================================
## 系统动量
## ============================================================

"""N 体系统总线动量"""
function total_linear_momentum(masses::Vector{Float64}, velocities::Vector{Vec3})
    p_total = Vec3()
    for i in 1:length(masses)
        p_total = p_total + masses[i] * velocities[i]
    end
    return p_total
end

"""N 体系统总角动量（绕原点）"""
function total_angular_momentum(masses::Vector{Float64}, positions::Vector{Vec3}, velocities::Vector{Vec3})
    L_total = Vec3()
    for i in 1:length(masses)
        L_total = L_total + angular_momentum(positions[i], velocities[i], masses[i])
    end
    return L_total
end

"""质心位置"""
function center_of_mass(masses::Vector{Float64}, positions::Vector{Vec3})
    M_total = sum(masses)
    cm = Vec3()
    for i in 1:length(masses)
        cm = cm + masses[i] * positions[i]
    end
    return cm / M_total
end

"""质心速度"""
function center_of_mass_velocity(masses::Vector{Float64}, velocities::Vector{Vec3})
    M_total = sum(masses)
    v_cm = Vec3()
    for i in 1:length(masses)
        v_cm = v_cm + masses[i] * velocities[i]
    end
    return v_cm / M_total
end

"""约化质量: μ = m1*m2 / (m1+m2)"""
reduced_mass(m1::Float64, m2::Float64) = m1 * m2 / (m1 + m2)

## ============================================================
## 冲量
## ============================================================

"""冲量 J = F*Δt（恒力）"""
impulse_constant(F::Vec3, dt::Float64) = F * dt

"""冲量 = Δp = m*(v_f - v_i)（动量变化）"""
impulse_from_momentum(m::Float64, v_i::Vec3, v_f::Vec3) = m * (v_f - v_i)

"""平均力 = 冲量 / 时间"""
average_force_from_impulse(J::Vec3, dt::Float64) = J / dt

"""冲量-动量定理验证"""
function check_impulse_momentum(m::Float64, v_before::Vec3, v_after::Vec3, J::Vec3; tol::Float64=1e-10)
    delta_p = m * (v_after - v_before)
    diff = norm(delta_p - J)
    return diff < tol, diff
end

## ============================================================
## 一维碰撞
## ============================================================

"""
一维完全弹性碰撞
返回 (v1_final, v2_final)
"""
function elastic_collision_1d(m1::Float64, v1::Float64, m2::Float64, v2::Float64)
    v1f = (m1 - m2) / (m1 + m2) * v1 + 2 * m2 / (m1 + m2) * v2
    v2f = 2 * m1 / (m1 + m2) * v1 + (m2 - m1) / (m1 + m2) * v2
    return v1f, v2f
end

"""
一维完全非弹性碰撞（粘在一起）
返回共同末速度
"""
function inelastic_collision_1d(m1::Float64, v1::Float64, m2::Float64, v2::Float64)
    v_common = (m1 * v1 + m2 * v2) / (m1 + m2)
    return v_common
end

"""
碰撞恢复系数: e = |v2f - v1f| / |v2i - v1i|
e=1 弹性, e=0 完全非弹性
"""
function coefficient_of_restitution(v1i::Float64, v2i::Float64, v1f::Float64, v2f::Float64)
    return abs(v2f - v1f) / abs(v2i - v1i)
end

"""
给定恢复系数 e 的一维碰撞
"""
function collision_with_restitution(m1::Float64, v1::Float64, m2::Float64, v2::Float64, e::Float64)
    v_cm = (m1 * v1 + m2 * v2) / (m1 + m2)
    # 在质心系中速度乘以 -e
    v1_cm = v1 - v_cm
    v2_cm = v2 - v_cm
    v1f_cm = -e * v1_cm
    v2f_cm = -e * v2_cm
    return v1f_cm + v_cm, v2f_cm + v_cm
end

## ============================================================
## 三维碰撞
## ============================================================

"""
三维完全弹性碰撞（光滑球）
假定碰撞时沿法线方向交换动量，切向不变
normal: 碰撞法线方向（从球1中心指向球2中心）
"""
function elastic_collision_3d(m1::Float64, v1::Vec3, m2::Float64, v2::Vec3, normal::Vec3)
    n = normalize(normal)
    # 分解为法向和切向分量
    v1n = dot(v1, n)
    v1t = v1 - v1n * n
    v2n = dot(v2, n)
    v2t = v2 - v2n * n

    # 法向弹性碰撞
    v1nf, v2nf = elastic_collision_1d(m1, v1n, m2, v2n)

    return v1t + v1nf * n, v2t + v2nf * n
end

## ============================================================
## 火箭方程（变质量系统）
## ============================================================

"""
齐奥尔科夫斯基火箭方程
Δv = ve * ln(m0 / mf)
ve: 排气速度, m0: 初始质量, mf: 最终质量
"""
function tsiolkovsky_delta_v(ve::Float64, m0::Float64, mf::Float64)
    return ve * log(m0 / mf)
end

"""
给定 Δv 需要的质量比
"""
function mass_ratio_for_delta_v(delta_v::Float64, ve::Float64)
    return exp(delta_v / ve)
end
