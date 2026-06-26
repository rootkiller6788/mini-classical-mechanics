# forces.jl — 力定律库
# 参考：MIT 8.012 Ch.3-4, Goldstein Ch.1-2

## ============================================================
## 万有引力
## ============================================================

"""
牛顿万有引力定律
F = -G*M*m / r² * r̂
返回作用在位置 r_obj 上的力（由位置 r_source、质量 M 的源产生）
"""
function newton_gravity(r_obj::Vec3, r_source::Vec3, M::Float64, m::Float64; G_val::Float64=G)
    dr = r_obj - r_source
    dist = norm(dr)
    if dist < 1e-10
        return Vec3()  # 避免奇点
    end
    return -G_val * M * m / dist^2 * normalize(dr)
end

"""
多个引力源叠加
sources: [(r_i, M_i), ...]  源的位置和质量列表
返回作用在 r_obj、质量 m 上的合力
"""
function gravity_nbody(r_obj::Vec3, m::Float64, sources::Vector{Tuple{Vec3,Float64}}; G_val::Float64=G)
    f_total = Vec3()
    for (r_src, M_src) in sources
        f_total = f_total + newton_gravity(r_obj, r_src, M_src, m; G_val=G_val)
    end
    return f_total
end

"""
均匀重力场（近地表近似）
F = m * g_vec（默认 -z 方向）
"""
function uniform_gravity(m::Float64; g::Float64=G_EARTH, direction::Vec3=Vec3(0,0,-1))
    return m * g * normalize(direction)
end

"""
引力加速度：a = F/m
"""
function gravity_acceleration(r_obj::Vec3, r_source::Vec3, M::Float64; G_val::Float64=G)
    dr = r_obj - r_source
    dist = norm(dr)
    if dist < 1e-10
        return Vec3()
    end
    return -G_val * M / dist^2 * normalize(dr)
end

## ============================================================
## 弹力（弹簧力）
## ============================================================

"""
胡克定律：F = -k*(r - r_eq)
r_eq 为平衡位置
"""
function hooke_force(r::Vec3, r_eq::Vec3, k::Float64)
    return -k * (r - r_eq)
end

"""
一维弹簧力
"""
function hooke_force_1d(x::Float64, x_eq::Float64, k::Float64)
    return -k * (x - x_eq)
end

## ============================================================
## 阻力/阻尼
## ============================================================

"""
线性阻力（Stokes阻力，适用于低速/小雷诺数）
F_drag = -b * v
"""
function linear_drag(v::Vec3, b::Float64)
    return -b * v
end

"""
平方阻力（适用于高速/大雷诺数）
F_drag = -c * |v| * v
"""
function quadratic_drag(v::Vec3, c::Float64)
    speed = norm(v)
    return -c * speed * v
end

"""
通用阻力模型：F = -(b*v + c*|v|*v)
"""
function combined_drag(v::Vec3, b::Float64, c::Float64)
    speed = norm(v)
    return -(b + c * speed) * v
end

## ============================================================
## 摩擦力
## ============================================================

"""
库仑摩擦力（滑动摩擦）
|f| = μ_k * |N|, 方向与速度相反
N 为法向力大小
"""
function kinetic_friction(v::Vec3, normal_force_magnitude::Float64, mu_k::Float64)
    speed = norm(v)
    if speed < 1e-10
        return Vec3()
    end
    return -mu_k * normal_force_magnitude * (v / speed)
end

"""
静摩擦判断
如果外力切向分量 <= μ_s * N，静摩擦力 = -外力切向分量
否则物体开始滑动
"""
function static_friction_max(normal_force_magnitude::Float64, mu_s::Float64)
    return mu_s * normal_force_magnitude
end

## ============================================================
## 洛伦兹力（电磁力，为后续电磁学模块预留接口）
## ============================================================

"""
洛伦兹力：F = q(E + v × B)
"""
function lorentz_force(q::Float64, E::Vec3, v::Vec3, B::Vec3)
    return q * (E + cross(v, B))
end

"""
纯磁场中带电粒子运动（回旋运动）
"""
function cyclotron_frequency(q::Float64, B_magnitude::Float64, m::Float64)
    return abs(q) * B_magnitude / m
end

function cyclotron_radius(v_perp::Float64, q::Float64, B_magnitude::Float64, m::Float64)
    return m * v_perp / (abs(q) * B_magnitude)
end

## ============================================================
## 约束力（法向力、张力）
## ============================================================

"""
斜面上的法向力
mg cos(theta), theta 为斜面倾角
"""
function normal_force_incline(m::Float64, theta::Float64; g::Float64=G_EARTH)
    return m * g * cos(theta)
end

"""
圆锥摆张力
水平分量提供向心力，垂直分量平衡重力
"""
function conical_pendulum_tension(m::Float64, v::Float64, r::Float64, phi::Float64; g::Float64=G_EARTH)
    # phi 是绳与竖直方向的夹角
    return m * g / cos(phi)
end

## ============================================================
## 常用力的组合
## ============================================================

"""
阻尼谐振子合力
F = -k*(r - r_eq) - b*v
"""
function damped_harmonic_force(r::Vec3, r_eq::Vec3, v::Vec3, k::Float64, b::Float64)
    return hooke_force(r, r_eq, k) + linear_drag(v, b)
end

"""
受迫阻尼谐振子
F = -k*(r - r_eq) - b*v + F_drive(t)
"""
function driven_damped_force(r::Vec3, r_eq::Vec3, v::Vec3, k::Float64, b::Float64,
                             driving_force::Function, t::Float64)
    return hooke_force(r, r_eq, k) + linear_drag(v, b) + driving_force(t)
end
