# constraints.jl — 约束系统：滑轮、斜面、张力
# 参考：MIT 8.012 Ch.7, Kleppner & Kolenkow Ch.3

## ============================================================
## 斜面
## ============================================================

"""
斜面上物体的加速度（沿斜面方向）
theta: 斜面倾角
mu: 动摩擦因数（0 = 无摩擦）
返回沿斜面向下的加速度大小
"""
function incline_acceleration(theta::Float64; mu::Float64=0.0, g::Float64=G_EARTH)
    return g * (sin(theta) - mu * cos(theta))
end

"""
斜面上物体不下滑的临界角（静摩擦）
mu_s: 静摩擦因数
"""
function angle_of_repose(mu_s::Float64)
    return atan(mu_s)
end

"""
斜面上的法向力
"""
function normal_force_incline(m::Float64, theta::Float64; g::Float64=G_EARTH)
    return m * g * cos(theta)
end

"""
斜面上物体的运动（数值模拟）
返回沿斜面方向的加速度矢量（3D，斜面在 xz 平面，沿 +x 向上）
"""
function incline_accel_3d(r::Vec3, v::Vec3, t::Float64, theta::Float64; mu::Float64=0.0, g::Float64=G_EARTH)
    # 斜面法向：(-sinθ, 0, cosθ)，沿斜面向下：(cosθ, 0, sinθ)
    g_vec = Vec3(0.0, 0.0, -g)
    normal = Vec3(-sin(theta), 0.0, cos(theta))
    # 法向分量被约束力抵消
    g_normal = dot(g_vec, normal) * normal
    g_tangential = g_vec - g_normal
    # 摩擦力沿速度反方向
    friction = Vec3()
    spd = norm(v)
    if spd > 1e-10
        friction = -mu * m * abs(dot(g_vec, normal)) * (v / spd)
    end
    return g_tangential + friction
end

## ============================================================
## 滑轮系统
## ============================================================

"""
Atwood 机（两个质量通过轻绳跨过无摩擦滑轮）
返回 (加速度 a₁, 加速度 a₂, 绳子张力 T)
a₁ > 0 表示 m₁ 向下加速
"""
function atwood_machine(m1::Float64, m2::Float64; g::Float64=G_EARTH)
    total_mass = m1 + m2
    a = (m1 - m2) * g / total_mass  # m1 向下为正
    T = 2 * m1 * m2 * g / total_mass
    return a, -a, T  # a1 = a, a2 = -a
end

"""
多滑轮系统（n 段绳，各段力相等）
动滑轮组：n 段绳 → 拉力放大 n 倍
"""
function pulley_system_mechanical_advantage(n_segments::Int, load_mass::Float64; g::Float64=G_EARTH)
    effort = load_mass * g / n_segments
    return effort
end

"""
滑轮系统的加速度（给定拉力 F_pull）
负载质量 m_load，拉力 F_pull，滑轮组绳段数 n
返回负载向上的加速度
"""
function pulley_acceleration(m_load::Float64, F_pull::Float64, n_segments::Int=1; g::Float64=G_EARTH)
    F_net = n_segments * F_pull - m_load * g
    return F_net / m_load
end

## ============================================================
## 张力分析
## ============================================================

"""
轻绳中的张力（质量忽略不计）
两端受力相等，张力沿绳方向
"""
function rope_tension(attached_mass::Float64, acceleration::Vec3, other_forces::Vec3=Vec3(); g::Float64=G_EARTH)
    # T + other_forces + mg = ma → T = ma - mg - other_forces
    g_vec = Vec3(0.0, 0.0, -g)
    return attached_mass * acceleration - attached_mass * g_vec - other_forces
end

"""
圆锥摆（conical pendulum）
绳长 L，摆锤质量 m，角速度 omega
返回 (绳与竖直方向夹角, 张力大小)
"""
function conical_pendulum_parameters(m::Float64, L::Float64, omega::Float64; g::Float64=G_EARTH)
    # T*cosθ = mg, T*sinθ = m*ω²*L*sinθ
    # → cosθ = g/(ω²*L)
    cos_theta = g / (omega^2 * L)
    if cos_theta > 1.0 || cos_theta < 0.0
        return NaN, NaN  # 角速度不够大
    end
    theta = acos(cos_theta)
    T = m * g / cos_theta
    return theta, T
end

"""
圆锥摆稳定条件：ω > √(g/L)
"""
function conical_pendulum_critical_omega(L::Float64; g::Float64=G_EARTH)
    return sqrt(g / L)
end

## ============================================================
## 曲面约束
## ============================================================

"""
质点在光滑曲面上运动的法向约束力
曲率半径 R，速度 v，质量 m
"""
function surface_normal_force(m::Float64, v::Float64, R::Float64, theta::Float64; g::Float64=G_EARTH)
    # N = mg*cosθ + mv²/R（向心加速度贡献）
    return m * (g * cos(theta) + v^2 / R)
end

"""
物体脱离曲面的条件：N ≤ 0
"""
function detachment_condition(m::Float64, v::Float64, R::Float64, theta::Float64; g::Float64=G_EARTH)
    N = surface_normal_force(m, v, R, theta; g=g)
    return N <= 0.0, N
end

"""
过山车/环形轨道最低速度（顶点）
R: 环半径
"""
function loop_the_loop_min_speed(R::Float64; g::Float64=G_EARTH)
    return sqrt(g * R)  # v²/R = g → 向心力 = 重力
end

## ============================================================
## 斜面 + 滑轮组合（经典习题）
## ============================================================

"""
斜面上的物体通过滑轮连接悬挂物
m1 在倾角 theta 的斜面上（摩擦因数 mu），通过轻绳跨过滑轮连接 m2（竖直悬挂）
返回 (系统加速度, 绳子张力), a>0 表示 m1 沿斜面向上
"""
function incline_pulley_system(m1::Float64, m2::Float64, theta::Float64; mu::Float64=0.0, g::Float64=G_EARTH)
    total = m1 + m2
    # 沿斜面方向：T - m1*g*sinθ - μ*m1*g*cosθ = m1*a
    # 竖直方向：m2*g - T = m2*a
    # → a = (m2*g - m1*g*sinθ - μ*m1*g*cosθ) / (m1+m2)
    a = (m2 - m1*sin(theta) - mu*m1*cos(theta)) * g / total
    T = m2 * (g - a)
    return a, T
end

"""
判断斜面-滑轮系统的运动方向
返回 :m1_up（m1 沿斜面向上）, :m1_down（m1 沿斜面向下）, :static（静止）
"""
function incline_pulley_direction(m1::Float64, m2::Float64, theta::Float64; mu_s::Float64=0.0, g::Float64=G_EARTH)
    # 无摩擦时纯重力比较
    force_m2 = m2 * g
    force_m1_down = m1 * g * sin(theta)

    if abs(force_m2 - force_m1_down) <= mu_s * m1 * g * cos(theta)
        return :static
    elseif force_m2 > force_m1_down
        return :m1_up
    else
        return :m1_down
    end
end

## ============================================================
## 非惯性系中的有效重力
## ============================================================

"""
加速参考系中的有效重力
a_frame: 参考系的加速度
"""
function effective_gravity(a_frame::Vec3; g::Float64=G_EARTH, g_direction::Vec3=Vec3(0,0,-1))
    g_real = g * normalize(g_direction)
    return g_real - a_frame
end

"""
加速上升/下降的电梯中的表观重量
a_elevator > 0 表示向上加速
"""
function apparent_weight_in_elevator(m::Float64, a_elevator::Float64; g::Float64=G_EARTH)
    return m * (g + a_elevator)  # a_elevator>0上升→更重; a_elevator<0下降→更轻
end

"""
车辆在弯道上的最佳倾斜角（无侧向摩擦）
v: 车速, R: 弯道半径
"""
function banked_curve_angle(v::Float64, R::Float64; g::Float64=G_EARTH)
    tan_theta = v^2 / (g * R)
    return atan(tan_theta)
end

"""
弯道上车辆的最大安全速度（给定倾斜角 theta 和摩擦系数 mu）
"""
function max_safe_speed_banked_curve(R::Float64, theta::Float64, mu::Float64; g::Float64=G_EARTH)
    # v_max = √(g*R * (tanθ + μ) / (1 - μ*tanθ))
    tan_theta = tan(theta)
    numerator = tan_theta + mu
    denominator = 1.0 - mu * tan_theta
    if denominator <= 0
        return Inf  # 不会滑出
    end
    return sqrt(g * R * numerator / denominator)
end
