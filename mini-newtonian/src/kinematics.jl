# kinematics.jl — 运动学：位置/速度/加速度的关系与轨迹计算
# 参考：MIT 8.012 Ch.1-2

## ============================================================
## 一维运动学
## ============================================================

"""匀速直线运动：r(t) = r0 + v*t"""
function uniform_motion(r0::Float64, v::Float64, t::Float64)
    return r0 + v * t
end

"""匀加速直线运动：r(t) = r0 + v0*t + 0.5*a*t²"""
function accelerated_motion(r0::Float64, v0::Float64, a::Float64, t::Float64)
    return r0 + v0 * t + 0.5 * a * t^2
end

"""匀加速末速度：v(t) = v0 + a*t"""
function accelerated_velocity(v0::Float64, a::Float64, t::Float64)
    return v0 + a * t
end

"""无时间公式：v² = v0² + 2a(r - r0)"""
function velocity_displacement(v0::Float64, a::Float64, dr::Float64)
    return sqrt(v0^2 + 2 * a * dr)
end

## ============================================================
## 三维运动学
## ============================================================

"""匀速直线运动（3D）"""
function uniform_motion(r0::Vec3, v::Vec3, t::Float64)
    return r0 + v * t
end

"""匀加速运动（3D）"""
function accelerated_motion(r0::Vec3, v0::Vec3, a::Vec3, t::Float64)
    return r0 + v0 * t + 0.5 * a * t^2
end

"""匀加速末速度（3D）"""
function accelerated_velocity(v0::Vec3, a::Vec3, t::Float64)
    return v0 + a * t
end

## ============================================================
## 抛体运动（均匀重力场）
## ============================================================

"""
抛体运动解析解
返回：(r(t), v(t))
g 默认指向 -z 方向
"""
function projectile_motion(r0::Vec3, v0::Vec3, t::Float64; g::Float64=G_EARTH)
    g_vec = Vec3(0.0, 0.0, -g)
    r = r0 + v0 * t + 0.5 * g_vec * t^2
    v = v0 + g_vec * t
    return r, v
end

"""抛体飞行时间（从高度 z0 到 z=0）"""
function projectile_flight_time(v0z::Float64, z0::Float64; g::Float64=G_EARTH)
    # z0 + v0z*t - 0.5*g*t² = 0
    discriminant = v0z^2 + 2 * g * z0
    if discriminant < 0
        return NaN  # 不会落地
    end
    return (v0z + sqrt(discriminant)) / g
end

"""抛体最大射程（初速 v0，角度 theta 弧度）"""
function projectile_range(v0::Float64, theta::Float64; g::Float64=G_EARTH, z0::Float64=0.0)
    v0x = v0 * cos(theta)
    v0z = v0 * sin(theta)
    t_flight = projectile_flight_time(v0z, z0; g=g)
    return v0x * t_flight
end

"""最佳抛射角（给定初速和目标水平距离）"""
function optimal_launch_angle(v0::Float64, range::Float64; g::Float64=G_EARTH)
    sin2theta = range * g / v0^2
    if abs(sin2theta) > 1.0
        return NaN  # 射程不可达
    end
    return 0.5 * asin(sin2theta)
end

## ============================================================
## 圆周运动
## ============================================================

"""匀速圆周运动位置"""
function circular_position(center::Vec3, radius::Float64, omega::Float64, t::Float64;
                           normal::Vec3=Vec3(0,0,1))
    # normal 是旋转轴方向
    axis = normalize(normal)
    # 在垂直于轴的平面内选取两个正交方向
    e1 = abs(axis.z) < 0.999 ? normalize(cross(axis, Vec3(0,0,1))) : normalize(cross(axis, Vec3(1,0,0)))
    e2 = cross(axis, e1)
    return center + radius * (cos(omega * t) * e1 + sin(omega * t) * e2)
end

"""匀速圆周运动速度"""
function circular_velocity(radius::Float64, omega::Float64, t::Float64;
                           normal::Vec3=Vec3(0,0,1))
    axis = normalize(normal)
    e1 = abs(axis.z) < 0.999 ? normalize(cross(axis, Vec3(0,0,1))) : normalize(cross(axis, Vec3(1,0,0)))
    e2 = cross(axis, e1)
    return radius * omega * (-sin(omega * t) * e1 + cos(omega * t) * e2)
end

"""向心加速度大小"""
centripetal_acceleration(v::Float64, radius::Float64) = v^2 / radius
centripetal_acceleration_omega(omega::Float64, radius::Float64) = omega^2 * radius

## ============================================================
## 坐标系变换
## ============================================================

"""
惯性系 → 旋转系变换
旋转系以角速度 omega_vec 绕原点旋转
"""
function inertial_to_rotating(r_inertial::Vec3, v_inertial::Vec3, omega_vec::Vec3)
    # 位置不变（瞬时共原点）
    r_rot = r_inertial
    # 速度变换：v_rot = v_inertial - ω × r
    v_rot = v_inertial - cross(omega_vec, r_inertial)
    return r_rot, v_rot
end

"""
旋转系中的表观加速度（含科里奥利力和离心力）
a_inertial 是真实外力产生的加速度
"""
function rotating_frame_acceleration(a_inertial::Vec3, v_rot::Vec3, r_rot::Vec3, omega_vec::Vec3)
    # a_rot = a_inertial - 2(ω × v_rot) - ω × (ω × r)
    coriolis = -2.0 * cross(omega_vec, v_rot)
    centrifugal = -cross(omega_vec, cross(omega_vec, r_rot))
    return a_inertial + coriolis + centrifugal
end

## ============================================================
## 相对运动（伽利略变换）
## ============================================================

"""伽利略位置变换"""
function galilean_position(r_prime::Vec3, v_frame::Vec3, t::Float64)
    return r_prime + v_frame * t
end

"""伽利略速度变换"""
function galilean_velocity(v_prime::Vec3, v_frame::Vec3)
    return v_prime + v_frame
end

## ============================================================
## 轨迹分析工具
## ============================================================

"""从轨迹中提取某个时间点的状态（线性插值）"""
function interpolate_state(traj::Trajectory, t::Float64)
    if t <= traj.ts[1]
        return traj.positions[1], traj.velocities[1]
    end
    if t >= traj.ts[end]
        return traj.positions[end], traj.velocities[end]
    end
    # 二分查找
    i = searchsortedlast(traj.ts, t)
    if i == length(traj.ts)
        return traj.positions[end], traj.velocities[end]
    end
    frac = (t - traj.ts[i]) / (traj.ts[i+1] - traj.ts[i])
    r = traj.positions[i] + frac * (traj.positions[i+1] - traj.positions[i])
    v = traj.velocities[i] + frac * (traj.velocities[i+1] - traj.velocities[i])
    return r, v
end

"""轨迹曲率半径（二维轨迹）"""
function curvature_radius(positions::Vector{Vec3}, velocities::Vector{Vec3})
    n = length(positions)
    radii = Vector{Float64}(undef, n)
    for i in 1:n
        v = norm(velocities[i])
        if i < n && i > 1
            # 用中心差分估算加速度
            dt = 1.0  # 假设等时间步长
            a_approx = (velocities[i+1] - velocities[i-1]) / (2 * dt)
            a_perp = norm(a_approx - dot(a_approx, velocities[i])/v^2 * velocities[i])
            radii[i] = v > 0 ? v^2 / max(a_perp, 1e-15) : Inf
        else
            radii[i] = Inf
        end
    end
    return radii
end
