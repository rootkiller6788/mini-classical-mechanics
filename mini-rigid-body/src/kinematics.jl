# kinematics.jl -- 欧拉角、旋转矩阵与角速度的互转
# 参考: Goldstein Ch.4, Landau Ch.6

using LinearAlgebra

# ============================================================
# 欧拉角 → 旋转矩阵 (ZXZ convention)
# ============================================================

"""
欧拉角 → 3x3 旋转矩阵 (ZXZ convention: R = Rz(ϕ)·Rx(θ)·Rz(ψ))

参数:
- phi: 进动角
- theta: 章动角
- psi: 自转角

返回: 3x3 正交旋转矩阵 R

参考: Goldstein (4.46)
"""
function euler_to_rotation(phi::Float64, theta::Float64, psi::Float64)
    cp, sp = cos(phi), sin(phi)
    ct, st = cos(theta), sin(theta)
    cpp, spp = cos(psi), sin(psi)

    return [
        cp*cpp - sp*ct*spp   -cp*spp - sp*ct*cpp   sp*st;
        sp*cpp + cp*ct*spp   -sp*spp + cp*ct*cpp  -cp*st;
        st*spp                st*cpp               ct
    ]
end

# ============================================================
# 旋转矩阵 → 欧拉角 (ZXZ convention)
# ============================================================

"""
3x3 旋转矩阵 → 欧拉角 (ZXZ convention)

参数:
- R: 3x3 正交旋转矩阵

返回: EulerAngles(phi, theta, psi)

参考: Goldstein (4.47)
"""
function rotation_to_euler(R::AbstractMatrix)
    # theta = arccos(R33)
    theta = acos(clamp(R[3,3], -1.0, 1.0))

    if abs(R[3,3] - 1.0) < 1e-15
        # theta ≈ 0: 仅有进动+自转的总和可确定
        phi = atan(R[1,1], -R[1,2])
        psi = 0.0
    elseif abs(R[3,3] + 1.0) < 1e-15
        # theta ≈ π: 仅进动-自转的差可确定
        phi = atan(-R[1,1], -R[1,2])
        psi = 0.0
    else
        phi = atan(R[1,3], -R[2,3])
        psi = atan(R[3,1], R[3,2])
    end

    return EulerAngles(mod(phi, 2π), theta, mod(psi, 2π))
end

# ============================================================
# 欧拉角速率 ↔ 角速度
# ============================================================

"""
欧拉角速率 → 体坐标系角速度

将 (φ̇, θ̇, ψ̇) 转换为体坐标系角速度分量 (ωx, ωy, ωz)

公式 (Goldstein 4.86-4.88):
  ωx = φ̇ sinθ sinψ + θ̇ cosψ
  ωy = φ̇ sinθ cosψ - θ̇ sinψ
  ωz = φ̇ cosθ + ψ̇
"""
function euler_rates_to_omega(phidot::Float64, thetadot::Float64, psidot::Float64,
                               theta::Float64, psi::Float64)
    st = sin(theta)
    sp = sin(psi)
    cp = cos(psi)
    return [
        phidot * st * sp + thetadot * cp,
        phidot * st * cp - thetadot * sp,
        phidot * cos(theta) + psidot
    ]
end

"""
体坐标系角速度 → 欧拉角速率

从 (ωx, ωy, ωz, θ, ψ) 计算 (φ̇, θ̇, ψ̇)

公式 (Goldstein 4.87 求逆):
  φ̇ = (ωx sinψ + ωy cosψ) / sinθ
  θ̇ = ωx cosψ - ωy sinψ
  ψ̇ = ωz - φ̇ cosθ

注意: θ≈0 时 sinθ→0, φ̇ 发散（万向锁）
"""
function omega_to_euler_rates(w::AbstractVector, theta::Float64, psi::Float64)
    st = sin(theta)
    sp = sin(psi)
    cp = cos(psi)

    if abs(st) < 1e-15
        # 万向锁: φ 和 ψ 无法单独确定，设 ψ̇=0
        phidot = (w[1]*sp + w[2]*cp) / max(st, 1e-15)
        thetadot = w[1]*cp - w[2]*sp
        psidot = 0.0
    else
        phidot = (w[1]*sp + w[2]*cp) / st
        thetadot = w[1]*cp - w[2]*sp
        psidot = w[3] - phidot * cos(theta)
    end

    return (phidot, thetadot, psidot)
end

"""
欧拉角增量 → 旋转矩阵增量

给定欧拉角小变化 (dφ, dθ, dψ)，计算无限小旋转矩阵变化
"""
function euler_increment_to_dR(dphi::Float64, dtheta::Float64, dpsi::Float64,
                                theta::Float64, psi::Float64)
    omega_in_body = euler_rates_to_omega(dphi, dtheta, dpsi, theta, psi)
    # 反对称矩阵 (角速度矩阵)
    wx, wy, wz = omega_in_body
    return [0 -wz wy; wz 0 -wx; -wy wx 0]
end

# ============================================================
# 角速度与旋转矩阵的微分关系
# ============================================================

"""
旋转矩阵的微分方程: dR/dt = R·[ω]ₓ

[ω]ₓ 为角速度 ω 的反对称矩阵 (so(3))

返回: dR/dt = R * omega_cross_matrix
"""
function rotation_derivative(R::AbstractMatrix, w::AbstractVector)
    wx, wy, wz = w
    omega_cross = [0 -wz wy; wz 0 -wx; -wy wx 0]
    return R * omega_cross
end

"""
角速度 → 反对称矩阵 (cross product matrix in so(3))

[ω]ₓ = [[0, -ωz, ωy], [ωz, 0, -ωx], [-ωy, ωx, 0]]
"""
function cross_matrix(w::AbstractVector)
    wx, wy, wz = w
    return [0 -wz wy; wz 0 -wx; -wy wx 0]
end

# ============================================================
# 旋转矩阵的 Cayley-Klein 参数 (可选)
# ============================================================

"""
旋转矩阵 → 四元数 (Hamilton convention)
"""
function rotation_to_quaternion(R::AbstractMatrix)
    trace = R[1,1] + R[2,2] + R[3,3]

    if trace > 0
        s = 0.5 / sqrt(trace + 1.0)
        w = 0.25 / s
        x = (R[3,2] - R[2,3]) * s
        y = (R[1,3] - R[3,1]) * s
        z = (R[2,1] - R[1,2]) * s
    elseif R[1,1] > R[2,2] && R[1,1] > R[3,3]
        s = 2.0 * sqrt(1.0 + R[1,1] - R[2,2] - R[3,3])
        w = (R[3,2] - R[2,3]) / s
        x = 0.25 * s
        y = (R[1,2] + R[2,1]) / s
        z = (R[1,3] + R[3,1]) / s
    elseif R[2,2] > R[3,3]
        s = 2.0 * sqrt(1.0 + R[2,2] - R[1,1] - R[3,3])
        w = (R[1,3] - R[3,1]) / s
        x = (R[1,2] + R[2,1]) / s
        y = 0.25 * s
        z = (R[2,3] + R[3,2]) / s
    else
        s = 2.0 * sqrt(1.0 + R[3,3] - R[1,1] - R[2,2])
        w = (R[2,1] - R[1,2]) / s
        x = (R[1,3] + R[3,1]) / s
        y = (R[2,3] + R[3,2]) / s
        z = 0.25 * s
    end

    return (w, x, y, z)
end

# ============================================================
# 常用旋转操作
# ============================================================

"绕 x 轴旋转 angle 弧度"
function rot_x(angle::Float64)
    c, s = cos(angle), sin(angle)
    return [1 0 0; 0 c -s; 0 s c]
end

"绕 y 轴旋转 angle 弧度"
function rot_y(angle::Float64)
    c, s = cos(angle), sin(angle)
    return [c 0 s; 0 1 0; -s 0 c]
end

"绕 z 轴旋转 angle 弧度"
function rot_z(angle::Float64)
    c, s = cos(angle), sin(angle)
    return [c -s 0; s c 0; 0 0 1]
end

"绕任意轴旋转（Rodrigues 公式）"
function rot_axis(axis::AbstractVector, angle::Float64)
    n = axis / norm(axis)
    c, s = cos(angle), sin(angle)
    nx, ny, nz = n
    return [
        c + nx*nx*(1-c)       nx*ny*(1-c) - nz*s    nx*nz*(1-c) + ny*s;
        ny*nx*(1-c) + nz*s    c + ny*ny*(1-c)       ny*nz*(1-c) - nx*s;
        nz*nx*(1-c) - ny*s    nz*ny*(1-c) + nx*s    c + nz*nz*(1-c)
    ]
end

export euler_to_rotation, rotation_to_euler
export euler_rates_to_omega, omega_to_euler_rates, euler_increment_to_dR
export rotation_derivative, cross_matrix, rotation_to_quaternion
export rot_x, rot_y, rot_z, rot_axis
