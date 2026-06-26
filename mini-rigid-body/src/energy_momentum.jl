# energy_momentum.jl -- 转动动能、角动量、稳定性分析
# 参考: Goldstein Ch.5, Landau Ch.6

using LinearAlgebra

# ============================================================
# 转动动能
# ============================================================

"""
转动动能 T = ½ ωᵀ·I·ω

参数:
- I: 惯性张量 (InertiaTensor 或已对角化)
- w: 角速度向量

公式 (Goldstein 5.9):
  T = ½ (Ixx ωx² + Iyy ωy² + Izz ωz² + 2Ixy ωxωy + 2Ixz ωxωz + 2Iyz ωyωz)

主轴系简化:
  T = ½ (I1 ω₁² + I2 ω₂² + I3 ω₃²)
"""
function rotational_kinetic_energy(I::InertiaTensor, w::AbstractVector)
    return 0.5 * dot(w, inertia_matrix(I) * w)
end

"""
主轴系转动动能（更快，免矩阵乘法）
"""
function rotational_KE_principal(I::InertiaTensor, w::AbstractVector)
    I1, I2, I3 = principal_moments(I)
    return 0.5 * (I1*w[1]^2 + I2*w[2]^2 + I3*w[3]^2)
end

"""
总动能 = 平动动能 + 转动动能

T_total = ½ M v_cm² + ½ ωᵀ·I_cm·ω
"""
function total_kinetic_energy(M::Float64, v_cm::AbstractVector,
                               I_cm::InertiaTensor, w::AbstractVector)
    T_trans = 0.5 * M * dot(v_cm, v_cm)
    T_rot = rotational_kinetic_energy(I_cm, w)
    return T_trans + T_rot, T_trans, T_rot
end

# ============================================================
# 角动量
# ============================================================

"""
刚体角动量 L = I·ω（体坐标系分量）

公式 (Goldstein 5.10):
  Lx = Ixx ωx + Ixy ωy + Ixz ωz
  Ly = Ixy ωx + Iyy ωy + Iyz ωz
  Lz = Ixz ωx + Iyz ωy + Izz ωz
"""
function angular_momentum_rigid(I::InertiaTensor, w::AbstractVector)
    return inertia_matrix(I) * w
end

"""
主轴系角动量（免矩阵乘法）: L = (I1 ω₁, I2 ω₂, I3 ω₃)
"""
function angular_momentum_principal(I::InertiaTensor, w::AbstractVector)
    I1, I2, I3 = principal_moments(I)
    return [I1*w[1], I2*w[2], I3*w[3]]
end

"""
总角动量 = 轨道角动量 + 自转角动量
L_total = r_cm × (M v_cm) + I_cm·ω
"""
function total_angular_momentum(r_cm::AbstractVector, M::Float64, v_cm::AbstractVector,
                                 I_cm::InertiaTensor, w::AbstractVector)
    L_orbital = M * LinearAlgebra.cross(r_cm, v_cm)
    L_spin = angular_momentum_rigid(I_cm, w)
    return L_orbital + L_spin, L_orbital, L_spin
end

# ============================================================
# 中间轴定理 (Tennis Racket Theorem / Dzhanibekov Effect)
# ============================================================

"""
中间轴不稳定性分析

分析绕主轴旋转的线性稳定性（小扰动线性化）

设 ω = ω₀ e_i + δω
- I1 > I2 > I3 (I2 为中间轴)
- 绕 I1 旋转: 稳定
- 绕 I2 旋转: 不稳定（指数发散）
- 绕 I3 旋转: 稳定

参数:
- I: 惯性张量
- axis: 要检测的主轴编号 (1, 2, 3)

返回: (:stable, growth_rate) 或 (:unstable, growth_rate)

公式: 绕主轴 e_k 旋转时，线性化系统 d(δω)/dt = A·δω
A 的特征值 λ = sqrt((I_k - I_i)(I_j - I_k) / (I_i I_j)) * ω₀

参考: Goldstein 5.46, Landau Ch.37
"""
function axis_stability(I::InertiaTensor, axis::Int, omega_mag::Float64=1.0)
    I1, I2, I3 = principal_moments(I)

    if axis == 1
        rate = omega_mag * sqrt(abs((I1 - I2)*(I3 - I1) / (I2 * I3)))
        return (:stable, 0.0)  # I1 最大 → 稳定
    elseif axis == 2
        rate = omega_mag * sqrt(abs((I2 - I1)*(I3 - I2) / (I1 * I3)))
        return (:unstable, rate)
    elseif axis == 3
        rate = omega_mag * sqrt(abs((I3 - I1)*(I2 - I3) / (I1 * I2)))
        return (:stable, 0.0)  # I3 最小 → 稳定
    else
        error("axis must be 1, 2, or 3")
    end
end

"""
稳定性分析字典（三轴通用判定）

返回: Dict(:I1=>"stable"/"unstable", :I2=>..., :I3=>...)

中间轴不稳定准则:
- 如果 I1 ≥ I2 ≥ I3 且 I2 严格位于中间 → 绕 I2 旋转不稳定
- 对称陀螺 (I1 = I2 或 I2 = I3) → 退化情况
"""
function stability_analysis(I::InertiaTensor)
    I1, I2, I3 = principal_moments(I)

    tol = 1e-12
    s1 = "stable"   # 最大轴总是稳定
    s3 = "stable"   # 最小轴总是稳定

    # I2 判定
    if abs(I1 - I2) < tol || abs(I2 - I3) < tol
        s2 = "degenerate"  # 对称陀螺
    elseif I2 > I1 || I2 < I3
        s2 = "degenerate"  # 内部检验
    else
        s2 = "unstable"    # 严格中间轴 → 不稳定
    end

    return Dict(:I1 => s1, :I2 => s2, :I3 => s3)
end

"""
网球拍定理演示: 绕中间轴每半周期翻转

分析绕中间轴的周期:
  T_flip ≈ log(1/δω₀) / growth_rate

其中 δω₀ 为初始扰动大小
"""
function flipping_period_estimate(I::InertiaTensor, omega_mag::Float64, delta_omega0::Float64)
    I1, I2, I3 = principal_moments(I)
    rate = axis_stability(I, 2, omega_mag)[2]
    if rate < 1e-15
        return Inf
    end
    return log(1.0 / delta_omega0) / rate
end

# ============================================================
# 动能椭球与角动量球
# ============================================================

"""
在单位角动量球上采样动能（产生 Binet 椭球交线）

返回: (ω1_grid, ω2_grid, ω3_grid, T_on_sphere)

T = L²/(2I3) 最大 (绕最小惯量轴), T = L²/(2I1) 最小 (绕最大惯量轴)
"""
function energy_on_momentum_sphere(I::InertiaTensor, L2::Float64; n_pts::Int=100)
    I1, I2, I3 = principal_moments(I)
    thetas = range(0, π, length=n_pts)
    phis = range(0, 2π, length=n_pts)
    T_map = zeros(n_pts, n_pts)

    for (i, th) in enumerate(thetas)
        st = sin(th)
        for (j, ph) in enumerate(phis)
            # ω 在球面上的参数化
            ω1 = sin(th) * cos(ph)
            ω2 = sin(th) * sin(ph)
            ω3 = cos(th)
            # 缩放到角动量球面: I1² ω1² + I2² ω2² + I3² ω3² = L2
            denom = sqrt((I1*ω1)^2 + (I2*ω2)^2 + (I3*ω3)^2)
            scale = sqrt(L2) / max(denom, 1e-15)
            ω1 *= scale; ω2 *= scale; ω3 *= scale
            T_map[i, j] = 0.5 * (I1*ω1^2 + I2*ω2^2 + I3*ω3^2)
        end
    end

    return thetas, phis, T_map
end

# ============================================================
# 力矩-功率关系
# ============================================================

"""
力矩做的功: dW = N·ω dt

参数:
- N: 力矩向量
- w: 角速度向量
- dt: 时间间隔
"""
function torque_power(N::AbstractVector, w::AbstractVector)
    return dot(N, w)
end

"""
力矩冲量: ∫N dt

用于估计角动量的变化
"""
function torque_impulse(N::AbstractVector, dt::Float64)
    return [N[1]*dt, N[2]*dt, N[3]*dt]
end

export rotational_kinetic_energy, rotational_KE_principal, total_kinetic_energy
export angular_momentum_rigid, angular_momentum_principal, total_angular_momentum
export axis_stability, stability_analysis, flipping_period_estimate
export energy_on_momentum_sphere, torque_power, torque_impulse
