# mission.jl — 轨道转移: Hohmann, 双椭圆, 平面变化, 引力辅助
# 参考: Vallado Ch.6, Bate Mueller White Ch.6

# ============================================================
# Hohmann 转移
# ============================================================

"""
Hohmann 转移轨道

从圆轨道 r1 转移到同平面圆轨道 r2 (r1 < r2, 向外转移)

转移轨道半长径: a_trans = (r1 + r2)/2
转移轨道离心率: e_trans = (r2 - r1)/(r2 + r1)

参数:
- r1: 初始轨道半径
- r2: 目标轨道半径
- mu: 引力参数

返回: (delta_v1, delta_v2, delta_v_total, transfer_time, a_trans)
"""
function hohmann_transfer(r1::Float64, r2::Float64; mu::Float64=G_SUN)
    v1 = sqrt(mu / r1)
    v2 = sqrt(mu / r2)

    a_trans = (r1 + r2) / 2.0
    e_trans = abs(r2 - r1) / (r2 + r1)

    # 转移轨道近心点速度
    v_trans_peri = sqrt(mu * (2.0/r1 - 1.0/a_trans))
    delta_v1 = abs(v_trans_peri - v1)

    # 转移轨道远心点速度
    v_trans_apo = sqrt(mu * (2.0/r2 - 1.0/a_trans))
    delta_v2 = abs(v2 - v_trans_apo)

    delta_v_total = delta_v1 + delta_v2
    transfer_time = π * sqrt(a_trans^3 / mu)  # 半个椭圆周期

    return (delta_v1=delta_v1, delta_v2=delta_v2,
            delta_v_total=delta_v_total, transfer_time=transfer_time,
            a_trans=a_trans, e_trans=e_trans)
end

"""
Hohmann 向内转移 (r1 > r2)
"""
function hohmann_transfer_inward(r1::Float64, r2::Float64; mu::Float64=G_SUN)
    # Hohmann 在数学上是对称的
    result = hohmann_transfer(r2, r1; mu=mu)
    return (delta_v1=result.delta_v2, delta_v2=result.delta_v1,
            delta_v_total=result.delta_v_total, transfer_time=result.transfer_time,
            a_trans=result.a_trans, e_trans=result.e_trans)
end

# ============================================================
# 双椭圆转移 (Bi-elliptic Transfer)
# ============================================================

"""
双椭圆转移

从圆轨道 r1 到 r2，经由中间点 r_b (> r1, r2)

步骤:
  1. 从 r1 到 r_b 的椭圆 (Δv₁)
  2. 在 r_b 变轨到从 r_b 到 r2 的椭圆 (Δv₂)
  3. 在 r2 圆化 (Δv₃)

参数:
- r1, r2: 起始和目标轨道半径
- r_b: 中间点半径 (> max(r1, r2))
- mu: 引力参数

返回: (Δv1, Δv2, Δv3, total_Δv, transfer_time)
"""
function bi_elliptic_transfer(r1::Float64, r2::Float64, r_b::Float64; mu::Float64=G_SUN)
    v1 = sqrt(mu / r1)
    v2 = sqrt(mu / r2)

    # 第一个椭圆: r1 → r_b
    a1 = (r1 + r_b) / 2.0
    v_peri1 = sqrt(mu * (2.0/r1 - 1.0/a1))
    v_apo1 = sqrt(mu * (2.0/r_b - 1.0/a1))
    delta_v1 = abs(v_peri1 - v1)

    # 第二个椭圆: r_b → r2
    a2 = (r_b + r2) / 2.0
    v_peri2 = sqrt(mu * (2.0/r_b - 1.0/a2))
    v_apo2 = sqrt(mu * (2.0/r2 - 1.0/a2))
    delta_v2 = abs(v_peri2 - v_apo1)

    # 圆化
    delta_v3 = abs(v2 - v_apo2)

    delta_v_total = delta_v1 + delta_v2 + delta_v3
    transfer_time = π * (sqrt(a1^3 / mu) + sqrt(a2^3 / mu))

    return (delta_v1=delta_v1, delta_v2=delta_v2, delta_v3=delta_v3,
            delta_v_total=delta_v_total, transfer_time=transfer_time)
end

# ============================================================
# 平面变化
# ============================================================

"""
简单平面变化 (在圆轨道上)

Δv = 2 v_circ sin(Δi/2)

参数:
- v: 轨道速度
- delta_i: 倾角变化 [rad]
"""
simple_plane_change_delta_v(v::Float64, delta_i::Float64) = 2.0 * v * sin(delta_i / 2.0)

"""
复合平面变化 (在远心点做平面变化可省 Δv)

Δv = 2 v_apo sin(Δi/2)
"""
function apoapsis_plane_change_delta_v(v_apo::Float64, delta_i::Float64)
    return 2.0 * v_apo * sin(delta_i / 2.0)
end

"""
组合转移+平面变化 (Hohmann + 平面变化在远心点)
"""
function hohmann_with_plane_change(r1::Float64, r2::Float64, delta_i::Float64;
                                    mu::Float64=G_SUN)
    hoh = hohmann_transfer(r1, r2; mu=mu)

    a_trans = hoh.a_trans
    v_apo = sqrt(mu * (2.0/r2 - 1.0/a_trans))
    delta_v_plane = 2.0 * v_apo * sin(delta_i / 2.0)

    # 远心点合速度变化: √(Δv² + Δv_plane²) (若同时执行)
    # 或分开执行: Δv + Δv_plane
    delta_v_combined = sqrt(hoh.delta_v2^2 + delta_v_plane^2)
    delta_v_total = hoh.delta_v1 + delta_v_combined

    return (delta_v1=hoh.delta_v1, delta_v_plane=delta_v_plane,
            delta_v_total=delta_v_total, delta_v2_original=hoh.delta_v2)
end

# ============================================================
# 引力辅助 (Gravity Assist / Swing-by)
# ============================================================

"""
引力辅助: 在行星参考系中的速度变化

V_∞_in: 进入速度 (相对行星)
V_∞_out: 离开速度 (相对行星)
|V_∞_in| = |V_∞_out| (能量在行星系守恒)

转向角: δ = 2 arcsin(1 / (1 + r_p V_∞²/μ))

参数:
- V_inf: 双曲线超速 [km/s]
- r_p: 近心距 [km]
- mu: 行星引力参数

返回: (delta_deg, V_out_vec)
"""
function gravity_assist_turn_angle(V_inf::Float64, r_p::Float64; mu::Float64=398600.5)
    e = 1.0 + r_p * V_inf^2 / mu
    delta = 2.0 * asin(1.0 / e)
    return delta
end

"""
引力辅助后的日心速度变化

在行星参考系:
  V_∞_in = V_sc_in - V_planet
  V_∞_out 通过转向角 δ 计算出
  V_sc_out = V_∞_out + V_planet

日心速度增益最大:
  ΔV_max = 2 V_planet sin(δ/2)
"""
function gravity_assist_delta_v(V_planet::Float64, V_inf::Float64, r_p::Float64;
                                 mu::Float64=398600.5)
    delta = gravity_assist_turn_angle(V_inf, r_p; mu=mu)
    return 2.0 * V_planet * sin(delta / 2.0)
end

"""
引力辅助后的最大日心速度变化 (最佳对齐)

当天体速度矢量和飞越后的超速矢量同向 → 最大增益
"""
function max_gravity_assist_gain(V_planet::Float64, V_inf::Float64)
    return 2.0 * V_planet * V_inf / V_planet  # 近似
end

# ============================================================
# Patched Conics (简化)
# ============================================================

"""
Patched Conics 近似: 计算行星际转移的总 Δv

方案:
  1. 逃逸出发行星: v_∞_depart
  2. Hohmann 日心转移: Δv_hohmann
  3. 捕获到目标行星: v_∞_arrive

参数:
- r_depart: 出发行星轨道半径
- r_arrive: 目标行星轨道半径
- mu_planet_depart: 出发行星 μ
- mu_planet_arrive: 目标行星 μ
- r_p_depart: 出发行星近心距
- r_p_arrive: 目标行星近心距
- mu_sun: 太阳引力参数

返回: 总 Δv 分解
"""
function patched_conics_transfer(r_depart::Float64, r_arrive::Float64,
                                  mu_planet_depart::Float64, mu_planet_arrive::Float64,
                                  r_p_depart::Float64, r_p_arrive::Float64;
                                  mu_sun::Float64=G_SUN)
    # Step 1: 逃逸出发行星
    v_depart_planet = sqrt(mu_sun / r_depart)
    v_inf_depart = sqrt(mu_sun / r_depart) * abs(sqrt(2.0 * r_arrive / (r_depart + r_arrive)) - 1.0)

    # 逃逸所需超速
    v_esc_depart = sqrt(2.0 * mu_planet_depart / r_p_depart + v_inf_depart^2)
    delta_v_depart = v_esc_depart - sqrt(mu_planet_depart / r_p_depart)

    # Step 2: Hohmann 日心转移
    hoh = hohmann_transfer(r_depart, r_arrive; mu=mu_sun)
    delta_v_hohmann = hoh.delta_v_total

    # Step 3: 捕获到目标行星
    v_arrive_planet = sqrt(mu_sun / r_arrive)
    v_inf_arrive = sqrt(mu_sun / r_arrive) * abs(1.0 - sqrt(2.0 * r_depart / (r_depart + r_arrise)))
    # Note: fixed variable name typo r_arrise → r_arrive
    v_inf_arrive = sqrt(mu_sun / r_arrive) * abs(1.0 - sqrt(2.0 * r_depart / (r_depart + r_arrive)))

    v_capture = sqrt(v_inf_arrive^2 + 2.0 * mu_planet_arrive / r_p_arrive)
    delta_v_arrive = v_capture - sqrt(mu_planet_arrive / r_p_arrive)

    delta_v_total = delta_v_depart + delta_v_hohmann + delta_v_arrive

    return (delta_v_depart=delta_v_depart, delta_v_hohmann=delta_v_hohmann,
            delta_v_arrive=delta_v_arrive, delta_v_total=delta_v_total,
            transfer_time=hoh.transfer_time)
end

# ============================================================
# 发射窗口 (简化)
# ============================================================

"""
相位角 (Phasing angle) 用于 Hohmann 转移

出发时目标行星需要领先的角度:
  φ = π - n_target · T_transfer

参数:
- r_depart: 出发轨道半径
- r_arrive: 目标轨道半径
- mu: 引力参数

返回: φ [rad]
"""
function phasing_angle_hohmann(r_depart::Float64, r_arrive::Float64; mu::Float64=G_SUN)
    a_trans = (r_depart + r_arrive) / 2.0
    T_transfer = π * sqrt(a_trans^3 / mu)
    n_target = sqrt(mu / r_arrive^3)
    phi = π - n_target * T_transfer
    return mod(phi + 2π, 2π)
end

"""
会合周期 (Synodic period)

两次相同相对位置的时间间隔
  1/T_syn = |1/T1 - 1/T2|
"""
function synodic_period(T1::Float64, T2::Float64)
    return 1.0 / abs(1.0/T1 - 1.0/T2)
end

synodic_period_from_radii(r1::Float64, r2::Float64; mu::Float64=G_SUN) =
    synodic_period(2π*sqrt(r1^3/mu), 2π*sqrt(r2^3/mu))

export hohmann_transfer, hohmann_transfer_inward, bi_elliptic_transfer
export simple_plane_change_delta_v, apoapsis_plane_change_delta_v
export hohmann_with_plane_change
export gravity_assist_turn_angle, gravity_assist_delta_v, max_gravity_assist_gain
export patched_conics_transfer, phasing_angle_hohmann, synodic_period, synodic_period_from_radii
