#!/usr/bin/env julia
# examples/non_inertial.jl — 非惯性参考系
# 傅科摆、科里奥利力、离心力、旋转系中的运动

include("../src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Non-Inertial Reference Frames")
    println("="^60)

    g_val = G_EARTH
    omega_earth = 2π / DAY  # 地球自转角速度 [rad/s]

    # ---- 1. 科里奥利力：落体偏东 ----
    println("\n--- 1. Coriolis Deflection — Eastward Drift of Falling Body ---")

    # 从高度 H 自由下落，在北纬 latitude 处
    H     = 100.0     # 下落高度 [m]
    lat   = deg2rad(42)  # 纬度（波士顿 ≈ 42°N）
    omega = omega_earth

    # 地球自转轴（近似静止系，忽略地球公转）
    # ω 矢量指向北极星方向
    omega_vec = Vec3(0.0, omega * cos(lat), omega * sin(lat))  # x=东, y=北, z=上

    r0 = Vec3(0.0, 0.0, H)
    v0 = Vec3(0.0, 0.0, 0.0)

    function falling_body_accel(r, v, t)
        # 重力 + 科里奥利力（离心力很小，忽略）
        gravity = Vec3(0.0, 0.0, -g_val)
        coriolis = -2.0 * cross(omega_vec, v)
        return gravity + coriolis
    end

    t_flight = sqrt(2 * H / g_val)
    dt = t_flight / 1000
    traj = solve_fixed_step(falling_body_accel, r0, v0, t_flight * 1.1, dt; method=:verlet, record_every=10)

    r_final = traj.positions[end]
    deflection_x = r_final.x * 100  # cm
    deflection_y = r_final.y * 100  # cm

    println("Fall height: $H m, latitude: $(round(rad2deg(lat))) deg N")
    println("Flight time: $(round(t_flight, digits=2)) s")
    println("Eastward deflection: $(round(deflection_x, digits=3)) cm")
    println("Southward deflection: $(round(deflection_y, digits=4)) cm (should be ~0)")

    # 理论值：Δx ≈ (2/3)*ω*cos(λ)*t_flight³*g 或等效公式
    # Δx ≈ (1/3)*ω*cos(λ)*g*t_flight³ 更准确
    theory_deflection = (1/3) * omega * cos(lat) * g_val * t_flight^3 * 100  # cm
    println("Theoretical eastward deflection: $(round(theory_deflection, digits=2)) cm")

    # ---- 2. 傅科摆 ----
    println("\n--- 2. Foucault Pendulum ---")
    L_foucault = 67.0        # 摆长 [m]（万神殿实际值）
    m_pendulum = 28.0        # 质量 [kg]（实际摆锤）
    theta0_f = deg2rad(10.0) # 初始偏角

    # 小角度近似下，摆动平面的旋转角速度 = ω*sin(λ)
    precession_rate = omega * sin(lat)  # rad/s
    precession_period_hours = 2π / precession_rate / 3600
    println("Pendulum length: $L_foucault m")
    println("Precession rate: $(round(rad2deg(precession_rate*3600), digits=4)) deg/hour")
    println("Precession period: $(round(precession_period_hours, digits=1)) hours")
    println("(Theoretical at pole: 24h, at Paris 48.9°N: ~32h)")

    function foucault_pendulum_accel(r, v, t)
        # 摆的位置约束在球面上：|r| = L
        # 小角度近似：r = (x, y, -L*cosθ) ≈ (x, y, -L)
        # 张力 ≈ mg（小角度）
        tension_magnitude = m_pendulum * g_val
        # 恢复力：-T*(x/L, y/L, 0) ≈ -mg/L * (x, y, 0)
        spring_constant = g_val / L_foucault

        restoring = Vec3(-spring_constant * r.x, -spring_constant * r.y, 0.0)
        coriolis = -2.0 * cross(omega_vec, v)

        return restoring + coriolis
    end

    # 在 xy 平面释放（初始在 x 方向）
    r0_f = Vec3(L_foucault * sin(theta0_f), 0.0, L_foucault * (1 - cos(theta0_f)))
    v0_f = Vec3(0.0, 0.0, 0.0)

    # 模拟 24 小时的一半看进动
    sim_time = 12 * 3600  # 12 hours
    dt_f = 0.1            # 0.1 s step

    traj_foucault = solve_fixed_step(foucault_pendulum_accel, r0_f, v0_f, sim_time, dt_f;
                                     method=:verlet, record_every=1000)

    n_recorded = length(traj_foucault.ts)

    # 分析摆动平面的偏转
    # 记录每次经过 x 轴正半轴的角度
    crossing_angles = Float64[]
    for i in 2:n_recorded
        p_prev = traj_foucault.positions[i-1]
        p_curr = traj_foucault.positions[i]
        # x>0 且 y 符号改变 → 经过 x 正半轴
        if p_curr.x > 0 && p_prev.y * p_curr.y < 0
            angle_val = atan(p_curr.y, p_curr.x)
            push!(crossing_angles, angle_val)
        end
    end

    if length(crossing_angles) >= 3
        # 拟合进动率
        angle_shift_per_hour = (crossing_angles[end] - crossing_angles[1]) / (sim_time / 3600)
        angle_shift_deg_per_hour = rad2deg(angle_shift_per_hour)
        println("\nObserved precession (from simulation):")
        println("  Angle shift: $(round(angle_shift_deg_per_hour, digits=3)) deg/hour")
        println("  Equivalent period: $(round(360/abs(angle_shift_deg_per_hour), digits=1)) hours")
        theory_deg_per_hour = rad2deg(precession_rate * 3600)
        println("  Theoretical: $(round(theory_deg_per_hour, digits=3)) deg/hour")
    end

    # ---- 3. 离心力：旋转太空栖息地 ----
    println("\n--- 3. Artificial Gravity (Rotating Space Habitat) ---")
    R_habitat = 100.0   # 半径 [m]
    g_target  = 9.8     # 目标人工重力

    # 需要的角速度：ω²R = g → ω = √(g/R)
    omega_habitat = sqrt(g_target / R_habitat)
    rpm = omega_habitat * 60 / (2π)
    println("Radius: $R_habitat m, Target gravity: $g_target m/s²")
    println("Required rotation: $(round(omega_habitat, digits=3)) rad/s = $(round(rpm, digits=2)) RPM")

    # 在旋转系中模拟一个"下落"的物体
    omega_hab_vec = Vec3(0.0, 0.0, omega_habitat)

    function habitat_accel(r, v, t)
        # 重力 → 离心力 → 科里奥利力
        centrifugal = -cross(omega_hab_vec, cross(omega_hab_vec, r))
        coriolis    = -2.0 * cross(omega_hab_vec, v)
        return centrifugal + coriolis
    end

    # 从"地面"释放一个物体
    r_ground = Vec3(R_habitat, 0.0, 0.0)  # "地面"在圆柱内壁
    r_drop   = Vec3(R_habitat - 1.0, 0.0, 0.0)  # 离地 1m
    v_drop   = Vec3(0.0, 0.0, 0.0)

    traj_habitat = solve_fixed_step(habitat_accel, r_drop, v_drop, 5.0, 0.01;
                                    method=:verlet, record_every=10)

    r_final_h = traj_habitat.positions[end]
    println("\nDrop test from 1m above 'ground' in rotating habitat:")
    println("  Initial position: $(r_drop)")
    println("  After 5s: $(r_final_h)")
    println("  Note: Coriolis deflects falling objects — apparent 'gravity' is not perfectly radial")

    # ---- 4. 地转偏向的日常影响 ----
    println("\n--- 4. Coriolis in Daily Life ---")

    # 浴缸排水涡旋（北半球逆时针）
    # 科里奥利参数 f = 2ω*sin(λ)
    f_param = 2 * omega * sin(lat)
    println("Coriolis parameter f at $(round(rad2deg(lat)))°N: $(round(f_param, digits=8)) s⁻¹")
    println("Rossby number for bathtub (U=0.1m/s, L=0.5m): Ro = U/(f*L) = $(round(0.1/(f_param*0.5), digits=1))")
    println("→ Ro >> 1: Coriolis negligible in bathtubs (drain shape dominates)")
    println("→ Ro ~ 1: Coriolis significant for weather systems (U~10m/s, L~1000km)")

    println("\n✅ non_inertial.jl done.")
end

main()
