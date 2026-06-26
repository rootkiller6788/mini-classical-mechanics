#!/usr/bin/env julia
# gyroscope.jl — 陀螺仪原理
# 参考: Goldstein Ch.5.8, Landau Ch.36
# 演示: 科里奥利力矩、进动率、单/双自由度陀螺仪

include("../src/RigidBody.jl"); using .RigidBody

function main()
    println("="^60)
    println("  Gyroscope Principles")
    println("="^60)

    # ==============================================
    # 陀螺仪基本参数
    # ==============================================
    # 高速旋转的圆盘
    M_rotor = 0.2      # 转子质量
    R_rotor = 0.05     # 转子半径
    I_spin = 0.5 * M_rotor * R_rotor^2  # 转子的自转转动惯量
    ω_spin = 200.0     # 自转角速度 (~1900 RPM)

    L_spin = I_spin * ω_spin
    println("\n转子参数:")
    println("  质量 M = $M_rotor kg, 半径 R = $R_rotor m")
    println("  自转惯量 I₃ = $(round(I_spin, 8)) kg·m²")
    println("  自转角速度 ω = $ω_spin rad/s ($(round(ω_spin*60/(2π))) RPM)")
    println("  自转角动量 L = $(round(L_spin, 6)) kg·m²/s")

    # ==============================================
    # 实验 1: 陀螺仪进动
    # ==============================================
    println("\n" * "-"^40)
    println("实验 1: 陀螺仪进动 (稳态)")

    # 在外力矩 N = mgℓ 下，进动角速度 Ω = N/L
    M_gimbal = 0.05   # 万向架质量
    l_arm = 0.1       # 力臂长度
    g = 9.81

    N_ext = M_gimbal * g * l_arm  # 重力矩
    Ω_prec = gyroscope_precession(N_ext, L_spin)
    T_prec = 2π / Ω_prec          # 进动周期

    println("  外力矩 N = $(round(N_ext, 6)) N·m")
    println("  进动角速率 Ω = $(round(Ω_prec, 4)) rad/s")
    println("  进动周期 T = $(round(T_prec, 2)) s")

    # ==============================================
    # 实验 2: 不同角动量下的进动率
    # ==============================================
    println("\n" * "-"^40)
    println("实验 2: 进动率 vs 自转角动量")

    println("  ω_spin (RPM)  |  L_spin      |  Ω (rad/s)  |  T_prec (s)")
    println("  " * "-"^55)
    for ω_test in [50, 100, 200, 400, 800]
        L_test = I_spin * ω_test
        Ω_test = gyroscope_precession(N_ext, L_test)
        T_test = 2π / max(Ω_test, 1e-15)
        rpm = round(ω_test * 60 / (2π))
        println("  $(rpad(rpm,13)) | $(rpad(round(L_test,6),12)) | $(rpad(round(Ω_test,6),11)) | $(round(T_test,2))")
    end

    # ==============================================
    # 实验 3: 科里奥利力矩 (强制旋转)
    # ==============================================
    println("\n" * "-"^40)
    println("实验 3: 科里奥利力矩 (陀螺仪感受的力矩)")

    # 陀螺仪被强制绕空间轴旋转
    Ω_forced = [0.0, 0.0, 0.5]  # 绕 z 轴旋转 0.5 rad/s
    L_vec = [0.0, 1.0, 0.0]     # 角动量沿 y 轴

    N_gyro = gyroscopic_torque(L_vec, Ω_forced)
    println("  L = $(L_vec), Ω_forced = $(Ω_forced)")
    println("  → N_gyro = $(round.(N_gyro, digits=6))")

    # 验证: N = L × Ω
    N_expected = cross(L_vec, Ω_forced)
    println("  期望: L × Ω = $(round.(N_expected, digits=6))")
    println("  匹配: $(maximum(abs.(N_gyro - N_expected)) < 1e-14)")

    # 各种配置
    println("\n  不同配置:")
    for (L_vec_cfg, Ω_cfg, desc) in [
        ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], "L∥x, Ω∥y → N∥z"),
        ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], "L∥z, Ω∥x → N∥y"),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 1.0], "L∥x, Ω∥z → N∥-y"),
    ]
        N_cfg = gyroscopic_torque(L_vec_cfg, Ω_cfg)
        println("    $desc: N = $(round.(N_cfg, digits=3))")
    end

    # ==============================================
    # 实验 4: 外力矩驱动下的刚体转动
    # ==============================================
    println("\n" * "-"^40)
    println("实验 4: 恒力矩加速刚体")

    I_cfg = InertiaTensor(1.0, 2.0, 0.5, 0.0, 0.0, 0.0)
    w0_still = [0.0, 0.0, 0.0]

    # 绕 z 轴的恒力矩 (对角化 I 中 Izz=0.5)
    constant_torque_z(t, w) = [0.0, 0.0, 0.5]

    times_acc, omegas_acc = simulate_rigid_body_with_torque(
        I_cfg, w0_still, constant_torque_z, 2.0, 0.05)

    println("  力矩 N = [0, 0, 0.5] N·m")
    println("  Izz = $(I_cfg.Izz) kg·m²")
    println("  理论: ωz(t) = Nz/Izz · t = $(0.5/0.5)·t = t")
    println("  模拟终值 ωz = $(round(omegas_acc[end][3], 4)) rad/s")
    println("  期望 ωz = $(round(2.0, 4)) rad/s")
    println("  误差 = $(round(abs(omegas_acc[end][3] - 2.0), 6))")

    # ==============================================
    # 应用: 陀螺罗经 (Gyrocompass)
    # ==============================================
    println("\n" * "-"^40)
    println("应用: 陀螺罗经原理")

    # 地球自转在陀螺仪上产生的力矩
    ω_earth = 2π / 86400  # 地球自转角速度 (~7.27e-5 rad/s)
    L_gyrocompass = 10.0  # 陀螺仪角动量

    # 在纬度 λ 处，陀螺仪感受到的进动力矩指向北
    for lat in [0, 30, 60, 90]
        λ = lat * π / 180
        # 地球自转向量在水平和垂直方向的分量
        ω_h = ω_earth * cos(λ)  # 水平分量
        ω_v = ω_earth * sin(λ)  # 垂直分量
        N_compass = L_gyrocompass * ω_h
        println("  纬度 $(rpad(lat,2))°: ω_h=$(round(ω_h*1e5,2))e-5 rad/s, 寻北力矩 ≈ $(round(N_compass,8)) N·m")
    end

    println("\n✅ gyroscope.jl 完成.")
end

function cross(a, b)
    return [a[2]*b[3] - a[3]*b[2],
            a[3]*b[1] - a[1]*b[3],
            a[1]*b[2] - a[2]*b[1]]
end

main()
