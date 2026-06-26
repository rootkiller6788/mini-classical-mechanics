#!/usr/bin/env julia
# examples/oscillator.jl — 谐振子：简谐、阻尼、受迫、共振

include("../src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Harmonic Oscillator — 谐振子系统")
    println("="^60)

    # ---- 1. 简谐振子 ----
    println("\n--- 1. Simple Harmonic Oscillator ---")
    k = 10.0    # 弹性系数
    m = 1.0     # 质量
    omega0 = sqrt(k / m)
    T0 = 2π / omega0
    x0 = 1.0    # 初始位移

    println("ω₀ = $(round(omega0, digits=3)) rad/s, T₀ = $(round(T0, digits=3)) s")

    function sho_accel(r, v, t)
        return Vec3(-omega0^2 * r.x, 0.0, 0.0)
    end

    r0 = Vec3(x0, 0.0, 0.0)
    v0 = Vec3(0.0, 0.0, 0.0)

    traj_sho = solve_fixed_step(sho_accel, r0, v0, 5*T0, 0.001; method=:verlet, record_every=10)

    # 解析解验证: x(t) = x0*cos(ω₀*t)
    x_analytic = x0 * cos(omega0 * traj_sho.ts[end])
    x_numeric  = traj_sho.positions[end].x
    println("t = $(round(traj_sho.ts[end], digits=3)): analytic x=$(round(x_analytic,digits=6)), numeric x=$(round(x_numeric,digits=6))")
    println("Error: $(abs(x_analytic - x_numeric))")

    # 能量守恒
    function sho_potential(r) ; return 0.5 * k * r.x^2 ; end
    energies, drift = energy_drift(traj_sho, m, sho_potential)
    println("Energy drift (max relative): $(round(maximum(abs.(drift)), digits=10))")

    # ---- 2. 阻尼谐振子（欠阻尼、临界阻尼、过阻尼） ----
    println("\n--- 2. Damped Harmonic Oscillator ---")
    b_under  = 0.5       # 欠阻尼: b < 2√(km)
    b_crit   = 2 * sqrt(k * m)  # 临界阻尼
    b_over   = 10.0      # 过阻尼

    for (label, b_val) in [("Underdamped", b_under), ("Critical", b_crit), ("Overdamped", b_over)]
        function damped_accel(r, v, t)
            return Vec3(-omega0^2 * r.x - (b_val/m) * v.x, 0.0, 0.0)
        end
        traj = solve_fixed_step(damped_accel, r0, v0, 5*T0, 0.001; method=:verlet, record_every=50)
        x_final = traj.positions[end].x
        println("  $label (b=$b_val): x(t=5T₀) = $(round(x_final, digits=6))")
    end

    # ---- 3. 受迫谐振子 + 共振 ----
    println("\n--- 3. Forced Oscillator & Resonance ---")
    F0 = 1.0  # 驱动力幅值
    b_res = 0.2  # 小阻尼

    function forced_accel(r, v, t, omega_drive)
        driving = F0 * cos(omega_drive * t) / m
        return Vec3(driving - omega0^2 * r.x - (b_res/m) * v.x, 0.0, 0.0)
    end

    r0_forced = Vec3(0.0, 0.0, 0.0)
    v0_forced = Vec3(0.0, 0.0, 0.0)

    # 扫频：测量稳态振幅
    omegas = [0.5, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.5, 2.0] .* omega0
    println("  ω/ω₀ | Amplitude")
    println("  " * "-"^25)

    for omega_d in omegas
        accel_func(r, v, t) = forced_accel(r, v, t, omega_d)
        traj_f = solve_fixed_step(accel_func, r0_forced, v0_forced, 20*T0, 0.005;
                                  method=:verlet, record_every=50)

        # 稳态振幅（取后一半的平均峰峰值）
        n_half = length(traj_f.ts) ÷ 2
        xs = [traj_f.positions[i].x for i in n_half:length(traj_f.ts)]
        amplitude = (maximum(xs) - minimum(xs)) / 2
        println("  $(round(omega_d/omega0, digits=2))   | $(round(amplitude, digits=4))")
    end

    # 理论共振频率: ω_res = √(ω₀² - b²/(2m²))
    omega_res_theory = sqrt(omega0^2 - b_res^2 / (2 * m^2))
    println("\nTheoretical resonance: ω_res/ω₀ = $(round(omega_res_theory/omega0, digits=4))")

    # ---- 4. 耦合谐振子（两质量三弹簧） ----
    println("\n--- 4. Coupled Oscillators ---")
    k_couple = 2.0  # 耦合弹簧系数
    m1, m2 = 1.0, 1.0
    k1, k2 = k, k  # 边弹簧

    function coupled_accel(r, v, t)
        x1, x2 = r.x, r.y
        a1 = (-(k1 + k_couple) * x1 + k_couple * x2) / m1
        a2 = (k_couple * x1 - (k2 + k_couple) * x2) / m2
        return Vec3(a1, a2, 0.0)
    end

    r0_coupled = Vec3(1.0, -0.5, 0.0)
    v0_coupled = Vec3(0.0, 0.0, 0.0)

    traj_coupled = solve_fixed_step(coupled_accel, r0_coupled, v0_coupled, 5*T0, 0.001;
                                    method=:verlet, record_every=10)
    println("Coupled oscillators simulated for $(length(traj_coupled.ts)) steps")

    # 简正模频率（理论）
    omega_sym = sqrt(k / m)  # 对称模: x1 = x2
    omega_anti = sqrt((k + 2*k_couple) / m)  # 反对称模: x1 = -x2
    println("Normal modes: ω_sym/ω₀ = $(round(omega_sym/omega0,digits=3)), ω_anti/ω₀ = $(round(omega_anti/omega0,digits=3))")

    println("\n✅ oscillator.jl done.")
end

main()
