#!/usr/bin/env julia
# examples/central_force.jl — 中心力问题 / 轨道力学
# 开普勒问题、有效势、轨道类型

include("../src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  Central Force Motion — 中心力与轨道力学")
    println("="^60)

    # ---- 1. 圆形轨道 ----
    println("\n--- 1. Circular orbit (Earth around Sun) ---")
    M_sun = M_SUN
    m_earth = M_EARTH
    r_orbit = AU  # 1 AU

    # 圆形轨道速度: v = sqrt(G*M/r)
    v_circ = sqrt(G * M_sun / r_orbit)
    T_circ = 2π * r_orbit / v_circ  # 轨道周期
    println("Circular orbit velocity: $(round(v_circ/1000, digits=2)) km/s")
    println("Orbital period: $(round(T_circ/DAY, digits=1)) days (actual: 365.25)")

    # 数值模拟
    r0 = Vec3(r_orbit, 0.0, 0.0)
    v0 = Vec3(0.0, v_circ, 0.0)

    function sun_gravity(r, v, t)
        return gravity_acceleration(r, Vec3(), M_sun)
    end

    dt = DAY  # 1 天步长
    t_end = 2 * YEAR
    traj = solve_fixed_step(sun_gravity, r0, v0, t_end, dt; method=:verlet, record_every=10)

    n = length(traj.ts)
    println("Simulated $(n) steps, final t = $(round(traj.ts[end]/DAY, digits=1)) days")

    # 检查轨道是否保持圆形
    distances = [norm(p) for p in traj.positions]
    d_min, d_max = minimum(distances), maximum(distances)
    drift_pct = (d_max - d_min) / r_orbit * 100
    println("Orbit radius: min=$(round(d_min/AU,digits=4)) AU, max=$(round(d_max/AU,digits=4)) AU")
    println("Radius variation: $(round(drift_pct, digits=4))% (should be << 1% for Verlet)")
    println("Expected radius: 1.0 AU")

    # ---- 2. 椭圆轨道 (e = 0.5) ----
    println("\n--- 2. Elliptical orbit (e = 0.5) ---")
    a = r_orbit  # 半长轴
    e = 0.5      # 离心率
    # 近地点速度: v_peri = sqrt(G*M*(1+e)/(a*(1-e)))
    r_peri = a * (1 - e)
    v_peri = sqrt(G * M_sun * (1 + e) / (a * (1 - e)))

    r0_ellip = Vec3(r_peri, 0.0, 0.0)
    v0_ellip = Vec3(0.0, v_peri, 0.0)

    traj_ellip = solve_fixed_step(sun_gravity, r0_ellip, v0_ellip, 2*YEAR, dt; method=:verlet, record_every=10)

    distances_ellip = [norm(p) for p in traj_ellip.positions]
    r_min_sim = minimum(distances_ellip)
    r_max_sim = maximum(distances_ellip)
    e_sim = (r_max_sim - r_min_sim) / (r_max_sim + r_min_sim)
    println("Simulated eccentricity: $(round(e_sim, digits=4)) (target: $e)")
    println("Perihelion: $(round(r_min_sim/AU,digits=4)) AU, Aphelion: $(round(r_max_sim/AU,digits=4)) AU")

    # 开普勒第三定律验证: T² ∝ a³
    # 找轨道周期（两次经过 x 轴正半轴的时间差）
    cross_times = Float64[]
    for i in 2:length(traj_ellip.positions)
        p_prev = traj_ellip.positions[i-1]
        p_curr = traj_ellip.positions[i]
        # 穿越 x 轴正半轴：y 符号改变且 x > 0
        if p_prev.y * p_curr.y < 0 && p_curr.x > 0
            t_cross = traj_ellip.ts[i-1] + (traj_ellip.ts[i] - traj_ellip.ts[i-1]) * abs(p_prev.y) / abs(p_curr.y - p_prev.y)
            push!(cross_times, t_cross)
            if length(cross_times) >= 3
                break
            end
        end
    end
    if length(cross_times) >= 3
        T_sim = cross_times[3] - cross_times[1]
        T_kepler = 2π * sqrt(a^3 / (G * M_sun))
        println("Orbital period (sim): $(round(T_sim/DAY,digits=1)) days")
        println("Orbital period (Kepler): $(round(T_kepler/DAY,digits=1)) days")
        println("Kepler 3rd law error: $(round(abs(T_sim - T_kepler)/T_kepler*100, digits=4))%")
    end

    # ---- 3. 逃逸轨道 (v > v_esc) ----
    println("\n--- 3. Escape trajectory ---")
    v_esc = sqrt(2 * G * M_sun / r_peri)  # 逃逸速度
    v_hyper = v_esc * 1.2  # 双曲线轨道

    r0_esc = Vec3(r_peri, 0.0, 0.0)
    v0_esc = Vec3(0.0, v_hyper, 0.0)

    traj_esc = solve_fixed_step(sun_gravity, r0_esc, v0_esc, 5*YEAR, dt*5; method=:verlet, record_every=50)
    r_final = norm(traj_esc.positions[end])
    println("Initial distance: $(round(r_peri/AU,digits=4)) AU")
    println("Final distance (after 5 yr): $(round(r_final/AU,digits=2)) AU — escaping!")

    # ---- 4. 有效势分析 ----
    println("\n--- 4. Effective potential ---")
    L_mag = norm(angular_momentum(r0, v0, m_earth))
    println("Angular momentum: $(round(L_mag, digits=2)) kg·m²/s")

    # U_eff(r) = -GMm/r + L²/(2mr²)
    function U_eff(r, M, m, L_val)
        return -G * M * m / r + L_val^2 / (2 * m * r^2)
    end

    rs = [0.5AU, 0.8AU, AU, 1.5AU, 2.0AU]
    println("r [AU] | U_eff [J]")
    println("-"^30)
    for r_val in rs
        u = U_eff(r_val, M_sun, m_earth, L_mag)
        println("$(round(r_val/AU,digits=1))     | $(round(u, digits=-20))")
    end

    println("\n✅ central_force.jl done.")
end

main()
