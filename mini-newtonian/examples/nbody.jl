#!/usr/bin/env julia
# examples/nbody.jl — N 体引力模拟
# 太阳系简化（太阳 + 地球 + 火星）或随机 N 体

include("../src/Newtonian.jl")
using .Newtonian

function main()
    println("="^60)
    println("  N-Body Gravitational Simulation")
    println("="^60)

    # ---- 1. 太阳-地球-火星 三体简化 ----
    println("\n--- 1. Sun-Earth-Mars (restricted 3-body) ---")

    # 使用天文单位制: G=G_AU, M in Msun, r in AU, t in years
    G_solar = G_AU

    masses = [1.0, 3.0e-6, 3.2e-7]  # Msun: Sun, Earth, Mars
    labels = ["Sun", "Earth", "Mars"]

    # 初始位置和速度（近似圆形轨道，共面）
    pos = Vec3[
        Vec3(0.0, 0.0, 0.0),            # Sun at origin
        Vec3(1.0, 0.0, 0.0),            # Earth at 1 AU
        Vec3(0.0, 1.524, 0.0),           # Mars at 1.524 AU
    ]
    vel = Vec3[
        Vec3(0.0, 0.0, 0.0),            # Sun (approx stationary)
        Vec3(0.0, 2π, 0.0),             # Earth: v = 2π AU/yr
        Vec3(-2π / sqrt(1.524), 0.0, 0.0),  # Mars: v = 2π/√a AU/yr
    ]

    st = NBodyState(masses, pos, vel)

    # N体加速度函数
    function nbody_accel(st::NBodyState)
        n = n_particles(st)
        acc = [Vec3() for _ in 1:n]
        for i in 1:n
            for j in 1:n
                if i != j
                    acc[i] = acc[i] + gravity_acceleration(
                        st.positions[i], st.positions[j], st.masses[j]; G_val=G_solar
                    )
                end
            end
        end
        return acc
    end

    dt = 0.001  # yr (~0.365 days)
    t_end = 3.0 # 3 years
    n_steps = Int(ceil(t_end / dt))

    println("Simulating $(n_steps) steps...")

    # 记录轨迹
    records_every = 100
    trajectories = [Trajectory() for _ in 1:3]

    for step in 1:n_steps
        if step % records_every == 1
            for i in 1:3
                record!(trajectories[i], st.t, st.positions[i], st.velocities[i])
            end
        end
        st = nbody_verlet_step(nbody_accel, st, dt)
    end

    # 分析
    for i in 1:3
        traj = trajectories[i]
        dists = [norm(p) for p in traj.positions]
        println("$(labels[i]): $(length(traj.ts)) records, " *
                "distance range: [$(round(minimum(dists),digits=4)), $(round(maximum(dists),digits=4))] AU")
    end

    # ---- 2. 随机 N 体坍缩（玩具模型） ----
    println("\n--- 2. Random N-body collapse (toy model, N=20) ---")
    N_random = 20
    rng_seed = 42
    # 简单随机初始化（不用Random包，用手写LCG）
    function simple_rand(seed)
        seed = (1103515245 * seed + 12345) & 0x7fffffff
        return seed, seed / 0x7fffffff
    end

    seed = rng_seed
    masses_rand = Float64[]
    pos_rand = Vec3[]
    vel_rand = Vec3[]

    for i in 1:N_random
        seed, rx = simple_rand(seed)
        seed, ry = simple_rand(seed)
        seed, rz = simple_rand(seed)
        seed, vx = simple_rand(seed)
        seed, vy = simple_rand(seed)
        seed, vz = simple_rand(seed)
        push!(masses_rand, 1.0)  # 等质量
        push!(pos_rand, Vec3((rx-0.5)*10, (ry-0.5)*10, (rz-0.5)*2))
        push!(vel_rand, Vec3((vx-0.5)*0.5, (vy-0.5)*0.5, (vz-0.5)*0.5))
    end

    st_rand = NBodyState(masses_rand, pos_rand, vel_rand)
    dt_rand = 0.01
    t_end_rand = 10.0

    # 能量记录
    energies_history = Float64[]
    E0 = nothing

    for step in 1:Int(ceil(t_end_rand / dt_rand))
        if step % 10 == 1
            T = total_kinetic_energy(st_rand.masses, st_rand.velocities)
            U = total_gravitational_potential(st_rand.masses, st_rand.positions)
            E = T + U
            if E0 === nothing
                E0 = E
            end
            push!(energies_history, E)
        end
        st_rand = nbody_verlet_step(nbody_accel, st_rand, dt_rand)
    end

    E_rel = [(e - E0) / abs(E0) for e in energies_history]
    println("N=$(N_random): initial E=$(round(E0,digits=2)), final E=$(round(energies_history[end],digits=2))")
    println("Energy drift (max relative): $(round(maximum(abs.(E_rel)), digits=8))")
    println("Velocity Verlet energy conservation: $(maximum(abs.(E_rel)) < 1e-4 ? "PASS" : "WARN")")

    println("\n✅ nbody.jl done.")
end

main()
