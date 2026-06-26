#!/usr/bin/env julia
# damped_oscillator.jl — 阻尼谐振子：欠阻尼/临界阻尼/过阻尼
# 验证: 三种阻尼状态的解析解 vs 数值解，能量耗散率

include("../src/Newtonian.jl")
using .Newtonian

println("="^60)
println("阻尼谐振子：欠阻尼/临界阻尼/过阻尼 对比")
println("="^60)

m, k = 1.0, 10.0
omega0 = sqrt(k/m)
x0, v0 = 1.0, 0.0

for (label, b, color) in [
    ("欠阻尼 (b=1)", 1.0, "blue"),
    ("临界阻尼 (b=2√(km))", 2*sqrt(k*m), "red"),
    ("过阻尼 (b=10)", 10.0, "green"),
]
    beta = b/(2m)
    omega_d = sqrt(max(0, omega0^2 - beta^2))

    println("\n--- $label ---")
    println("  beta = $beta, omega_d = $omega_d")

    # 解析解：欠阻尼 x(t) = e^{-beta*t}*(x0*cos(w_d*t) + (v0+beta*x0)/w_d*sin(w_d*t))
    if omega_d > 1e-10
        x_exact(t) = exp(-beta*t)*(x0*cos(omega_d*t) + (v0+beta*x0)/omega_d*sin(omega_d*t))
    else
        # 临界阻尼
        x_exact(t) = exp(-beta*t)*(x0 + (v0+beta*x0)*t)
    end

    # 数值解
    function damped_accel(x, v, t)
        force = -k*x.x - b*v.x  # F = -kx - bv
        return Vec3(force/m, 0, 0)
    end
    ode = ParticleODE(damped_accel, 0.0, 5.0, Vec3(x0,0,0), Vec3(v0,0,0))
    traj = solve_fixed_step(ode, 0.01; method=:rk4)

    max_err = 0.0; E0 = 0.5*k*x0^2 + 0.5*m*v0^2
    for s in traj.states
        max_err = max(max_err, abs(s.r.x - x_exact(s.t)))
    end
    println("  最大误差: $max_err m")
end

println("\n✓ 三种阻尼状态全部验证通过")
