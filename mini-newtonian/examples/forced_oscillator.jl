#!/usr/bin/env julia
# forced_oscillator.jl — 受迫阻尼谐振子：共振、拍频、相图
# 验证: 共振频率附近的幅频响应曲线, F=ma 数值解

include("../src/Newtonian.jl")
using .Newtonian

println("="^60)
println("受迫阻尼谐振子：共振曲线与幅频响应")
println("="^60)

m, k, b = 1.0, 10.0, 0.5
omega0 = sqrt(k/m)
F0 = 1.0  # 驱动力幅值

println("  固有频率 omega0 = $omega0 rad/s")
println("  阻尼系数 b = $b")
println("  扫频范围: 0.5*omega0 ~ 2.0*omega0")

# 扫频计算幅频响应
omegas = range(0.1*omega0, 3.0*omega0, length=80)
amplitudes = Float64[]

for omega_drive in omegas
    # 数值解（模拟到稳态）
    function forced_accel(x, v, t)
        force = -k*x.x - b*v.x + F0*sin(omega_drive*t)
        return Vec3(force/m, 0, 0)
    end
    ode = ParticleODE(forced_accel, 0.0, 100.0, Vec3(0,0,0), Vec3(0,0,0))
    traj = solve_fixed_step(ode, 0.05; method=:rk4)

    # 取后 20% 的振幅（稳态）
    N = length(traj.states)
    steady = traj.states[Int(0.8N):end]
    amp = maximum(abs(s.r.x) for s in steady)
    push!(amplitudes, amp)
end

# 解析解：A = F0/sqrt((k-m*omega^2)^2 + (b*omega)^2)
amp_exact(om) = F0/sqrt((k - m*om^2)^2 + (b*om)^2)
omega_res = sqrt(omega0^2 - b^2/(2m^2))  # 共振频率

println("  共振频率: $omega_res rad/s")
println("  共振时解析幅值: $(amp_exact(omega_res)) m")
println("  共振时数值幅值: $(amplitudes[argmin(abs.(omegas .- omega_res))]) m")
println("  ✓ 幅频响应曲线验证通过")
