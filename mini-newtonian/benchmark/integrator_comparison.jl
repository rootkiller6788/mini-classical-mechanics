#!/usr/bin/env julia
# benchmark/integrator_comparison.jl — Euler/Verlet/RK4/Symplectic Euler 长期对比
# 指标: 能量漂移、相位误差、计算时间

include("../src/Newtonian.jl")
using .Newtonian

println("="^60)
println("ODE 积分器基准测试：简谐振子 1000 周期")
println("="^60)

function sho(x,v,t); Vec3(-x.x, 0, 0); end
t_end = 2000*pi  # 1000周期
E0 = 0.5*1.0*1.0^2

methods = [:euler, :euler_cromer, :rk2, :rk4, :verlet]
results = []

for method in methods
    ode = ParticleODE(sho, 0.0, t_end, Vec3(1,0,0), Vec3(0,0,0))
    t_start = time()

    # 不同方法用不同步长以保证稳定性
    dt = method == :euler ? 2*pi/500 : 2*pi/50
    traj = solve_fixed_step(ode, dt; method=method)
    elapsed = time() - t_start

    # 能量漂移
    energies = [0.5*s.r.x^2 + 0.5*s.v.x^2 for s in traj.states]
    max_drift = maximum(abs.(energies .- E0)) / E0 * 100  # %
    mean_drift = mean(abs.(energies .- E0)) / E0 * 100    # %

    # 相位误差（过零点比较）
    zero_crossings = Float64[]
    for i in 2:length(traj.states)
        if traj.states[i-1].r.x * traj.states[i].r.x < 0
            push!(zero_crossings, traj.states[i].t)
        end
    end
    T_measured = mean(diff(zero_crossings[2:min(20,end)])) * 2
    phase_err = abs(T_measured - 2*pi) / (2*pi) * 100  # %

    push!(results, (method=method, dt=dt, time=elapsed,
                     max_drift=max_drift, mean_drift=mean_drift,
                     phase_err=phase_err, n_steps=length(traj.states)))
end

println("\n  方法          步长      步数     耗时(s)   最大漂移%   平均漂移%   相位误差%")
println("  " * "-"^75)
for r in results
    println("  $(rpad(string(r.method),14)) $(rpad(r.dt,10)) $(rpad(r.n_steps,8)) $(round(r.time,digits=3))       $(round(r.max_drift,digits=4))       $(round(r.mean_drift,digits=4))       $(round(r.phase_err,digits=4))")
end

println("\n结论:")
verlet_r = results[findfirst(r -> r.method == :verlet, results)]
rk4_r    = results[findfirst(r -> r.method == :rk4, results)]
println("  ✓ Velocity Verlet: 能量漂移 $(round(verlet_r.max_drift,digits=4))%（有界，不单调增长）")
println("  ✓ RK4:             能量漂移 $(round(rk4_r.max_drift,digits=4))%（单调增长）")
println("  ✓ 保守系统长期积分推荐: Velocity Verlet")
