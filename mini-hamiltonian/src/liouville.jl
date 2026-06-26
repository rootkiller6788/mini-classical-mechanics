# liouville.jl — Liouville theorem & phase space volume conservation
# Goldstein Ch.9

"""Generate ensemble cloud around (q0,p0) using Box-Muller (no external RNG)."""
function generate_ensemble(q0::Vector{Float64}, p0::Vector{Float64}, n_particles::Int; sigma_q=0.01, sigma_p=0.01)
    n=length(q0); ensemble=Vector{Tuple{Vector{Float64},Vector{Float64}}}(undef,n_particles)
    seed=12345
    lcg()=(global seed=(1103515245*seed+12345)&0x7fffffff; seed/0x7fffffff)
    function box_muller(); u1=lcg();u2=lcg(); return sqrt(-2*log(u1+1e-300))*cos(2pi*u2); end
    for k in 1:n_particles
        ensemble[k]=(q0+sigma_q*[box_muller() for _ in 1:n], p0+sigma_p*[box_muller() for _ in 1:n])
    end
    return ensemble
end

"""Liouville verification: evolve ensemble and track volume proxy."""
function verify_liouville(sys, q0::Vector{Float64}, p0::Vector{Float64}, t_end::Float64, dt::Float64; n_particles=100)
    ensemble=generate_ensemble(q0,p0,n_particles); V0=phase_volume_proxy(ensemble)
    volumes=[V0]; times=[0.0]; t=0.0; step=0
    rec_every=max(1,Int(ceil(t_end/dt/20)))
    while t<t_end
        h=min(dt,t_end-t)
        for k in 1:n_particles
            q,p=ensemble[k]; qn,pn=stormer_verlet_step(sys,q,p,h); ensemble[k]=(qn,pn)
        end
        t+=h;step+=1
        if step%rec_every==0; push!(volumes,phase_volume_proxy(ensemble)); push!(times,t); end
    end
    return times,volumes,(volumes[end]-V0)/(abs(V0)+1e-300)
end

"""Phase space density evolution: d(rho)/dt = -{rho, H} (Liouville equation)."""
function evolve_density(rho0::Function, sys, t_end::Float64, dt::Float64, q_grid::Vector{Float64}, p_grid::Vector{Float64})
    nq,np=length(q_grid),length(p_grid)
    rho=zeros(nq,np)
    for i in 1:nq,j in 1:np; rho[i,j]=rho0([q_grid[i]],[p_grid[j]]); end
    n_steps=Int(ceil(t_end/dt))
    for _ in 1:n_steps
        rho_new=copy(rho)
        for i in 2:nq-1,j in 2:np-1
            q=[q_grid[i]];p=[p_grid[j]]
            dH_dq,dH_dp=sys.grad_H(q,p)
            drho_dq=(rho[i+1,j]-rho[i-1,j])/(q_grid[2]-q_grid[1])
            drho_dp=(rho[i,j+1]-rho[i,j-1])/(p_grid[2]-p_grid[1])
            rho_new[i,j]=rho[i,j]-dt*(dH_dp[1]*drho_dq-dH_dq[1]*drho_dp)
        end
        rho=rho_new
    end
    return rho
end

"""Gibbs entropy: S = -kB * integral(rho*log(rho) dq dp). For Hamiltonian flow, dS/dt = 0."""
function gibbs_entropy(rho::Matrix{Float64}, dq::Float64, dp::Float64; kB=1.0)
    S=0.0
    for i in axes(rho,1),j in axes(rho,2)
        if rho[i,j]>0; S-=kB*rho[i,j]*log(rho[i,j])*dq*dp; end
    end
    return S
end

"""Kozlov's theorem: for integrable systems, phase space foliates into invariant tori."""
function invariant_torus_check(sys, traj::PhaseTrajectory, actions::Function)
    J_vals=[actions(q,p) for (q,p) in zip(traj.qs,traj.ps)]
    return maximum(J_vals)-minimum(J_vals)
end

"""相空间混合 (mixing): 初始区域随时间拉伸折叠. 检测: 每对粒子间距的分布."""
function mixing_proxy(ensemble::Vector{Tuple{Vector{Float64},Vector{Float64}}})
    n = length(ensemble)
    if n < 2; return 0.0; end
    phase_pts = [vcat(q,p) for (q,p) in ensemble]
    dists = Float64[]
    for i in 1:n-1, j in i+1:n
        push!(dists, norm(phase_pts[i] - phase_pts[j]))
    end
    return mean(dists), std(dists)
end

"""遍历性测试: 长时间平均 = 系综平均 (对观测函数 f)."""
function ergodicity_test(sys, q0::Vector{Float64}, p0::Vector{Float64}, t_end::Float64, dt::Float64, observable::Function)
    q, p = copy(q0), copy(p0)
    time_avg = 0.0; n_steps = Int(ceil(t_end/dt))
    for step in 1:n_steps
        q_new, p_new = stormer_verlet_step(sys, q, p, dt)
        time_avg += observable(q, p)*dt
        q, p = q_new, p_new
    end
    time_avg /= t_end
    return time_avg
end

"""庞加莱回归定理: 有限体积守恒系统中，几乎所有轨道无限次回到初始点附近."""
function poincare_recurrence_check(sys, q0::Vector{Float64}, p0::Vector{Float64}, t_max::Float64, dt::Float64, eps_val::Float64)
    q, p = copy(q0), copy(p0); t = 0.0
    while t < t_max
        q_new, p_new = stormer_verlet_step(sys, q, p, dt); t += dt
        if norm(vcat(q_new,q0) .- vcat(q0,p0)) < eps_val
            return t  # 回归时间
        end
        q, p = q_new, p_new
    end
    return Inf  # 未回归
end

"""Kolmogorov-Sinai 熵 (KS-entropy): 相空间信息产生的速率.
   对于混沌系统 KS>0; 对于可积系统 KS=0."""
function ks_entropy_estimate(sys, ensemble::Vector{Tuple{Vector{Float64},Vector{Float64}}}, t_eval::Float64, dt::Float64)
    n = length(ensemble)
    if n < 10; return 0.0; end
    # 用 Lyapunov 指数之和的近似: KS ≈ Σ λ_i^+
    total_uncertainty = 0.0
    for (q,p) in ensemble[1:min(10,n)]
        q_new, p_new = stormer_verlet_step(sys, q, p, dt)
        q_pert = q + 1e-6*randn(length(q))
        p_pert = p + 1e-6*randn(length(p))
        qp_new, pp_new = stormer_verlet_step(sys, q_pert, p_pert, dt)
        d = norm(vcat(qp_new-q_new, pp_new-p_new))
        total_uncertainty += log(d/1e-6)
    end
    return max(0.0, total_uncertainty / (min(10,n) * dt))
end
