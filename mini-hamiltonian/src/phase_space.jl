# phase_space.jl — 相空间结构与分析工具
# 参考：Goldstein Ch.8-9

struct PhasePoint
    n::Int; q::Vector{Float64}; p::Vector{Float64}
end
PhasePoint(q::Vector{Float64}, p::Vector{Float64}) = PhasePoint(length(q), q, p)

mutable struct PhaseTrajectory
    ts::Vector{Float64}; qs::Vector{Vector{Float64}}; ps::Vector{Vector{Float64}}; H_vals::Vector{Float64}
end
PhaseTrajectory() = PhaseTrajectory(Float64[], Vector{Float64}[], Vector{Float64}[], Float64[])

function record!(traj::PhaseTrajectory, t::Float64, q::Vector{Float64}, p::Vector{Float64}, H_val::Float64)
    push!(traj.ts, t); push!(traj.qs, q); push!(traj.ps, p); push!(traj.H_vals, H_val)
end

"""Phase space volume element (proxy via covariance determinant)."""
function phase_volume_proxy(ensemble::Vector{Tuple{Vector{Float64},Vector{Float64}}})
    n_p=length(ensemble); n_d=length(ensemble[1][1])
    qm=sum([e[1] for e in ensemble])/n_p; pm=sum([e[2] for e in ensemble])/n_p
    qs=max([norm(q-qm) for q in [e[1] for e in ensemble]])
    ps=max([norm(p-pm) for p in [e[2] for e in ensemble]])
    return qs*ps
end

"""Poincare section: crossings of surface_func(q,p)=0."""
function poincare_section(traj::PhaseTrajectory, surface_func::Function)
    points=PhasePoint[]
    for i in 2:length(traj.ts)
        s_prev=surface_func(traj.qs[i-1],traj.ps[i-1])
        s_curr=surface_func(traj.qs[i],traj.ps[i])
        if s_prev*s_curr<0
            frac=abs(s_prev)/(abs(s_prev)+abs(s_curr))
            push!(points,PhasePoint(traj.qs[i-1]+frac*(traj.qs[i]-traj.qs[i-1]), traj.ps[i-1]+frac*(traj.ps[i]-traj.ps[i-1])))
        end
    end
    return points
end

"""Phase space vector field (qdot, pdot) on a grid for visualization."""
function phase_vector_field(sys, q_range::Tuple{Float64,Float64}, p_range::Tuple{Float64,Float64}, n_grid::Int=20)
    qs=range(q_range[1],q_range[2],length=n_grid)
    ps=range(p_range[1],p_range[2],length=n_grid)
    field=zeros(n_grid,n_grid,2)
    for (iq,q) in enumerate(qs), (ip,p) in enumerate(ps)
        dH_dq,dH_dp=sys.grad_H([q],[p])
        field[iq,ip,1]=dH_dp[1]; field[iq,ip,2]=-dH_dq[1]
    end
    return qs,ps,field
end

"""Fixed points: where (qdot,pdot)=(0,0), i.e. grad_H=0. Newton-Raphson search."""
function find_fixed_points(sys, q_range::Tuple{Float64,Float64}, p_range::Tuple{Float64,Float64}; n_seeds=10, tol=1e-10)
    seeds=[]; qs=range(q_range[1],q_range[2],length=n_seeds); ps=range(p_range[1],p_range[2],length=n_seeds)
    for q in qs, p in ps; push!(seeds,[q,p]); end
    fixed=Float64[]
    for (q0,p0) in seeds
        q,p=q0,p0
        for _ in 1:50
            dH_dq,dH_dp=sys.grad_H([q],[p])
            if abs(dH_dp[1])<tol && abs(dH_dq[1])<tol
                push!(fixed,q); break
            end
            # Newton step (simplified: gradient descent)
            q+=0.1*dH_dp[1]; p-=0.1*dH_dq[1]
        end
    end
    return unique([round(f,digits=8) for f in fixed])
end

"""Stability of fixed point: eigenvalues of Jacobian of Hamilton's equations."""
function fixed_point_stability(sys, q_fixed::Float64, p_fixed::Float64; eps_val=1e-6)
    # Jacobian J = [∂²H/∂p∂q  ∂²H/∂p²; -∂²H/∂q²  -∂²H/∂q∂p]
    dH_dq_p,dH_dp_p=sys.grad_H([q_fixed+eps_val],[p_fixed])
    dH_dq_m,dH_dp_m=sys.grad_H([q_fixed-eps_val],[p_fixed])
    d2H_dq2=(dH_dq_p[1]-dH_dq_m[1])/(2eps_val)
    d2H_dpdq=(dH_dp_p[1]-dH_dp_m[1])/(2eps_val)
    dH_dq_p2,dH_dp_p2=sys.grad_H([q_fixed],[p_fixed+eps_val])
    dH_dq_m2,dH_dp_m2=sys.grad_H([q_fixed],[p_fixed-eps_val])
    d2H_dqdp=(dH_dq_p2[1]-dH_dq_m2[1])/(2eps_val)
    d2H_dp2=(dH_dp_p2[1]-dH_dp_m2[1])/(2eps_val)
    J=[d2H_dpdq d2H_dp2; -d2H_dq2 -d2H_dqdp]
    evals=eigvals(J)
    return evals, maximum(abs.(real(evals)))<1e-10 ? :elliptic : :hyperbolic
end

"""Moment map: J: phase space → R^k for symmetry group action."""
function momentum_map(momentum_funcs::Vector{Function}, q::Vector{Float64}, p::Vector{Float64})
    return [f(q,p) for f in momentum_funcs]
end
