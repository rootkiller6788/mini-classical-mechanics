# poisson.jl — Poisson bracket algebra & canonical invariants
# Goldstein Ch.9

"""Poisson bracket {f,g} with numerical gradients."""
function poisson_bracket(f::Function, g::Function, q::Vector{Float64}, p::Vector{Float64}; eps_val::Float64=1e-6)
    n=length(q); result=0.0
    for i in 1:n
        qp=copy(q);qp[i]+=eps_val;qm=copy(q);qm[i]-=eps_val
        pp=copy(p);pp[i]+=eps_val;pm=copy(p);pm[i]-=eps_val
        result+=((f(qp,p)-f(qm,p))*(g(q,pp)-g(q,pm)) - (f(q,pp)-f(q,pm))*(g(qp,p)-g(qm,p)))/(4eps_val^2)
    end
    return result
end

"""Fundamental Poisson bracket errors: {qi,qj},{pi,pj},{qi,pj}-delta_ij."""
function fundamental_poisson_brackets(n::Int, q::Vector{Float64}, p::Vector{Float64})
    err_qq=zeros(n,n);err_pp=zeros(n,n);err_qp=zeros(n,n)
    for i in 1:n,j in 1:n
        qi(qq,pp)=qq[i];qj(qq,pp)=qq[j];pi(qq,pp)=pp[i];pj(qq,pp)=pp[j]
        err_qq[i,j]=abs(poisson_bracket(qi,qj,q,p));err_pp[i,j]=abs(poisson_bracket(pi,pj,q,p))
        err_qp[i,j]=abs(poisson_bracket(qi,pj,q,p)-(i==j?1.0:0.0))
    end
    return err_qq,err_pp,err_qp
end

"""Test if f is a constant of motion: {f,H}=0."""
function is_constant_of_motion(f::Function, H::Function, q::Vector{Float64}, p::Vector{Float64})
    return abs(poisson_bracket(f,H,q,p))<1e-8
end

"""Angular momentum z-component: Lz = x*py - y*px."""
function angular_momentum_poisson_z(q::Vector{Float64}, p::Vector{Float64})
    return q[1]*p[2]-q[2]*p[1]
end

"""Verify Jacobi identity: {f,{g,h}} + {g,{h,f}} + {h,{f,g}} = 0."""
function verify_jacobi_identity(f::Function, g::Function, h::Function, q::Vector{Float64}, p::Vector{Float64})
    gh(qq,pp)=poisson_bracket(g,h,qq,pp);hf(qq,pp)=poisson_bracket(h,f,qq,pp);fg(qq,pp)=poisson_bracket(f,g,qq,pp)
    return poisson_bracket(f,gh,q,p)+poisson_bracket(g,hf,q,p)+poisson_bracket(h,fg,q,p)
end

"""Equation of motion via Poisson: df/dt = {f, H}."""
time_derivative_via_poisson(f::Function,H::Function,q::Vector{Float64},p::Vector{Float64}) = poisson_bracket(f,H,q,p)

"""Poisson bracket in canonical coordinates (analytic form)."""
function canonical_poisson(f::Function,g::Function,n::Int,q::Vector{Float64},p::Vector{Float64};eps_val=1e-6)
    return poisson_bracket(f,g,q,p,eps_val=eps_val)
end

"""Lie bracket of two Hamiltonian vector fields: [X_f, X_g] = X_{-{f,g}}."""
function lie_bracket_hamiltonian(f::Function,g::Function,H::Function,q::Vector{Float64},p::Vector{Float64})
    pb_fg=poisson_bracket(f,g,q,p)
    pb_negfg_H=poisson_bracket((qq,pp)->-pb_fg,H,q,p)
    commutator = poisson_bracket(f,(qq,pp)->poisson_bracket(g,H,qq,pp),q,p) - poisson_bracket(g,(qq,pp)->poisson_bracket(f,H,qq,pp),q,p)
    return commutator, pb_negfg_H, abs(commutator-pb_negfg_H)<1e-6
end

"""Poisson bracket in Darboux coordinates (where symplectic form omega = sum dq_i wedge dp_i)."""
function darboux_poisson_matrix(n::Int)
    J=zeros(2n,2n)
    for i in 1:n; J[i,n+i]=1.0; J[n+i,i]=-1.0; end
    return J
end

"""Gradient of a function on phase space."""
function phase_gradient(f::Function,q::Vector{Float64},p::Vector{Float64};eps_val=1e-6)
    n=length(q);grad=zeros(2n)
    for i in 1:n
        qp=copy(q);qp[i]+=eps_val;qm=copy(q);qm[i]-=eps_val
        grad[i]=(f(qp,p)-f(qm,p))/(2eps_val)
        pp=copy(p);pp[i]+=eps_val;pm=copy(p);pm[i]-=eps_val
        grad[n+i]=(f(q,pp)-f(q,pm))/(2eps_val)
    end
    return grad
end

"""Poisson bracket geometric: {f,g} = ∇f^T · J · ∇g where J = [[0,I],[-I,0]]."""
function poisson_bracket_geometric(f::Function, g::Function, q::Vector{Float64}, p::Vector{Float64}; eps_val=1e-6)
    grad_f = phase_gradient(f, q, p; eps_val=eps_val)
    grad_g = phase_gradient(g, q, p; eps_val=eps_val)
    n = length(q)
    result = 0.0
    for i in 1:n
        result += grad_f[i]*grad_g[n+i] - grad_f[n+i]*grad_g[i]
    end
    return result
end

"""Canonical Poisson bracket for quadratic Hamiltonians: {q_i q_j, p_k p_l} etc."""
function quadratic_poisson_example(q::Vector{Float64}, p::Vector{Float64})
    n = length(q)
    f(qq,pp) = qq[1]^2
    g(qq,pp) = pp[1]^2
    pb = poisson_bracket(f, g, q, p)
    return pb  # Should be 4 q1 p1
end

"""角动量 Poisson 括号: {L_i, L_j} = ε_{ijk} L_k (SO(3) 代数)."""
function angular_momentum_algebra_check(q::Vector{Float64}, p::Vector{Float64}; tol=1e-8)
    n = div(length(q), 3)
    Lx(qq,pp) = sum(qq[3*(i-1)+2]*pp[3*(i-1)+3] - qq[3*(i-1)+3]*pp[3*(i-1)+2] for i in 1:n)
    Ly(qq,pp) = sum(qq[3*(i-1)+3]*pp[3*(i-1)+1] - qq[3*(i-1)+1]*pp[3*(i-1)+3] for i in 1:n)
    Lz(qq,pp) = sum(qq[3*(i-1)+1]*pp[3*(i-1)+2] - qq[3*(i-1)+2]*pp[3*(i-1)+1] for i in 1:n)
    pb_xy = poisson_bracket(Lx, Ly, q, p)
    pb_yz = poisson_bracket(Ly, Lz, q, p)
    pb_zx = poisson_bracket(Lz, Lx, q, p)
    return (LX_LY=pb_xy - Lz(q,p), LY_LZ=pb_yz - Lx(q,p), LZ_LX=pb_zx - Ly(q,p),
            valid = all(abs.([pb_xy-Lz(q,p), pb_yz-Lx(q,p), pb_zx-Ly(q,p)]) .< tol))
end

"""Casimir 函数: 与所有函数 Poisson 交换 = 0 (如 L² = Lx²+Ly²+Lz²)."""
function is_casimir(C::Function, q::Vector{Float64}, p::Vector{Float64}; n_test=10, tol=1e-8)
    n = length(q)
    for _ in 1:n_test
        f_test(qq,pp) = sum(randn(n).*qq) + sum(randn(n).*pp)
        if abs(poisson_bracket(C, f_test, q, p)) > tol
            return false
        end
    end
    return true
end

"""Poisson 括号与正则变换的相容性: {f,g}_{(q,p)} = {f∘φ, g∘φ}_{(Q,P)}."""
function canonical_invariance_of_poisson(f::Function, g::Function, ct, q::Vector{Float64}, p::Vector{Float64}; tol=1e-8)
    pb_old = poisson_bracket(f, g, q, p)
    Q_func(qq,pp) = ct.forward(qq,pp)[1]; P_func(qq,pp) = ct.forward(qq,pp)[2]
    f_new(qq,pp) = f(Q_func(qq,pp), P_func(qq,pp))
    g_new(qq,pp) = g(Q_func(qq,pp), P_func(qq,pp))
    pb_new = poisson_bracket(f_new, g_new, q, p)
    return abs(pb_old - pb_new) < tol
end
