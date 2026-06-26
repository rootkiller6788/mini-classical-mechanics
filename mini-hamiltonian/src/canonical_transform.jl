# canonical_transform.jl — Canonical transformations & generating functions
# Goldstein Ch.9

struct CanonicalTransform
    name::String; transform::Function; inverse::Function; generating_function::Function
end

"""F2(q,P) = q·P generates identity: Q=q, P unchanged."""
function identity_transform(n::Int)
    F2(q,P)=dot(q,P)
    CanonicalTransform("identity",(q,p)->(copy(q),copy(p)),(Q,P)->(copy(Q),copy(P)),F2)
end

"""Point transform Q=f(q), P=(df/dq)^{-T}·p."""
function point_transform(f::Function, inv_f::Function, jac_f::Function, n::Int)
    CanonicalTransform("point",(q,p)->(f(q),inv(jac_f(q)')*p),(Q,P)->(inv_f(Q),jac_f(inv_f(Q))'*P),(q,Q)->0.0)
end

"""Scale: Q=mu*q, P=nu*p with mu*nu=1."""
function scale_transform(mu::Float64, n::Int)
    nu=1/mu
    CanonicalTransform("scale_$mu",(q,p)->(mu*q,nu*p),(Q,P)->(Q/mu,P/nu),(q,Q)->0.0)
end

"""Exchange: Q=p, P=-q."""
function exchange_transform(n::Int)
    CanonicalTransform("exchange",(q,p)->(copy(p),-copy(q)),(Q,P)->(-copy(P),copy(Q)),(q,Q)->dot(q,Q))
end

"""Verify canonical condition via fundamental Poisson brackets."""
function verify_canonical(ct::CanonicalTransform, q_test::Vector{Float64}, p_test::Vector{Float64}; tol=1e-8)
    Q_test,P_test=ct.transform(q_test,p_test); n=length(q_test)
    err_qq,err_pp,err_qp=fundamental_poisson_brackets(n,Q_test,P_test)
    return max(maximum(err_qq),maximum(err_pp),maximum(err_qp))<tol, maximum([err_qq;err_pp;err_qp])
end

"""Symplectic Jacobian check: J_new = M·J_old·M^T where M = d(Q,P)/d(q,p)."""
function symplectic_jacobian_check(ct::CanonicalTransform, q_test::Vector{Float64}, p_test::Vector{Float64}; eps_val=1e-6)
    n=length(q_test); M=zeros(2n,2n)
    Q0,P0=ct.transform(q_test,p_test)
    for j in 1:n
        qp=copy(q_test);qp[j]+=eps_val;qm=copy(q_test);qm[j]-=eps_val
        Qp,Pp=ct.transform(qp,p_test);Qm,Pm=ct.transform(qm,p_test)
        for i in 1:n
            M[i,j]=(Qp[i]-Qm[i])/(2eps_val); M[n+i,j]=(Pp[i]-Pm[i])/(2eps_val)
        end
        pp=copy(p_test);pp[j]+=eps_val;pm=copy(p_test);pm[j]-=eps_val
        Qp2,Pp2=ct.transform(q_test,pp);Qm2,Pm2=ct.transform(q_test,pm)
        for i in 1:n
            M[i,n+j]=(Qp2[i]-Qm2[i])/(2eps_val); M[n+i,n+j]=(Pp2[i]-Pm2[i])/(2eps_val)
        end
    end
    J_old=darboux_poisson_matrix(n)
    return M*J_old*M', J_old, maximum(abs.(M*J_old*M'-J_old))
end

"""Generating function F2(q,P) type: given q and p, compute (Q,P) via F2."""
function solve_F2_transform(F2::Function, grad_q_F2::Function, grad_P_F2::Function, q::Vector{Float64}, p::Vector{Float64}; tol=1e-10)
    n=length(q)
    # p = grad_q F2(q,P) → solve for P
    P_guess=copy(p)  # initial guess
    for _ in 1:20
        F_val=grad_q_F2(q,P_guess)-p
        if norm(F_val)<tol; break; end
        # dF/dP = grad_P(grad_q F2) — use numerical Hessian
        eps_val=1e-6; H=zeros(n,n)
        for j in 1:n; Pp=copy(P_guess);Pp[j]+=eps_val;Pm=copy(P_guess);Pm[j]-=eps_val
            H[:,j]=(grad_q_F2(q,Pp)-grad_q_F2(q,Pm))/(2eps_val); end
        P_guess-=H\F_val
    end
    Q=grad_P_F2(q,P_guess)
    return Q,P_guess
end

"""Compose two canonical transforms: CT3 = CT2 ∘ CT1."""
function compose_transforms(ct1::CanonicalTransform, ct2::CanonicalTransform)
    CanonicalTransform("$(ct1.name)∘$(ct2.name)",
        (q,p)->ct2.transform(ct1.transform(q,p)...),
        (Q,P)->ct1.inverse(ct2.inverse(Q,P)...),
        (q,Q)->0.0)
end
