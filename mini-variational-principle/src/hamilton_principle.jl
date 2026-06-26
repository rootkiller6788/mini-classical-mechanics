# hamilton_principle.jl — Hamilton变分原理: 最小作用量, 相空间变分, 多辛结构
# 参考: Goldstein Ch.2, Marsden & West 2001

"相空间作用量: S[q,p] = ∫ [pᵀq̇ - H(q,p)] dt"
function phase_space_action(H::Function, q::Function, qdot::Function, p::Function, t_span::Tuple{Float64,Float64}; n=500)
    t0,tf=t_span; dt=(tf-t0)/n; S=0.0
    for i in 1:n; t=t0+(i-0.5)*dt
        S+=(dot(p(t),qdot(t))-H(q(t),p(t)))*dt
    end; S
end

"变分导数在相空间: δS = 0 → q̇ = ∂H/∂p, ṗ = -∂H/∂q"
function phase_space_el_residual(H::Function, grad_H_q::Function, grad_H_p::Function, q::Vector{Float64}, p::Vector{Float64}, qdot::Vector{Float64}, pdot::Vector{Float64})
    n=length(q); res=zeros(2n)
    res[1:n]=qdot-grad_H_p(q,p)
    res[n+1:2n]=pdot+grad_H_q(q,p)
    res
end

"离散 Hamilton 原理 (变分积分器)"
function discrete_hamilton_principle(L_d::Function, qs::Vector{Vector{Float64}}, h::Float64)
    # ∂/∂q_k [L_d(q_{k-1},q_k) + L_d(q_k,q_{k+1})] = 0
    n=length(qs); residuals=Float64[]
    for k in 2:n-1
        dL1=(L_d(qs[k-1],qs[k]+1e-6,h)-L_d(qs[k-1],qs[k]-1e-6,h))/(2e-6)
        dL2=(L_d(qs[k],qs[k+1]+1e-6,h)-L_d(qs[k],qs[k+1]-1e-6,h))/(2e-6)
        push!(residuals,norm(dL1+dL2))
    end; residuals
end

"辛 2-形式检查: ω = Σ dq_i ∧ dp_i, 在演化下守恒"
function symplectic_2form(dp::Matrix{Float64}, dq::Matrix{Float64})
    # ω(ξ,η) = ξᵀ J η
    n=size(dp,1)
    dXi=[dq;dp]; J=[zeros(n,n) Matrix{Float64}(I,n,n); -Matrix{Float64}(I,n,n) zeros(n,n)]
    dXi'*J*dXi
end

"动量映射 (momentum map) 对对称群"
function momentum_map(G::Matrix{Float64}, q::Vector{Float64}, p::Vector{Float64})
    # J(q,p) = pᵀ G q  (线性对称性生成子)
    dot(p, G*q)
end

"对称性约化 (Marsden-Weinstein): 在动量映射 level set 上约化"
function symplectic_reduction(J::Float64, mu::Float64, q::Vector{Float64}, p::Vector{Float64})
    # 约束: J(q,p) = μ
    abs(J-mu) < 1e-10
end

"Hamilton 向量场 X_H = J^{-1}·∇H"
function hamiltonian_vector_field(grad_H::Vector{Float64}, n::Int)
    X = zeros(2n)
    X[1:n] = grad_H[n+1:2n]
    X[n+1:2n] = -grad_H[1:n]
    X
end

"Poisson 流形上的 Casimir: C 与所有 f 满足 {C,f}=0"
casimir_check(C::Function, H::Function, q::Vector{Float64}, p::Vector{Float64}) = poisson_bracket(C,H,q,p)

"多辛 PDE 离散 (Bridges 1997)"
function multisymplectic_residual(K::Matrix{Float64}, M::Matrix{Float64}, z::Vector{Float64}, dz_dt::Vector{Float64}, dz_dx::Vector{Float64})
    K*dz_dt + M*dz_dx
end

export phase_space_action, phase_space_el_residual
export discrete_hamilton_principle, symplectic_2form
export momentum_map, symplectic_reduction
export hamiltonian_vector_field, casimir_check, multisymplectic_residual
