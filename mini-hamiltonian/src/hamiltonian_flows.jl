# hamiltonian_flows.jl — Hamilton流形与几何: 辛流形, 动量映射, Lie-Poisson
# 参考: Arnold, Marsden & Ratiu

using LinearAlgebra

"辛梯度: X_f = J ∇f"
function symplectic_gradient(f::Function, q::Vector{Float64}, p::Vector{Float64}; eps_val=1e-6)
    n=length(q); X=zeros(2n)
    for i in 1:n
        qp,qm=copy(q),copy(q); qp[i]+=eps_val; qm[i]-=eps_val
        pp,pm=copy(p),copy(p); pp[i]+=eps_val; pm[i]-=eps_val
        X[i]=(f(q,pp)-f(q,pm))/(2*eps_val)  # ∂f/∂p_i
        X[n+i]=-(f(qp,p)-f(qm,p))/(2*eps_val)  # -∂f/∂q_i
    end; X
end

"辛形式: ω(ξ,η) = ξᵀ J η"
function symplectic_form(xi::Vector{Float64}, eta::Vector{Float64}, n::Int)
    J=[zeros(n,n) Matrix{Float64}(I,n,n); -Matrix{Float64}(I,n,n) zeros(n,n)]
    dot(xi, J*eta)
end

"辛矩阵 J"
symplectic_matrix(n::Int) = [zeros(n,n) Matrix{Float64}(I,n,n); -Matrix{Float64}(I,n,n) zeros(n,n)]

"正则变换的辛条件: Mᵀ J M = J"
function is_symplectic_matrix(M::Matrix{Float64}; tol=1e-10)
    n=size(M,1)÷2; J=symplectic_matrix(n)
    norm(M'*J*M - J) < tol
end

"动量映射 J: P → g* (对Lie群G的作用)"
function momentum_map_lie(group_action::Function, q::Vector{Float64}, p::Vector{Float64})
    # 对于SO(3): J = q × p (角动量)
    n=div(length(q),3); L=zeros(3)
    for i in 1:n
        qi=q[3i-2:3i]; pi=p[3i-2:3i]
        L+=cross(qi,pi)
    end; L
end

"Lie-Poisson 括号: {F,G}(μ) = -μ·[∇F,∇G]"
function lie_poisson_bracket(F::Function, G::Function, mu::Vector{Float64}, structure_constants::Array{Float64,3}; eps_val=1e-6)
    n=length(mu); gradF=zeros(n); gradG=zeros(n)
    for i in 1:n
        mp=copy(mu);mp[i]+=eps_val;mm=copy(mu);mm[i]-=eps_val
        gradF[i]=(F(mp)-F(mm))/(2*eps_val); gradG[i]=(G(mp)-G(mm))/(2*eps_val)
    end
    result=0.0
    for i in 1:n,j in 1:n,k in 1:n
        result-=mu[i]*structure_constants[i,j,k]*gradF[j]*gradG[k]
    end; result
end

"SO(3) 结构常数: [e_i,e_j] = ε_{ijk} e_k"
so3_structure_constants() = begin C=zeros(3,3,3)
    C[1,2,3]=1; C[2,3,1]=1; C[3,1,2]=1; C[2,1,3]=-1; C[3,2,1]=-1; C[1,3,2]=-1; C end

"Euler 方程的 Lie-Poisson 形式: ṁ = m × (I⁻¹ m)"
function euler_lie_poisson(m::Vector{Float64}, I_inv::Vector{Float64})
    omega = I_inv .* m
    cross(m, omega)
end

"Casimir 函数: SO(3) 的 C(μ) = |μ|²"
so3_casimir(mu::Vector{Float64}) = dot(mu,mu)

"Kostant-Kirillov 辛形式: ω_μ(ad*_ξ μ, ad*_η μ) = -μ·[ξ,η]"
function kostant_kirillov_form(mu::Vector{Float64}, xi::Vector{Float64}, eta::Vector{Float64}, C::Array{Float64,3})
    result=0.0
    for i in 1:length(mu); result-=mu[i]*dot(C[i,:,:]*xi,eta); end
    result
end

export symplectic_gradient, symplectic_form, symplectic_matrix, is_symplectic_matrix
export momentum_map_lie, lie_poisson_bracket, so3_structure_constants
export euler_lie_poisson, so3_casimir, kostant_kirillov_form
