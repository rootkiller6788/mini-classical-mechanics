# field_theory.jl — 场论变分原理: Ginzburg-Landau, 相场, 弹性场
# 参考: Landau Vol.5, Ginzburg & Landau 1950

using LinearAlgebra

# ============================================================
# Ginzburg-Landau 自由能
# ============================================================
"""
Ginzburg-Landau 自由能泛函 (1D):
F[ψ] = ∫ [a|ψ|² + b|ψ|⁴/2 + c|∂ψ/∂x|²] dx

a = α(T-Tc), b,c > 0

Euler-Lagrange: c ∂²ψ/∂x² - aψ - b|ψ|²ψ = 0
"""
function ginzburg_landau_free_energy(psi::Function, dpsidx::Function, a::Float64, b::Float64, c::Float64, x_range::Tuple{Float64,Float64}; n=200)
    x0,x1=x_range; dx=(x1-x0)/n; F=0.0
    for i in 1:n
        x=x0+(i-0.5)*dx; p=psi(x); dp=dpsidx(x)
        F += (a*p^2 + 0.5*b*p^4 + c*dp^2)*dx
    end; F
end

"""
Ginzburg-Landau 方程的数值解 (松弛法)
"""
function solve_ginzburg_landau(a::Float64, b::Float64, c::Float64, L::Float64, bc_left::Float64, bc_right::Float64; nx=500, n_iter=5000, omega=1.5)
    dx=L/(nx-1); psi=zeros(nx)
    for i in 1:nx; psi[i]=0.5*(bc_left+bc_right); end
    psi[1]=bc_left; psi[nx]=bc_right
    for _ in 1:n_iter
        psi_old=copy(psi)
        for i in 2:nx-1
            psi[i]=(1-omega)*psi[i]+omega*( (psi[i+1]+psi[i-1])/dx^2 )/( 2/dx^2 + a/c + b*psi[i]^2/c )
        end
        norm(psi-psi_old)<1e-10 && break
    end; range(0,L,length=nx), psi
end

"""
相干长度: ξ = √(c/|a|). 描述序参量的空间变化尺度
"""
gl_coherence_length(a::Float64, c::Float64) = sqrt(c/abs(a))

"""
临界磁场: Hc² = 4π a²/b (type-I超导体热力学临界场)
"""
gl_critical_field(a::Float64, b::Float64) = sqrt(4*pi*a^2/b)

# ============================================================
# 相场模型 (Phase Field)
# ============================================================
"""
Allen-Cahn 方程: ∂φ/∂t = -M δF/δφ = M(ε²∇²φ - f'(φ))

F[φ] = ∫ [ε²/2|∇φ|² + f(φ)] dx
f(φ) = (φ²-1)²/4 (双井势)
"""
function allen_cahn_rhs(phi::Vector{Float64}, epsilon::Float64, M::Float64, dx::Float64)
    n=length(phi); dphi=zeros(n)
    for i in 2:n-1
        laplacian = (phi[i+1]-2*phi[i]+phi[i-1])/dx^2
        f_prime = phi[i]*(phi[i]^2 - 1)  # ∂/∂φ[(φ²-1)²/4] = φ(φ²-1)
        dphi[i] = M*(epsilon^2*laplacian - f_prime)
    end; dphi
end

"""
Cahn-Hilliard 方程 (保守序参量): ∂φ/∂t = ∇·[M∇(δF/δφ)]
"""
function cahn_hilliard_rhs(phi::Vector{Float64}, epsilon::Float64, M::Float64, dx::Float64)
    n=length(phi); mu=zeros(n); dphi=zeros(n)
    for i in 2:n-1
        f_prime = phi[i]*(phi[i]^2 - 1)
        mu[i] = -epsilon^2*(phi[i+1]-2*phi[i]+phi[i-1])/dx^2 + f_prime
    end
    for i in 3:n-2
        dphi[i] = M*(mu[i+1]-2*mu[i]+mu[i-1])/dx^2
    end; dphi
end

"界面能: γ = ∫ ε²/2|∇φ|² + f(φ) dx"
function interface_energy(phi::Vector{Float64}, epsilon::Float64, dx::Float64)
    n=length(phi); E=0.0
    for i in 2:n
        grad = (phi[i]-phi[i-1])/dx; phi_mid = 0.5*(phi[i]+phi[i-1])
        f_val = (phi_mid^2 - 1)^2/4
        E += (0.5*epsilon^2*grad^2 + f_val)*dx
    end; E
end

# ============================================================
# 弹性场变分原理
# ============================================================
"Helmholtz 自由能: Ψ(ε,T) = ½ε:C:ε - (3λ+2μ)α(T-T₀)trε"
function helmholtz_free_energy(strain::Vector{Float64}, C::Matrix{Float64}, alpha::Float64, delta_T::Float64, lam::Float64, mu::Float64)
    0.5*dot(strain, C*strain) - (3*lam+2*mu)*alpha*delta_T*sum(strain[1:3])
end

"Gibbs 自由能: G(σ,T) = -½σ:S:σ"
function gibbs_free_energy(stress::Vector{Float64}, S::Matrix{Float64})
    -0.5*dot(stress, S*stress)
end

# ============================================================
# 电磁场变分原理
# ============================================================
"""
电磁场作用量 (Lorentz 不变):
S[A_μ] = -¼∫ F_{μν}F^{μν} d⁴x + ∫ j^μ A_μ d⁴x

Euler-Lagrange → Maxwell 方程: ∂_μ F^{μν} = j^ν
"""
function em_action_density(E::Vector{Float64}, B::Vector{Float64}, rho::Float64, J::Vector{Float64}, A::Vector{Float64}, phi::Float64)
    # -¼ F² + j·A
    -0.25*(dot(E,E) - dot(B,B)) - rho*phi + dot(J, A)
end

export ginzburg_landau_free_energy, solve_ginzburg_landau, gl_coherence_length, gl_critical_field
export allen_cahn_rhs, cahn_hilliard_rhs, interface_energy
export helmholtz_free_energy, gibbs_free_energy
