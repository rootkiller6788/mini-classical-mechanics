# numerical_variational.jl — 变分问题的数值方法: FEM, 谱方法, 优化
# 参考: Reddy, Johnson "Numerical Solution of PDEs by FEM"

using LinearAlgebra

"FEM 1D 刚度矩阵: K_e = ∫ Bᵀ D B dx"
function fem_stiffness_1d(x_nodes::Vector{Float64}, D_func::Function)
    n=length(x_nodes); K=zeros(n,n)
    for e in 1:n-1
        h=x_nodes[e+1]-x_nodes[e]; x_mid=(x_nodes[e]+x_nodes[e+1])/2
        D=D_func(x_mid)
        Ke=(D/h)*[1 -1; -1 1]
        K[e:e+1,e:e+1]+=Ke
    end; K
end

"FEM 质量矩阵: M_e = ∫ Nᵀ ρ N dx"
function fem_mass_1d(x_nodes::Vector{Float64}, rho_func::Function)
    n=length(x_nodes); M=zeros(n,n)
    for e in 1:n-1
        h=x_nodes[e+1]-x_nodes[e]; x_mid=(x_nodes[e]+x_nodes[e+1])/2
        rho=rho_func(x_mid)
        Me=(rho*h/6)*[2 1; 1 2]
        M[e:e+1,e:e+1]+=Me
    end; M
end

"线性基函数在单元上的积分"
function integrate_over_element(f::Function, h::Float64, n_gauss::Int=2)
    # 2-点 Gauss 积分
    gauss_pts=[-1/sqrt(3),1/sqrt(3)]; gauss_w=[1.0,1.0]
    result=0.0
    for (xi,w) in zip(gauss_pts,gauss_w)
        x_mapped=0.5*h*(xi+1)  # [0,h]
        result+=w*f(x_mapped)
    end; 0.5*h*result
end

"Rayleigh-Ritz 法: min (1/2 uᵀKu - uᵀf)"
function rayleigh_ritz_solve(K::Matrix{Float64}, f::Vector{Float64}, bc_dofs::Dict{Int,Float64})
    free_dofs=setdiff(1:length(f),collect(keys(bc_dofs)))
    K_ff=K[free_dofs,free_dofs]; f_f=f[free_dofs]
    K_fb=K[free_dofs,collect(keys(bc_dofs))]; u_b=[bc_dofs[i] for i in sort(collect(keys(bc_dofs)))]
    u_f=K_ff\(f_f-K_fb*u_b)
    u=zeros(length(f)); u[free_dofs]=u_f; u[collect(keys(bc_dofs))]=u_b
    u
end

"特征值问题: K u = λ M u (简正模)"
function generalized_eigenvalue(K::Matrix{Float64}, M::Matrix{Float64}, n_modes::Int=5)
    evals,evecs=eigen(K,M)
    idx=sortperm(real.(evals))[1:min(n_modes,length(evals))]
    real.(evals[idx]), evecs[:,idx]
end

"谱 Galerkin 法 (用 Fourier 基)"
function spectral_galerkin_fourier(operator::Function, n_modes::Int, L::Float64)
    K=zeros(n_modes,n_modes)
    for i in 1:n_modes, j in 1:n_modes
        # ∫₀ˡ sin(iπx/L) L[sin(jπx/L)] dx
        K[i,j]=quadrature_sin_product(operator,i,j,L)
    end
    K
end

function quadrature_sin_product(operator::Function, i::Int, j::Int, L::Float64; n_pts=100)
    dx=L/n_pts; result=0.0
    for k in 1:n_pts
        x=(k-0.5)*dx
        result+=sin(i*pi*x/L)*operator(j,x,L)*dx
    end; result
end

"有限差分 (变分推导的离散化)"
function variational_fd_laplacian(n::Int, dx::Float64)
    K=zeros(n,n)
    for i in 2:n-1; K[i,i-1:i+1]=[-1,2,-1]/dx^2; end
    K[1,1:2]=[1,-1]/dx^2; K[n,n-1:n]=[-1,1]/dx^2; K
end

"变分不等式 (障碍问题) 的投影算法"
function obstacle_problem_solve(K::Matrix{Float64}, f::Vector{Float64}, psi::Vector{Float64}; n_iter=1000)
    n=length(f); u=zeros(n)
    for _ in 1:n_iter
        for i in 1:n
            u[i]=max(psi[i], sum(K[i,j]*u[j] for j in 1:n if j!=i)/K[i,i])
        end
    end; u
end

export fem_stiffness_1d, fem_mass_1d, integrate_over_element
export rayleigh_ritz_solve, generalized_eigenvalue
export spectral_galerkin_fourier, variational_fd_laplacian, obstacle_problem_solve
# 附加函数 — 变分数值方法扩展
# 自适应网格, 误差估计, 变分多尺度

"Zienkiewicz-Zhu 误差估计 (超收敛应力恢复)"
function zienkiewicz_zhu_error_estimate(sigma_fe::Vector{Float64}, sigma_recovered::Vector{Float64})
    n=length(sigma_fe); eta=0.0; U=0.0
    for e in 1:n
        diff=sigma_fe[e]-sigma_recovered[e]
        eta+=diff^2; U+=sigma_recovered[e]^2
    end; sqrt(eta)/sqrt(max(U,1e-15))
end

"h-细化 (单元细分) 准则"
function h_refinement_criterion(error::Float64, tolerance::Float64, h::Float64, p::Int)
    error > tolerance ? h*0.5 : h
end

"p-细化 (阶次提升) 准则"
function p_refinement_criterion(smoothness::Float64)
    smoothness > 0.8 ? 2 : 1
end

"变分多尺度 (VMS) 稳定化参数: τ = h²/(4κ)"
vms_stabilization_parameter(h::Float64, kappa::Float64) = h^2/(4*kappa)

"Hughes 变分多尺度 (对流-扩散)"
function vms_advection_diffusion(Pe::Float64, h::Float64)
    tau = h/(2*max(Pe,1e-10))  # SUPG 参数
    tau
end

"Morse 指数估计 (共轭梯度法)"
function morse_index_estimate(H::Matrix{Float64})
    evals = eigvals(H)
    count(evals .< -1e-10)
end

export zienkiewicz_zhu_error_estimate, h_refinement_criterion, p_refinement_criterion
export vms_stabilization_parameter, vms_advection_diffusion, morse_index_estimate
"变分问题的 Newton-Kantorovich 线性化"
function newton_kantorovich_step(F::Function, dF::Function, u::Vector{Float64}; tol=1e-8, max_iter=20)
    for _ in 1:max_iter
        R = F(u); J = dF(u)
        du = J \ (-R); u += du
        norm(du) < tol && break
    end; u
end

"自适应变分积分器 (根据局部误差调整步长)"
function adaptive_variational_stepsize(error_est::Float64, tol::Float64, h::Float64)
    safety = 0.9
    h_new = h * safety * (tol/max(error_est,1e-15))^(1/3)
    clamp(h_new, 0.1*h, 2.0*h)
end

export newton_kantorovich_step, adaptive_variational_stepsize
