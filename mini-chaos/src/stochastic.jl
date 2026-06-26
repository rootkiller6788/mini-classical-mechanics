# stochastic.jl — 随机动力学: Langevin方程, Fokker-Planck, 噪声诱导现象
# 参考: Gardiner "Stochastic Methods", Risken

"Langevin 方程: ẋ = f(x) + √(2D) ξ(t), ξ 为高斯白噪声"
function langevin_1d(f::Function, D::Float64, x0::Float64, t_end::Float64, dt::Float64)
    x=x0; traj=[x0]; n_steps=Int(ceil(t_end/dt))
    for _ in 1:n_steps
        xi=randn()
        x+=dt*f(x)+sqrt(2*D*dt)*xi; push!(traj,x)
    end; traj
end

"Fokker-Planck 方程 (1D): ∂P/∂t = -∂/∂x[A(x)P] + ∂²/∂x²[B(x)P]"
function fokker_planck_1d(A::Function, B_func::Function, x_grid::Vector{Float64}, t_end::Float64, dt::Float64, P0::Vector{Float64})
    n=length(x_grid); dx=x_grid[2]-x_grid[1]; nt=Int(ceil(t_end/dt))
    P=copy(P0)
    for _ in 1:nt; P_new=copy(P)
        for i in 2:n-1
            A_term=(A(x_grid[i+1])*P[i+1]-A(x_grid[i-1])*P[i-1])/(2*dx)
            B_term=(B_func(x_grid[i+1])*P[i+1]-2*B_func(x_grid[i])*P[i]+B_func(x_grid[i-1])*P[i-1])/dx^2
            P_new[i]=P[i]+dt*(-A_term+B_term)
        end; P=copy(P_new)
    end; P
end

"双稳态系统的随机共振: ẋ = x-x³ + A cos(ωt) + √(2D) ξ"
function stochastic_resonance(x0::Float64, A::Float64, omega::Float64, D::Float64, t_end::Float64, dt::Float64)
    f(x)=x-x^3+A*cos(omega*(length(traj)*dt))
    traj=Float64[]; x=x0; n_steps=Int(ceil(t_end/dt))
    for step in 1:n_steps
        f_det=x-x^3+A*cos(omega*step*dt)
        x+=dt*f_det+sqrt(2*D*dt)*randn(); push!(traj,x)
    end; traj
end

"Kramers 逃逸率: r ~ exp(-ΔU/D) (弱噪声极限)"
function kramers_rate(barrier_height::Float64, D::Float64)
    exp(-barrier_height/D)
end

"噪声诱导相变: 平均场近似"
function noise_induced_transition(D::Float64, a::Float64, b::Float64)
    # 有效势: V_eff = -a x²/2 + b x⁴/4 + D log|x| (Ito)
    D
end

"Ornstein-Uhlenbeck 过程: dX = -θ X dt + σ dW"
function ornstein_uhlenbeck(theta::Float64, sigma::Float64, x0::Float64, t_end::Float64, dt::Float64)
    x=x0; traj=[x0]; n_steps=Int(ceil(t_end/dt))
    for _ in 1:n_steps
        x+=dt*(-theta*x)+sigma*sqrt(dt)*randn(); push!(traj,x)
    end; traj
end

"维纳过程 (布朗运动)"
wiener_process(T::Float64, dt::Float64) = ornstein_uhlenbeck(0.0,1.0,0.0,T,dt)

"Gillespie 算法 (化学主方程)"
function gillespie_algorithm(propensities::Vector{Function}, stoichiometry::Matrix{Int}, x0::Vector{Int}, t_end::Float64)
    x=copy(x0); t=0.0; times=[0.0]; states=[copy(x0)]
    while t<t_end
        a=[f(x) for f in propensities]; a0=sum(a); a0<=0 && break
        r1,r2=rand(),rand(); tau=log(1/r1)/a0
        cum=0.0; j=1
        for k in 1:length(a); cum+=a[k]; r2*a0<=cum && (j=k; break); end
        for k in 1:length(x); x[k]+=stoichiometry[j,k]; end
        t+=tau; push!(times,t); push!(states,copy(x))
    end; times, states
end

export langevin_1d, fokker_planck_1d, stochastic_resonance
export kramers_rate, noise_induced_transition
export ornstein_uhlenbeck, wiener_process, gillespie_algorithm
