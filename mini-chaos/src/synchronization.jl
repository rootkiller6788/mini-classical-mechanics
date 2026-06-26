# synchronization.jl — 混沌同步: Pecora-Carroll, Kuramoto模型, 相同步
# 参考: Pikovsky Rosenblum Kurths, Strogatz

"Pecora-Carroll 驱动-响应同步"
function pecorra_carroll_sync(f_drive::Function, f_coupled::Function, x0::Vector{Float64}, y0::Vector{Float64}, t_end::Float64, dt::Float64)
    x=copy(x0); y=copy(y0); sync_errors=Float64[]
    n_steps=Int(ceil(t_end/dt))
    for _ in 1:n_steps
        x=rk4_flow_step(f_drive,x,dt); y=rk4_flow_step(s->f_coupled(s,x),y,dt)
        push!(sync_errors,norm(x[1:min(length(x),length(y))]-y))
    end; sync_errors
end

"Kuramoto 模型: θ̇_i = ω_i + (K/N) Σ sin(θ_j-θ_i)"
function kuramoto_model(theta::Vector{Float64}, omega::Vector{Float64}, K::Float64)
    N=length(theta); dtheta=zeros(N)
    for i in 1:N; s=0.0
        for j in 1:N; s+=sin(theta[j]-theta[i]); end
        dtheta[i]=omega[i]+K/N*s
    end; dtheta
end

"Kuramoto 序参量: r = |(1/N) Σ exp(iθ_j)|"
function kuramoto_order_parameter(theta::Vector{Float64})
    N=length(theta); z=complex(0.0,0.0)
    for th in theta; z+=exp(im*th); end
    abs(z)/N, angle(z)
end

"临界耦合强度: K_c = 2/(π g(0)) (g为频率分布密度)"
function kuramoto_critical_coupling(g0::Float64)
    2/(pi*g0)
end

"Lorenz 系统间的同步误差 (两个耦合 Lorenz)"
function coupled_lorenz_sync(sigma=10.0, rho=28.0, beta=8/3, coupling=1.0, t_end=50.0, dt=0.01)
    x1=[1.0,1.0,1.0]; x2=[1.1,0.9,0.9]
    errors=Float64[]
    for _ in 1:Int(ceil(t_end/dt))
        dx1=lorenz(x1,sigma,rho,beta); dx2=lorenz(x2,sigma,rho,beta)
        dx2[1]+=coupling*(x1[1]-x2[1])
        x1+=dt*dx1; x2+=dt*dx2
        push!(errors,norm(x1-x2))
    end; errors
end

"相同步: 检测 m:n 锁相 (|mΔφ₁-nΔφ₂| < const)"
function phase_locking_detection(phase1::Vector{Float64}, phase2::Vector{Float64}, m::Int, n::Int; tol=0.1)
    delta=(m*phase1 .- n*phase2) .% (2pi)
    std(delta) < tol
end

"Rössler 系统的相位 (从YZ平面投影)"
function rossler_phase(y::Float64, z::Float64)
    atan(y,z)
end

"同步流形稳定性 (Lyapunov 条件Lyapunov指数)"
function synchronization_stability_lyapunov(f_drive::Function, f_response::Function, x0::Vector{Float64}, y0::Vector{Float64}, t_end::Float64, dt::Float64)
    # 横向 Lyapunov 指数 < 0 → 稳定同步
    conditional_lyapunov(f_drive,f_response,x0,y0,t_end,dt)
end

"Chimera 态检测: 部分同步, 部分非同步"
function chimera_state_detection(theta::Vector{Float64})
    # 计算局部序参量: 相邻粒子的小团簇
    n=length(theta); local_order=zeros(n); r_neighbor=3
    for i in 1:n
        z=complex(0.0); count=0
        for j in max(1,i-r_neighbor):min(n,i+r_neighbor)
            z+=exp(im*theta[j]); count+=1
        end; local_order[i]=abs(z)/count
    end
    # chimera: 某些区域 r~1 (同步), 某些区域 r~0 (去同步)
    mean(local_order), std(local_order)
end

export pecorra_carroll_sync, kuramoto_model, kuramoto_order_parameter
export kuramoto_critical_coupling, coupled_lorenz_sync
export phase_locking_detection, rossler_phase
export synchronization_stability_lyapunov, chimera_state_detection
