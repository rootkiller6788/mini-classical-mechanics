# control.jl — 混沌控制: OGY方法, 延迟反馈, 参数扰动
# 参考: Ott Grebogi Yorke 1990, Pyragas 1992

"OGY (Ott-Grebogi-Yorke) 方法: 在不动点附近线性化, 用参数扰动稳定"
function ogy_control(map_func::Function, x0::Float64, p_nominal::Float64, p_range::Tuple{Float64,Float64}, n_iter::Int)
    x=map_func(x0,p_nominal)
    for _ in 1:n_iter
        # 估计不动点附近的 Jacobian
        x_fp=find_fixed_point(r->map_func(1.0,r), p_nominal, x0)  # 不动点
        # 用 OGY 反馈: δp = -K·(x - x_fp)
        K=0.5; delta_x=x-x_fp
        p_corrected=p_nominal+clamp(K*delta_x,p_range[1]-p_nominal,p_range[2]-p_nominal)
        x=map_func(x,p_corrected)
    end; x
end

"Pyragas 延迟反馈控制 (Time-Delayed Feedback)"
function pyragas_control(f_ode::Function, x0::Vector{Float64}, t_end::Float64, dt::Float64, K::Float64, tau::Float64)
    n=length(x0); x=copy(x0); x_tau=copy(x0)
    delay_steps=Int(round(tau/dt)); history=[copy(x0) for _ in 1:delay_steps+1]
    traj=[copy(x0)]
    for step in 1:Int(ceil(t_end/dt))
        x_tau=history[max(1,step-delay_steps+1)]
        dx=f_ode(x)+K*(x-x_tau)
        x=x+dt*dx; push!(traj,copy(x))
        history[mod1(step,delay_steps+1)]=copy(x)
    end; traj
end

"参数共振抑制 (parametric resonance suppression)"
function parametric_chaos_control(map_func::Function, x0::Float64, p0::Float64, amplitude::Float64, omega::Float64, n_iter::Int)
    x=x0; traj=[x0]
    for i in 1:n_iter
        p=p0+amplitude*sin(omega*i)
        x=map_func(x,p); push!(traj,x)
    end; traj
end

"线性反馈控制: u = -K(x - x_target)"
function linear_feedback_control(f_ode::Function, x0::Vector{Float64}, x_target::Vector{Float64}, K::Float64, t_end::Float64, dt::Float64)
    x=copy(x0)
    for _ in 1:Int(ceil(t_end/dt))
        dx=f_ode(x)-K*(x-x_target); x+=dt*dx
    end; x
end

"自适应控制 (参数估计 + 反馈)"
function adaptive_control(f_ode::Function, x0::Vector{Float64}, t_end::Float64, dt::Float64, gamma::Float64)
    x=copy(x0); n=length(x0)
    K_hat=0.0; x_ref=1.0
    for _ in 1:Int(ceil(t_end/dt))
        e=x[1]-x_ref
        K_hat += dt*gamma*e^2
        dx=f_ode(x)-(K_hat+1.0)*(x[1]-x_ref)
        x+=dt*dx
    end; x,K_hat
end

"稳定流形计算 (Stable manifold of a saddle)"
function stable_manifold(map_func::Function, x_saddle::Vector{Float64}, p::Float64, n_points::Int)
    n=length(x_saddle); J=zeros(n,n)
    f(x)=map_func(x,p)
    for i in 1:n; xp=copy(x_saddle);xp[i]+=1e-5;xm=copy(x_saddle);xm[i]-=1e-5
        J[:,i]=(f(xp)-f(xm))/(2e-5)
    end
    evals,evecs=eigen(J); stable_idx=argmin(abs.(evals))
    v_stable=real(evecs[:,stable_idx])
    manifold=[x_saddle+v_stable*0.01*(i-n_points/2)/(n_points/2) for i in 1:n_points]
    manifold
end

export ogy_control, pyragas_control, parametric_chaos_control
export linear_feedback_control, adaptive_control, stable_manifold
