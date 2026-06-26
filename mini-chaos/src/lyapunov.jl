# lyapunov.jl — Lyapunov 指数: 两粒子法, Benettin全谱, 条件指数
# 参考: Benettin 1980, Wolf 1985, Strogatz Ch.9

using LinearAlgebra

"两粒子法 — 最大 Lyapunov 指数 (连续系统)"
function lyapunov_exponent(f::Function, x0::Vector{Float64}, t_end::Float64, dt::Float64; d0=1e-8)
    n=length(x0); xA=copy(x0); xB=x0+d0*randn_manual(n)
    n_steps=Int(ceil(t_end/dt)); log_sum=0.0; count=0
    for _ in 1:n_steps
        xAn=rk4_flow_step(f,xA,dt); xBn=rk4_flow_step(f,xB,dt)
        d=norm(xBn-xAn)
        if d>1e-15; log_sum+=log(d/d0); count+=1; xB=xAn+(d0/d)*(xBn-xAn); end
        xA=xAn
    end
    count>0 ? log_sum/(count*dt) : 0.0
end

"Wolf 算法 — 最大 Lyapunov 指数 (更稳定)"
function lyapunov_wolf(f::Function, x0::Vector{Float64}, t_end::Float64, dt::Float64; d0=1e-8, evol_steps=10)
    n=length(x0); x=copy(x0)
    pert = randn_manual(n); pert = pert/norm(pert)*d0
    x_pert = x + pert
    log_sum = 0.0; count = 0; total_steps = Int(ceil(t_end/dt))
    for step in 1:total_steps
        x = rk4_flow_step(f, x, dt)
        x_pert = rk4_flow_step(f, x_pert, dt)
        if step % evol_steps == 0
            d = norm(x_pert - x)
            if d > 1e-15
                log_sum += log(d/d0); count += 1
                x_pert = x + (x_pert - x)*(d0/d)
            end
        end
    end
    count>0 ? log_sum/(count*evol_steps*dt) : 0.0
end

"Benettin 全 Lyapunov 谱"
function lyapunov_spectrum_benettin(f::Function, x0::Vector{Float64}, t_end::Float64, dt::Float64; n_transient=500, n_renorm=5)
    n=length(x0); x=copy(x0)
    W=Matrix{Float64}(I,n,n); sums=zeros(n); count=0
    n_steps=Int(ceil(t_end/dt))
    for step in 1:n_steps
        x_new=rk4_flow_step(f,x,dt)
        for j in 1:n
            pert=W[:,j]*1e-8
            xp=rk4_flow_step(f,x+pert,dt); xm=rk4_flow_step(f,x-pert,dt)
            W[:,j]=(xp-xm)/(2e-8)
        end
        if step%n_renorm==0
            Q,R=qr(W); W=Matrix(Q)
            if step>n_transient
                for j in 1:n; sums[j]+=log(abs(R[j,j])); end; count+=1
            end
        end
        x=x_new
    end
    lyaps=[s/(count*dt*n_renorm) for s in sums]; sort!(lyaps,rev=true)
    sum_pos=sum(lyaps[lyaps.>0])
    kyd=lyaps[1]>0 ? 1+(lyaps[1]+(length(lyaps)>1?lyaps[2]:-1))/abs(length(lyaps)>1?lyaps[2]:1) : 0.0
    LyapunovSpectrum(lyaps,sum_pos,kyd)
end

"条件 Lyapunov 指数 (驱动-响应同步)"
function conditional_lyapunov(f_drive::Function, f_response::Function, x0::Vector{Float64}, y0::Vector{Float64}, t_end::Float64, dt::Float64; d0=1e-8)
    n=length(y0); yA=copy(y0); yB=y0+d0*randn_manual(n)
    log_sum=0.0; count=0; n_steps=Int(ceil(t_end/dt))
    for _ in 1:n_steps
        x_drive=rk4_flow_step(f_drive,x0,dt); x0=x_drive
        yA_new=rk4_flow_step(s->f_response(s,x_drive),yA,dt)
        yB_new=rk4_flow_step(s->f_response(s,x_drive),yB,dt)
        d=norm(yB_new-yA_new)
        if d>1e-15; log_sum+=log(d/d0); count+=1; yB=yA_new+(d0/d)*(yB_new-yA_new); end
        yA=yA_new
    end
    count>0 ? log_sum/(count*dt) : 0.0
end

"Logistic map Lyapunov: λ=lim(1/n)Σlog|r(1-2x_n)|"
function logistic_lyapunov(r::Float64; n_iter=10000, n_transient=500)
    x=0.5; for _ in 1:n_transient; x=logistic_map(x,r); end
    lam=0.0
    for _ in 1:n_iter; x=logistic_map(x,r); lam+=log(abs(r*(1-2x))); end
    lam/n_iter
end

"Lyapunov vs r (Logistic map 参数扫描)"
function logistic_lyapunov_scan(r_range::AbstractRange; n_iter=5000)
    [logistic_lyapunov(r; n_iter=n_iter) for r in r_range]
end

"Hénon map Lyapunov 对"
function henon_lyapunovs(a=1.4, b=0.3; n_iter=100000, n_transient=1000)
    x,y=0.5,0.5; for _ in 1:n_transient; x,y=henon_map(x,y,a,b); end
    lam1=lam2=0.0
    for _ in 1:n_iter
        J=[-2a*x 1; b 0]; Q,R=qr(J)
        lam1+=log(abs(R[1,1])); lam2+=log(abs(R[2,2])); x,y=henon_map(x,y,a,b)
    end
    lam1/n_iter, lam2/n_iter
end

"手动 LCG 随机数 (零外部依赖)"
function randn_manual(n::Int)
    seed=12345; lcg()=(seed=(1103515245*seed+12345)&0x7fffffff; seed/0x7fffffff)
    [sqrt(-2*log(lcg()+1e-300))*cos(2pi*lcg()) for _ in 1:n]
end

export lyapunov_exponent, lyapunov_wolf, lyapunov_spectrum_benettin
export conditional_lyapunov, logistic_lyapunov, logistic_lyapunov_scan, henon_lyapunovs
