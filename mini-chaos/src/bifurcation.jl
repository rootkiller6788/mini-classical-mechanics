# bifurcation.jl — 分岔分析: Logistic分岔图, Feigenbaum常数, 周期检测, 分岔分类
# 参考: Strogatz Ch.3, Ch.8, Kuznetsov "Elements of Applied Bifurcation Theory"

# ============================================================
# Logistic 分岔图
# ============================================================
"Logistic map 分岔图数据 (r范围, 去瞬态)"
function logistic_bifurcation(r_range::AbstractRange, n_transient=500, n_plot=100)
    rs,xs=Float64[],Float64[]
    for r in r_range; x=0.5
        for _ in 1:n_transient; x=logistic_map(x,r); end
        for _ in 1:n_plot; x=logistic_map(x,r); push!(rs,r); push!(xs,x); end
    end
    rs, xs
end

"通用 1D 映射的分岔图"
function bifurcation_diagram(f::Function, param_range::AbstractRange, x0::Float64; n_transient=500, n_plot=100)
    params,x_vals=Float64[],Float64[]
    for p in param_range; x=x0
        for _ in 1:n_transient; x=f(x,p); end
        for _ in 1:n_plot; x=f(x,p); push!(params,p); push!(x_vals,x); end
    end
    params, x_vals
end

"2D映射分岔图 (取某个变量的极值)"
function bifurcation_diagram_2d(f_map::Function, param_range::AbstractRange, x0::Vector{Float64}; n_transient=500, n_plot=100, var_idx=1)
    params,y_vals=Float64[],Float64[]
    for p in param_range; x=copy(x0)
        for _ in 1:n_transient; x=f_map(x,p); end
        for _ in 1:n_plot; x=f_map(x,p); push!(params,p); push!(y_vals,x[var_idx]); end
    end
    params, y_vals
end

# ============================================================
# 周期检测
# ============================================================
"检测离散映射的周期 (返回周期数, 0=未检测到)"
function detect_period_discrete(f::Function, x0::Float64; n_iter=2000, tol=1e-6)
    x=x0; visited=Float64[]
    for _ in 1:n_iter; x=f(x)
        for (j,v) in enumerate(visited); abs(x-v)<tol && return length(visited)-j+1; end
        push!(visited,x)
    end
    0
end

"检测在参数 r 处的周期"
function find_period_at_r(r::Float64, f::Function; n_iter=2000, n_transient=1000)
    x=0.5
    for _ in 1:n_transient; x=f(x,r); end
    detect_period_discrete(p->f(p,r), x; n_iter=n_iter)
end

"检测连续系统周期 (Poincare 截面)"
function detect_period_continuous(traj::Vector{Vector{Float64}}; tol=1e-4)
    x_last=traj[end]
    for i in 1:length(traj)-1
        if norm(traj[i]-x_last) < tol; return length(traj)-i; end
    end
    0
end

# ============================================================
# Feigenbaum 常数
# ============================================================
const FEIGENBAUM_DELTA = 4.669201609102990
const FEIGENBAUM_ALPHA = 2.502907875095892

"找周期倍增分岔点序列"
function find_bifurcation_points(r_start=2.8, n_periods=5; f=logistic_map)
    r_vals=Float64[]; r=r_start; dr=0.01
    for p in 1:n_periods
        while true
            per=find_period_at_r(r,f); per==2^p && (push!(r_vals,r); break); r+=dr
        end; dr/=5
    end
    r_vals
end

"从分岔点序列估计 Feigenbaum δ"
function estimate_feigenbaum(r_vals::Vector{Float64})
    deltas=Float64[]
    for i in 2:length(r_vals)-1
        push!(deltas, (r_vals[i]-r_vals[i-1])/(r_vals[i+1]-r_vals[i]))
    end
    deltas
end

"估计 Feigenbaum α (分岔宽度比)"
function estimate_feigenbaum_alpha(x_bif::Vector{Float64})
    alphas=Float64[]
    for i in 2:length(x_bif)-1
        push!(alphas, (x_bif[i]-x_bif[i-1])/(x_bif[i+1]-x_bif[i]))
    end
    alphas
end

# ============================================================
# 分岔分类
# ============================================================
"Saddle-node: 1个不动点出现→消失; Pitchfork: 1→3; Hopf: 极限环"
function classify_bifurcation(dx_dr::Float64, d2x_dx2::Float64, d2x_drdx::Float64)
    if abs(dx_dr) < 1e-10 && abs(d2x_dx2) > 1e-10
        :saddle_node
    elseif abs(dx_dr) < 1e-10 && abs(d2x_dx2) < 1e-10 && abs(d2x_drdx) > 1e-10
        :transcritical
    elseif abs(dx_dr) < 1e-10 && abs(d2x_dx2) < 1e-10 && abs(d2x_drdx) < 1e-10
        :pitchfork
    else
        :hopf
    end
end

"1D 系统的不动点稳定性 (f(x,r)=0): λ = ∂f/∂x"
function fixed_point_stability(f::Function, x::Float64, r::Float64; dx=1e-6)
    df_dx = (f(x+dx,r)-f(x-dx,r))/(2*dx)
    abs(df_dx) < 1.0 ? :stable : :unstable
end

"不动点求解 (1D): f(x,r)=0, Newton法"
function find_fixed_point(f::Function, r::Float64, x0::Float64; tol=1e-10, max_iter=50)
    x=x0
    for _ in 1:max_iter
        fx=f(x,r); dfx=(f(x+1e-6,r)-f(x-1e-6,r))/(2e-6)
        abs(dfx) < 1e-15 && break
        dx = -fx/dfx; x += dx; abs(dx) < tol && break
    end
    x
end

# ============================================================
# 相图辅助
# ============================================================
"线性稳定性矩阵 (流系统)"
function jacobian_flow(f::Function, x::Vector{Float64}; dx=1e-6)
    n=length(x); J=zeros(n,n); f0=f(x)
    for j in 1:n
        xp=copy(x); xp[j]+=dx; fp=f(xp)
        for i in 1:n; J[i,j]=(fp[i]-f0[i])/dx; end
    end
    J
end

"平衡点类型 (从Jacobian特征值)"
function equilibrium_type(eigenvalues::Vector{ComplexF64})
    re=real.(eigenvalues); im=imag.(eigenvalues)
    if all(re.<0) && all(im.==0)
        :stable_node
    elseif all(re.<0) && any(im.!=0)
        :stable_focus
    elseif all(re.>0) && all(im.==0)
        :unstable_node
    elseif all(re.>0) && any(im.!=0)
        :unstable_focus
    elseif any(re.>0) && any(re.<0)
        :saddle
    elseif all(abs.(re).<1e-10) && any(im.!=0)
        :center
    else
        :degenerate
    end
end

export logistic_bifurcation, bifurcation_diagram, bifurcation_diagram_2d
export detect_period_discrete, find_period_at_r, detect_period_continuous
export FEIGENBAUM_DELTA, FEIGENBAUM_ALPHA
export find_bifurcation_points, estimate_feigenbaum, estimate_feigenbaum_alpha
export classify_bifurcation, fixed_point_stability, find_fixed_point
export jacobian_flow, equilibrium_type
