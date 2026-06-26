# embedding.jl — 吸引子重构: 时间延迟嵌入, 假近邻, 互信息
# 参考: Takens 1981, Kantz & Schreiber

"时间延迟嵌入: X(t) = [x(t), x(t-τ), x(t-2τ), ..., x(t-(m-1)τ)]"
function time_delay_embedding(signal::Vector{Float64}, tau::Int, m::Int)
    n=length(signal); N=n-(m-1)*tau; points=Vector{Vector{Float64}}(undef,N)
    for i in 1:N; pts=Float64[]
        for j in 0:m-1; push!(pts,signal[i+j*tau]); end
        points[i]=pts
    end; points
end

"假近邻法 (False Nearest Neighbors) 确定嵌入维数 m"
function false_nearest_neighbors(signal::Vector{Float64}, tau::Int, m_max::Int; R_tol=10.0, A_tol=2.0)
    fnn_ratio=Float64[]
    for m in 1:m_max
        points=time_delay_embedding(signal,tau,m); N=length(points)
        fnn_count=0
        for i in 1:N
            dists=[norm(points[i]-points[j]) for j in 1:N]; dists[i]=Inf
            nn_idx=argmin(dists); dist_m=dists[nn_idx]
            dist_m1=norm([points[i];signal[i+m*tau]]-[points[nn_idx];signal[nn_idx+m*tau]])
            if dist_m > 0 && (dist_m1/dist_m > R_tol || dist_m1/dist_m > A_tol)
                fnn_count+=1
            end
        end
        push!(fnn_ratio,fnn_count/N)
    end; fnn_ratio
end

"平均互信息 (Average Mutual Information) 确定延迟 τ"
function average_mutual_information(signal::Vector{Float64}, tau::Int, n_bins::Int=20)
    n=length(signal)-tau; min_x, max_x = minimum(signal), maximum(signal)
    bin_width = (max_x-min_x)/n_bins
    hist2d=zeros(Int,n_bins,n_bins)
    for i in 1:n
        bx=min(n_bins,Int(floor((signal[i]-min_x)/bin_width))+1)
        by=min(n_bins,Int(floor((signal[i+tau]-min_x)/bin_width))+1)
        hist2d[bx,by]+=1
    end; P=hist2d./n; AMI=0.0
    px=sum(P,dims=2)[:]; py=sum(P,dims=1)[:]
    for i in 1:n_bins,j in 1:n_bins
        P[i,j]>0 && (AMI+=P[i,j]*log(P[i,j]/(px[i]*py[j])))
    end; AMI
end

"首个 AMI 极小值作为最佳延迟 τ"
function find_optimal_tau(signal::Vector{Float64}, tau_max::Int=50)
    amis=[average_mutual_information(signal,t) for t in 1:tau_max]
    for t in 2:tau_max-1
        amis[t-1]>amis[t] && amis[t]<amis[t+1] && return t
    end; argmin(amis)
end

"关联维数 (Grassberger-Procaccia)"
function correlation_dimension(points::Vector{Vector{Float64}}, epsilons::Vector{Float64})
    N=length(points); Cs=Float64[]
    for eps_val in epsilons; C=0.0; count=0
        for i in 1:N-1
            for j in i+1:N
                norm(points[i]-points[j]) < eps_val && (C+=1.0)
            end
        end; push!(Cs,2*C/(N*(N-1)))
    end
    x=log.(epsilons); y=log.(Cs)
    n=length(x); (n*sum(x.*y)-sum(x)*sum(y))/(n*sum(x.^2)-sum(x)^2)
end

"替代数据 (Surrogate data) 检验非线性"
function random_shuffle_surrogate(signal::Vector{Float64})
    s=copy(signal); shuffle!(s); s
end

"AAFT (Amplitude Adjusted Fourier Transform) 替代数据"
function aaft_surrogate(signal::Vector{Float64})
    n=length(signal); s_sorted=sort(signal)
    gaussian=randn(n); g_sorted=sort(gaussian)
    ranks=[searchsortedfirst(s_sorted,x) for x in signal]
    surrogate=[g_sorted[r] for r in ranks]
    surrogate
end

export time_delay_embedding, false_nearest_neighbors
export average_mutual_information, find_optimal_tau
export correlation_dimension, random_shuffle_surrogate, aaft_surrogate
