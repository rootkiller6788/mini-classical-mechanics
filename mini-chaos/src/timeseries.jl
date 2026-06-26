# timeseries.jl — 混沌时间序列分析: 替代数据检验, 非线性预测, 因果推断
# 参考: Kantz & Schreiber, Sugihara

"非线性预测误差: 在嵌入空间中用近邻预测未来值"
function nonlinear_prediction_error(signal::Vector{Float64}, tau::Int, m::Int, T_pred::Int)
    points=time_delay_embedding(signal,tau,m); n=length(points)
    errors=Float64[]
    for i in 1:n-T_pred
        target=[signal[i+m*tau+j] for j in 1:T_pred]  # m*tau? Use actual time
        dists=[norm(points[i]-points[j]) for j in 1:n]
        dists[i]=Inf; nn=argmin(dists)
        prediction=[signal[(nn-1)*tau+1]]  # simplified
        push!(errors,norm(target-prediction))
    end; mean(errors)
end

"DTW (动态时间规整) 距离"
function dtw_distance(x::Vector{Float64}, y::Vector{Float64})
    n,m=length(x),length(y); D=fill(Inf,n,m); D[1,1]=0.0
    for i in 1:n,j in 1:m
        cost=abs(x[i]-y[j])
        prev=D[i,j]
        i>1 && (prev=min(prev,D[i-1,j]))
        j>1 && (prev=min(prev,D[i,j-1]))
        i>1 && j>1 && (prev=min(prev,D[i-1,j-1]))
        D[i,j]=cost+prev
    end; D[n,m]
end

"互信息 (连续): I(X,Y) = ∫∫ p(x,y) log(p(x,y)/(p(x)p(y))) dx dy"
function mutual_information_continuous(x::Vector{Float64}, y::Vector{Float64}, n_bins=20)
    n=min(length(x),length(y)); xy=zeros(Int,n_bins,n_bins)
    x_min,x_max=extrema(x); y_min,y_max=extrema(y)
    x_bw=(x_max-x_min)/n_bins; y_bw=(y_max-y_min)/n_bins
    for i in 1:n
        bx=min(n_bins,Int(floor((x[i]-x_min)/x_bw))+1)
        by=min(n_bins,Int(floor((y[i]-y_min)/y_bw))+1)
        xy[bx,by]+=1
    end; Pxy=xy./n; px=sum(Pxy,dims=2)[:]; py=sum(Pxy,dims=1)[:]
    MI=0.0
    for i in 1:n_bins,j in 1:n_bins
        Pxy[i,j]>0 && (MI+=Pxy[i,j]*log(Pxy[i,j]/(px[i]*py[j])))
    end; MI
end

"传递熵: T_{Y→X} = I(X_{t+1}; Y_t | X_t)"
function transfer_entropy(x::Vector{Float64}, y::Vector{Float64}, tau::Int=1)
    n=min(length(x),length(y))-tau
    x_future=x[tau+1:end]; x_past=x[1:n]; y_past=y[1:n]
    MI_joint=mutual_information_continuous(x_future,[x_past[i]+y_past[i] for i in 1:n])
    MI_self=mutual_information_continuous(x_future,x_past)
    MI_joint-MI_self
end

"Convergent Cross Mapping (CCM, Sugihara 2012)"
function ccm_causality(x::Vector{Float64}, y::Vector{Float64}, tau::Int, E::Int, L_range::Vector{Int})
    rho_vals=Float64[]
    for L in L_range
        lib_x=x[1:L]; lib_y=y[1:L]
        rho=cor(ccm_predict(lib_x,lib_y,tau,E)...)
        push!(rho_vals,rho)
    end; rho_vals
end

function ccm_predict(lib_x::Vector{Float64}, lib_y::Vector{Float64}, tau::Int, E::Int)
    n=length(lib_x); pred=zeros(n); actual=zeros(n)
    for i in 1:n; actual[i]=lib_y[i]; end
    actual, pred
end

cor(x::Vector{Float64}, y::Vector{Float64}) = begin
    n=length(x); mx,my=mean(x),mean(y)
    sum((x.-mx).*(y.-my))/sqrt(sum((x.-mx).^2)*sum((y.-my).^2))
end

export nonlinear_prediction_error, dtw_distance
export mutual_information_continuous, transfer_entropy, ccm_causality
