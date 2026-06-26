# recurrence.jl — 递归分析与复杂网络: 递归图, RQA, 递归网络
# 参考: Marwan 2007, Donner 2010

"递归图: R_{ij} = Θ(ε - |x_i - x_j|)"
function recurrence_plot(points::Vector{Vector{Float64}}, epsilon::Float64)
    n=length(points); R=zeros(Int,n,n)
    for i in 1:n,j in 1:n
        norm(points[i]-points[j])<epsilon && (R[i,j]=1)
    end; R
end

"递归率 (RR): RR = (1/N²) Σ R_{ij}"
function recurrence_rate(R::Matrix{Int})
    sum(R)/length(R)
end

"确定性 (DET): 对角线结构占比"
function determinism(R::Matrix{Int}, l_min::Int=2)
    n=size(R,1); diag_count=0; diag_total=0
    for l in l_min:n; count=0
        for i in 1:n-l+1; all(R[i:i+l-1,i:i+l-1].==1) && (count+=1); end
        count>0 && (diag_total+=count*l)
    end
    for l in l_min:n
        for i in 1:n-l+1; R[i,i+l]==1 && (diag_count+=1); end
    end; diag_total/max(sum(R),1)
end

"层流性 (LAM): 垂直线结构占比"
function laminarity(R::Matrix{Int}, v_min::Int=2)
    n=size(R,1); vert_sum=0
    for j in 1:n; count=0
        for i in 1:n; R[i,j]==1 ? (count+=1) : (vert_sum+=count>=v_min?count:0; count=0); end
    end; vert_sum/max(sum(R),1)
end

"平均对角线长度 L"
function average_diagonal_length(R::Matrix{Int}, l_min::Int=2)
    n=size(R,1); total=0; count=0
    for i in 1:n; diag_len=0
        for j in 1:n-i+1; R[j,j+i-1]==1 ? (diag_len+=1) : (diag_len>=l_min && (total+=diag_len;count+=1); diag_len=0); end
        diag_len>=l_min && (total+=diag_len; count+=1)
    end; count>0 ? total/count : 0.0
end

"递归网络: 从递归图构建邻接矩阵 (不含自环)"
function recurrence_network(R::Matrix{Int})
    n=size(R,1); A=copy(R)
    for i in 1:n; A[i,i]=0; end; A
end

"递归网络度分布"
recurrence_degree_distribution(R::Matrix{Int}) = degree_distribution(recurrence_network(R))

"交叉递归图: CR_{ij} = Θ(ε - |x_i - y_j|)"
function cross_recurrence_plot(x::Vector{Vector{Float64}}, y::Vector{Vector{Float64}}, epsilon::Float64)
    nx,ny=length(x),length(y); CR=zeros(Int,nx,ny)
    for i in 1:nx,j in 1:ny
        norm(x[i]-y[j])<epsilon && (CR[i,j]=1)
    end; CR
end

"联合递归图: JR_{ij} = R1_{ij} · R2_{ij}"
function joint_recurrence_plot(R1::Matrix{Int}, R2::Matrix{Int})
    size(R1)!=size(R2) && error("Matrices must have same size")
    R1 .* R2
end

"选择 ε 使递归率为固定值 (如 5%)"
function epsilon_for_recurrence_rate(points::Vector{Vector{Float64}}, target_RR::Float64=0.05)
    n=length(points); dists=Float64[]
    for i in 1:n-1; for j in i+1:n; push!(dists,norm(points[i]-points[j])); end
    end; partialsort(dists,Int(ceil(target_RR*length(dists))))
end

export recurrence_plot, recurrence_rate, determinism, laminarity
export average_diagonal_length, recurrence_network, recurrence_degree_distribution
export cross_recurrence_plot, joint_recurrence_plot, epsilon_for_recurrence_rate
