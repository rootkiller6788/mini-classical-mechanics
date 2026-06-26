# networks.jl — 耦合网络与时空混沌: CML, 小世界, 无标度网络
# 参考: Kaneko, Watts-Strogatz, Barabasi-Albert

"Coupled Map Lattice (CML): x_i(t+1) = (1-ε)f(x_i(t)) + ε/2[f(x_{i-1})+f(x_{i+1})]"
function coupled_map_lattice(f::Function, x0::Vector{Float64}, epsilon::Float64, n_iter::Int)
    n=length(x0); x=copy(x0); lattice=zeros(n_iter+1,n)
    lattice[1,:]=x0
    for t in 1:n_iter
        x_new=zeros(n)
        for i in 1:n
            left=mod1(i-1,n); right=mod1(i+1,n)
            x_new[i]=(1-epsilon)*f(x[i])+0.5*epsilon*(f(x[left])+f(x[right]))
        end; x=x_new; lattice[t+1,:]=x
    end; lattice
end

"耦合强度临界值 (空间均匀定态失稳)"
function cml_critical_coupling(f_prime::Float64)
    1-f_prime
end

"Watts-Strogatz 小世界网络生成"
function watts_strogatz_network(n::Int, k::Int, p::Float64)
    # ring lattice with k neighbors per side, rewiring probability p
    adj=zeros(Int,n,n)
    for i in 1:n
        for j in 1:k
            neighbor=mod1(i+j,n)
            if rand()<p
                new_neighbor=rand(1:n)
                while new_neighbor==i || adj[i,new_neighbor]==1
                    new_neighbor=rand(1:n)
                end; adj[i,new_neighbor]=1; adj[new_neighbor,i]=1
            else; adj[i,neighbor]=1; adj[neighbor,i]=1
            end
        end
    end; adj
end

"Barabasi-Albert 无标度网络 (优先连接)"
function barabasi_albert_network(n::Int, m::Int)
    adj=zeros(Int,n,n)
    # 初始完全图 m0=m+1
    for i in 1:m+1
        for j in i+1:m+1; adj[i,j]=1; adj[j,i]=1; end
    end; degrees=sum(adj,dims=2)[:]
    for new_node in m+2:n
        targets=Int[]
        total_deg=sum(degrees)
        while length(targets)<m
            for candidate in 1:new_node-1
                candidate in targets && continue
                rand() < degrees[candidate]/total_deg && push!(targets,candidate)
                length(targets)==m && break
            end
        end
        for t in targets; adj[new_node,t]=1; adj[t,new_node]=1; end
        degrees=sum(adj,dims=2)[:]
    end; adj
end

"网络度分布"
function degree_distribution(adj::Matrix{Int})
    deg=sum(adj,dims=2)[:]
    deg_count=Dict{Int,Int}()
    for d in deg; deg_count[d]=get(deg_count,d,0)+1; end
    deg_count
end

"聚类系数"
function clustering_coefficient(adj::Matrix{Int})
    n=size(adj,1); C=zeros(n)
    for i in 1:n
        neighbors=findall(x->x>0,adj[i,:])
        k=length(neighbors); k<2 && continue
        triangles=0
        for a in 1:k-1
            for b in a+1:k; triangles+=adj[neighbors[a],neighbors[b]]; end
        end; C[i]=2*triangles/(k*(k-1))
    end; mean(C)
end

"Kuramoto 模型 + 网络拓扑"
function kuramoto_on_network(theta::Vector{Float64}, omega::Vector{Float64}, adj::Matrix{Int}, K::Float64)
    n=length(theta); dtheta=zeros(n); deg=sum(adj,dims=2)[:]
    for i in 1:n; s=0.0
        for j in 1:n; adj[i,j]>0 && (s+=sin(theta[j]-theta[i])); end
        dtheta[i]=omega[i]+K/deg[i]*s
    end; dtheta
end

export coupled_map_lattice, cml_critical_coupling
export watts_strogatz_network, barabasi_albert_network
export degree_distribution, clustering_coefficient, kuramoto_on_network
