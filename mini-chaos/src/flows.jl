# flows.jl — 连续时间混沌系统: Lorenz, Rossler, Chua, Duffing, Sprott
# 参考: Strogatz Ch.9, Sprott "Chaos and Time-Series Analysis"

using LinearAlgebra

"Lorenz: dx/dt=σ(y-x), dy/dt=x(ρ-z)-y, dz/dt=xy-βz"
function lorenz(state::Vector{Float64}, sigma=10.0, rho=28.0, beta=8/3)
    x,y,z=state[1],state[2],state[3]
    [sigma*(y-x), x*(rho-z)-y, x*y-beta*z]
end

"Rossler: dx/dt=-y-z, dy/dt=x+ay, dz/dt=b+z(x-c)"
function rossler(state::Vector{Float64}, a=0.2, b=0.2, c=5.7)
    x,y,z=state[1],state[2],state[3]
    [-y-z, x+a*y, b+z*(x-c)]
end

"Chua电路: dx/dt=α(y-x-f(x)), dy/dt=x-y+z, dz/dt=-βy"
function chua(state::Vector{Float64}, alpha=15.6, beta=28.0, m0=-1.143, m1=-0.714)
    x,y,z=state[1],state[2],state[3]
    fx=m1*x+0.5*(m0-m1)*(abs(x+1)-abs(x-1))
    [alpha*(y-x-fx), x-y+z, -beta*y]
end

"Duffing: ẋ=v, v̇=-δv+βx-αx³+γcos(ωt), ṫ=1"
function duffing(state::Vector{Float64}, delta=0.2, beta=1.0, alpha=1.0, gamma=0.3, omega=1.0)
    x,v,t=state[1],state[2],state[3]
    [v, -delta*v+beta*x-alpha*x^3+gamma*cos(omega*t), 1.0]
end

"Van der Pol: ẋ=v, v̇=μ(1-x²)v-x+Fcos(ωt)"
function forced_vdp(state::Vector{Float64}, mu=3.0, F=1.2, omega=1.0)
    x,v,t=state[1],state[2],state[3]
    [v, mu*(1-x^2)*v-x+F*cos(omega*t), 1.0]
end

"Sprott B: dx/dt=yz, dy/dt=x-y, dz/dt=1-xy"
sprott_b(state::Vector{Float64}) = let x,y,z=state[1],state[2],state[3]; [y*z, x-y, 1-x*y] end

"Sprott C: dx/dt=yz, dy/dt=x-y, dz/dt=1-x²"
sprott_c(state::Vector{Float64}) = let x,y,z=state[1],state[2],state[3]; [y*z, x-y, 1-x^2] end

"Chen系统: dx/dt=a(y-x), dy/dt=(c-a)x-xz+cy, dz/dt=xy-bz"
function chen_system(state::Vector{Float64}, a=35.0, b=3.0, c=28.0)
    x,y,z=state[1],state[2],state[3]; [a*(y-x), (c-a)*x-x*z+c*y, x*y-b*z]
end

"Simplest quadratic chaotic flow: dx/dt=y, dy/dt=z, dz/dt=-x+yz"
function simplest_quadratic(state::Vector{Float64})
    x,y,z=state[1],state[2],state[3]; [y, z, -x+y*z]
end

"Hénon-Heiles (Hamiltonian): ẋ=p_x, ẏ=p_y, ṗ_x=-x-2xy, ṗ_y=-y-x²+y²"
function henon_heiles(state::Vector{Float64})
    x,y,px,py=state[1],state[2],state[3],state[4]
    [px, py, -x-2*x*y, -y-x^2+y^2]
end

"4D hyperchaos (Rössler): dx/dt=-y-z, dy/dt=x+0.25y+w, dz/dt=3+xz, dw/dt=-0.5z+0.05w"
function hyperchaos_rossler(state::Vector{Float64})
    x,y,z,w=state[1],state[2],state[3],state[4]
    [-y-z, x+0.25*y+w, 3+x*z, -0.5*z+0.05*w]
end

"RK4 步"
function rk4_flow_step(f::Function, state::Vector{Float64}, dt::Float64)
    half=dt/2; k1=f(state); k2=f(state+half*k1)
    k3=f(state+half*k2); k4=f(state+dt*k3)
    state+(dt/6)*(k1+2k2+2k3+k4)
end

"RK4 轨迹积分"
function integrate_flow(f::Function, x0::Vector{Float64}, t_end::Float64, dt::Float64)
    n=Int(ceil(t_end/dt)); traj=[copy(x0)]; x=copy(x0)
    for _ in 1:n; x=rk4_flow_step(f,x,dt); push!(traj,copy(x)); end
    traj
end

"Poincaré 截面: 记录每次穿过截面的事件"
function poincare_section(f::Function, x0::Vector{Float64}, t_end::Float64, dt::Float64, section_dim::Int=2)
    traj=integrate_flow(f,x0,t_end,dt)
    points=Vector{Float64}[]
    for i in 2:length(traj)
        if traj[i-1][section_dim]*traj[i][section_dim] < 0
            push!(points, traj[i])
        end
    end
    points
end

export lorenz, rossler, chua, duffing, forced_vdp, sprott_b, sprott_c
export chen_system, simplest_quadratic, henon_heiles, hyperchaos_rossler
export rk4_flow_step, integrate_flow, poincare_section
