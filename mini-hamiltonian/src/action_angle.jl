# action_angle.jl — Action-angle variables for integrable systems
# Goldstein Ch.10, Arnold Ch.10

struct HarmonicOscillatorActionAngle; m::Float64; omega0::Float64; end
struct PendulumActionAngle; m::Float64; L::Float64; g::Float64; end

"""SHO: J=E/omega, energy_from_action: E=omega*J."""
action_from_energy_ho(aa::HarmonicOscillatorActionAngle, E::Float64) = E/aa.omega0
energy_from_action_ho(aa::HarmonicOscillatorActionAngle, J::Float64) = aa.omega0*J

function action_to_qp_ho(aa::HarmonicOscillatorActionAngle, J::Float64, theta::Float64)
    amp=sqrt(2J/(aa.m*aa.omega0)); return amp*sin(theta), aa.m*aa.omega0*amp*cos(theta)
end

function qp_to_action_angle_ho(aa::HarmonicOscillatorActionAngle, q::Float64, p::Float64)
    E=p^2/(2aa.m)+0.5*aa.m*aa.omega0^2*q^2; return E/aa.omega0, atan(aa.m*aa.omega0*q,p)
end

pendulum_energy(theta::Float64, p::Float64, aa::PendulumActionAngle) = p^2/(2*aa.m*aa.L^2)-aa.m*aa.g*aa.L*cos(theta)

"""Numerical action: J = (1/2pi) contour_integral p dq along energy surface."""
function numerical_action(U::Function, m::Float64, E::Float64, q_range::Tuple{Float64,Float64}; n_points=2000)
    q_min,q_max=q_range; dq=(q_max-q_min)/(n_points-1); J=0.0
    for i in 1:n_points
        q=q_min+(i-1)*dq; Uv=U(q)
        if Uv<E; p=sqrt(2*m*(E-Uv)); w=(i==1||i==n_points)?0.5:1.0; J+=w*p*dq; end
    end
    return J/pi
end

"""Frequency: omega(E)=dE/dJ via finite difference."""
function numerical_frequency(U::Function, m::Float64, E::Float64, q_range::Tuple{Float64,Float64}; dE=1e-3)
    J1=numerical_action(U,m,E-dE/2,q_range); J2=numerical_action(U,m,E+dE/2,q_range)
    return dE/(J2-J1)
end

"""Action for the pendulum: piecewise for libration and rotation."""
function pendulum_action(aa::PendulumActionAngle, E::Float64; n_points=2000)
    sep=aa.m*aa.g*aa.L
    if E<-sep; return 0.0
    elseif E<sep  # libration
        theta_max=acos(-E/(aa.m*aa.g*aa.L))
        U(q)= -aa.m*aa.g*aa.L*cos(q)
        return numerical_action(U,aa.m*aa.L^2,E,(-theta_max,theta_max);n_points=n_points)
    else  # rotation
        U(q)= -aa.m*aa.g*aa.L*cos(q)
        return numerical_action(U,aa.m*aa.L^2,E,(-pi,pi);n_points=n_points)
    end
end

"""Toda lattice action variables (3-particle, classical integrable system)."""
function toda_lattice_hamiltonian(q::Vector{Float64}, p::Vector{Float64}; alpha=1.0)
    return 0.5*sum(p.^2)+exp(q[2]-q[1])+exp(q[3]-q[2])+exp(alpha*(q[1]-q[3]))
end

"""Check if a 1D system is integrable (always true for 1 DOF)."""
is_integrable_1d(H::Function) = true

"""Arnold-Liouville: n-DOF integrable if n involutive integrals exist."""
function arnold_liouville_check(integrals::Vector{Function}, q::Vector{Float64}, p::Vector{Float64}; tol=1e-8)
    m=length(integrals)
    for i in 1:m,j in i+1:m
        if abs(poisson_bracket(integrals[i],integrals[j],q,p))>tol; return false; end
    end
    return true
end

"""Frequency vector omega(J) = grad H0(J) for integrable system."""
function frequency_vector(H0::Function, J::Vector{Float64}; dJ=1e-6)
    n=length(J); omega=zeros(n)
    for i in 1:n; Jp=copy(J);Jp[i]+=dJ;Jm=copy(J);Jm[i]-=dJ; omega[i]=(H0(Jp)-H0(Jm))/(2dJ); end
    return omega
end

"""Resonance condition: k·omega = 0 for integer vector k."""
function is_resonant(omega::Vector{Float64}, k::Vector{Int}; tol=1e-8)
    return abs(dot(Float64.(k),omega))<tol
end

"""无限深势阱的 action: V=0 for |q|<a, V=∞ otherwise. J = (2a/π)√(2mE)."""
function infinite_well_action(m::Float64, a::Float64, E::Float64)
    J = (2*a/pi) * sqrt(2*m*E)
    return J
end

"""无限深势阱的频率: ω(E) = π² ħ²/(2ma²) → classical ω = π√(2E/m)/(2a)."""
function infinite_well_frequency(m::Float64, a::Float64, E::Float64)
    ω = pi * sqrt(2*E/m) / (2*a)
    return ω
end

"""Poschl-Teller 势: V(q) = -V0 sech²(α q). Action: J = √(2mV0)/α - √(-2mE)/α."""
function poschl_teller_action(m::Float64, V0::Float64, alpha::Float64, E::Float64)
    if E >= 0
        return Inf  # 无界运动
    end
    J = (sqrt(2*m*V0) - sqrt(-2*m*E)) / alpha
    return max(J, 0.0)
end

"""开普勒问题的径向作用量: J_r = -√(-2mE) + √(μ m²/(-2E))."""
function kepler_radial_action(m::Float64, mu::Float64, E::Float64, L::Float64)
    if E >= 0
        return Inf  # 无界
    end
    k = mu * m  # GMm
    a = -k / (2*E)
    J_r = -L + k * sqrt(m/(2*abs(E)))
    return max(J_r, 0.0)
end

"""绝热不变量: 在缓慢参数变化下 J 近似守恒。模拟参数 λ(t) 缓慢变化."""
function adiabatic_invariant_simulation(H_func::Function, grad_H::Function, q0::Float64, p0::Float64, lambda_func::Function, t_end::Float64, dt::Float64)
    q, p = q0, p0; J_vals = Float64[]
    n_steps = Int(ceil(t_end/dt))
    for step in 1:n_steps
        t = step*dt; lambda = lambda_func(t)
        H(qq,pp) = H_func(qq,pp,lambda)
        grad(qq,pp) = grad_H(qq,pp,lambda)
        dHdq, dHdp = grad(q, p)
        q += dt*dHdp; p -= dt*dHdq
        if step % 10 == 0
            push!(J_vals, numerical_action(q->0.5*q^2, 1.0, H(q,p), (q-1,q+1)))
        end
    end
    return J_vals
end

"""KAM 定理: 检查无理频率比 (Diophantine 条件)."""
function is_diophantine(omega::Vector{Float64}, k_max::Int=5, gamma::Float64=0.01)
    for k1 in -k_max:k_max, k2 in -k_max:k_max
        k1==0 && k2==0 && continue
        if abs(k1*omega[1] + k2*omega[2]) < gamma/(abs(k1)+abs(k2))^2
            return false
        end
    end
    return true
end
