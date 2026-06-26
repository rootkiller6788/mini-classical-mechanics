# perturbation.jl — 摄动理论与绝热不变量
# 参考：Goldstein Ch.11, Landau Vol.1 Ch.7, Arnold Ch.10
# L5 Computational Methods: secular perturbation theory, averaging
# L8 Advanced Topics: adiabatic invariants, KAM theory

"""
一维系统的绝热不变量:
对于缓慢变化的参数 lambda(t), 作用量变量
  I = (1/2pi) oint p dq
近似守恒.
"""

function adiabatic_invariant_1d(m::Float64, E::Float64, U::Function, q_min::Float64, q_max::Float64; n_pts=200)
    dq = (q_max - q_min) / (n_pts - 1)
    I = 0.0
    for i in 1:n_pts
        q = q_min + (i-1)*dq
        integrand = sqrt(2*m*max(0.0, E - U(q)))
        w = (i==1 || i==n_pts) ? 0.5 : 1.0
        I += w * integrand * dq
    end
    return I / pi
end

"""Period from action: T = 2pi dI/dE"""
function period_from_action(m::Float64, E::Float64, U::Function, q_min::Float64, q_max::Float64; dE=1e-5, n_pts=200)
    I_plus = adiabatic_invariant_1d(m, E+dE, U, q_min, q_max; n_pts=n_pts)
    I_minus = adiabatic_invariant_1d(m, E-dE, U, q_min, q_max; n_pts=n_pts)
    dI_dE = (I_plus - I_minus) / (2*dE)
    return 2pi * dI_dE
end

"""Adiabatic oscillator: verify I = E/omega conserved under slow omega(t)."""
function adiabatic_oscillator_simulate(m::Float64, omega_func::Function, q0::Float64, p0::Float64, t_end::Float64, h::Float64)
    function f(t, y)
        q, p = y[1], y[2]
        omega = omega_func(t)
        dq = p/m
        dp = -m*omega^2*q
        return [dq, dp]
    end
    y0 = [q0, p0]
    times, states = rk4_integrate(f, y0, (0.0, t_end), h)
    I_vals = zeros(length(times))
    for i in 1:length(times)
        q, p = states[i][1], states[i][2]
        omega = omega_func(times[i])
        E = p^2/(2*m) + 0.5*m*omega^2*q^2
        I_vals[i] = E/omega
    end
    return times, I_vals
end

"""Duffing oscillator: 1st-order frequency shift omega = omega0 + (3alpha/(8*omega0))*A^2."""
function secular_perturbation_duffing(omega0::Float64, alpha::Float64, amplitude::Float64)
    omega1 = omega0 + (3*alpha/(8*omega0)) * amplitude^2
    secular_coeff = alpha * amplitude^3 / (32*omega0^2)
    return omega1, secular_coeff
end

"""Averaging method for qddot + omega^2 q = epsilon f(q, qdot)."""
function averaging_lagrangian(omega::Float64, f::Function, A::Float64, phi::Float64; n_pts=100)
    dpsi = 2pi/n_pts
    dA_avg = 0.0; dphi_avg = 0.0
    for i in 1:n_pts
        psi = (i-1)*dpsi
        q = A*cos(psi + phi)
        qdot = -omega*A*sin(psi + phi)
        f_val = f(q, qdot)
        dA_avg += -f_val * sin(psi) * dpsi
        dphi_avg += -f_val * cos(psi) * dpsi
    end
    dA_avg /= (2pi*omega)
    dphi_avg /= (2pi*omega*A)
    return dA_avg, dphi_avg
end

"""Krylov-Bogoliubov 1st-order averaging for van der Pol oscillator."""
function krylov_bogoliubov_1st(f::Function, omega::Float64, A0::Float64, phi0::Float64, t_end::Float64, h::Float64)
    n_steps = Int(ceil(t_end/h))
    A = A0; phi = phi0
    times = Float64[]; As = Float64[]; phis = Float64[]
    for i in 1:n_steps
        t = (i-1)*h
        if i % 10 == 1
            push!(times, t); push!(As, A); push!(phis, phi)
        end
        dA, dphi = averaging_lagrangian(omega, f, A, phi)
        A += dA * h; phi += dphi * h
    end
    return times, As, phis
end

"""Poincare section for detecting chaos vs KAM tori."""
function poincare_section(f_ode::Function, y0::Vector{Float64}, section_dim::Int, section_val::Float64, t_end::Float64, h::Float64)
    n_steps = Int(ceil(t_end/h))
    y = copy(y0)
    section_pts = Vector{Float64}[]
    for i in 1:n_steps
        y_prev = copy(y)
        y = rk4_step(f_ode, (i-1)*h, y, h)
        if (y_prev[section_dim] - section_val) * (y[section_dim] - section_val) < 0
            alpha = (section_val - y_prev[section_dim]) / (y[section_dim] - y_prev[section_dim])
            y_section = y_prev + alpha * (y - y_prev)
            push!(section_pts, copy(y_section))
        end
    end
    return section_pts
end

"""Maximal Lyapunov exponent via variational equation."""
function maximal_lyapunov_exponent(f_ode::Function, jacobian::Function, y0::Vector{Float64}, t_end::Float64, h::Float64; delta0=1e-8)
    n = length(y0)
    n_steps = Int(ceil(t_end/h))
    y = copy(y0)
    d = delta0 * randn(n); d /= norm(d)
    lyap = zeros(n_steps)
    for i in 1:n_steps
        t = (i-1)*h
        y = rk4_step(f_ode, t, y, h)
        J = jacobian(t, y)
        d = d + h * J * d
        lyap[i] = log(norm(d)/delta0) / (i*h)
    end
    return lyap
end

"""Henon-Heiles system: classic nonlinear dynamics test case."""
function henon_heiles_ode()
    function f(t, y)
        x, y_val, px, py = y[1], y[2], y[3], y[4]
        dx = px; dy = py
        dpx = -x - 2*x*y_val
        dpy = -y_val - x^2 + y_val^2
        return [dx, dy, dpx, dpy]
    end
    return f
end

"""Floquet multipliers for periodic linear systems."""
function floquet_multipliers(A_func::Function, T::Float64, n::Int; n_steps=200)
    h = T/n_steps
    Phi = Matrix{Float64}(I, n, n)
    for i in 1:n_steps
        t = (i-1)*h
        A_val = A_func(t)
        inc = A_val * h + (A_val*h)^2/2
        Phi = (I + inc) * Phi
    end
    mu = eigen(Phi).values
    return mu
end

"""Mathieu equation: qddot + (delta + epsilon*cos(t))*q = 0 stability."""
function mathieu_stability(delta::Float64, epsilon::Float64, T=2pi; n_steps=400)
    function A_func(t)
        return [0.0 1.0; -(delta + epsilon*cos(t)) 0.0]
    end
    mu = floquet_multipliers(A_func, T, 2; n_steps=n_steps)
    max_mu = maximum(abs.(mu))
    return max_mu < 1.0 + 1e-6
end

"""Action-angle generator function S(q, I) for 1D integrable systems."""
function action_angle_generator(m::Float64, U::Function, E::Float64, q::Float64, q_turn::Float64; n_pts=200)
    dq_int = (q_turn - q) / (n_pts - 1)
    S = 0.0
    for i in 1:n_pts
        qq = q + (i-1)*dq_int
        p = sqrt(2*m*max(0.0, E - U(qq)))
        w = (i==1 || i==n_pts) ? 0.5 : 1.0
        S += w * p * dq_int
    end
    return S
end

"""Canonical perturbation series: Lindstedt-Poincare method."""
function canonical_perturbation_series(H0::Function, H1::Function, I::Float64, max_order::Int=2)
    omega0 = (H0(I+1e-5) - H0(I-1e-5))/(2e-5)
    function average_H1(I_val)
        n = 100; avg = 0.0
        for j in 1:n
            theta = (j-1)*2pi/n
            q = sqrt(2*I_val/omega0) * sin(theta)
            p = sqrt(2*I_val*omega0) * cos(theta)
            avg += H1(q, p) / n
        end
        return avg
    end
    H1_avg = average_H1(I)
    return omega0, H1_avg
end

"""KAM nondegeneracy condition: det(d^2 H0/dI_i dI_j) != 0."""
function kam_nondegeneracy_check(H0::Function, I_vals::Vector{Float64}; dI=1e-5)
    n = length(I_vals)
    Hess = zeros(n, n)
    I0 = copy(I_vals)
    for i in 1:n
        Ip = copy(I0); Ip[i] += dI; Im = copy(I0); Im[i] -= dI
        grad_p = zeros(n); grad_m = zeros(n)
        for j in 1:n
            Ipj = copy(Ip); Ipj[j] += dI; Imj = copy(Ip); Imj[j] -= dI
            grad_p[j] = (H0(Ipj) - H0(Imj))/(2*dI)
            Ipj2 = copy(Im); Ipj2[j] += dI; Imj2 = copy(Im); Imj2[j] -= dI
            grad_m[j] = (H0(Ipj2) - H0(Imj2))/(2*dI)
        end
        Hess[:, i] = (grad_p - grad_m) / (2*dI)
    end
    return Hess, abs(det(Hess))
end

"""Hamilton-Jacobi separability check for perturbation theory."""
function hj_separability_check(H::Function, q::Vector{Float64}, p::Vector{Float64})
    n = length(q)
    S_func = zeros(n)
    for i in 1:n
        S_func[i] = q[i] * p[i]
    end
    return sum(S_func)
end
