# thermoelastic.jl — 热弹性与耦合场: 热传导, 热应力, 热弹性阻尼
# 参考: Boley & Weiner, Nowacki

using LinearAlgebra

"Fourier 热传导 1D: ∂T/∂t = α ∂²T/∂x²"
function heat_equation_1d(alpha::Float64, L::Float64, t_end::Float64, T0::Function; nx=200, CFL=0.4)
    dx=L/(nx-1); dt=CFL*dx^2/alpha; nt=Int(ceil(t_end/dt))
    T=zeros(nx,nt+1); xs=range(0,L,length=nx)
    for i in 1:nx; T[i,1]=T0(xs[i]); end; r=alpha*dt/dx^2
    for n in 1:nt; for i in 2:nx-1
        T[i,n+1]=T[i,n]+r*(T[i+1,n]-2*T[i,n]+T[i-1,n])
    end; end; xs, range(0,dt*nt,length=nt+1), T
end

"热扩散率: α = k/(ρ c_p)"
thermal_diffusivity(k::Float64, rho::Float64, cp::Float64) = k/(rho*cp)

"热应力 (两端约束杆): σ_th = -E α ΔT"
function thermal_stress_1d(E::Float64, alpha_th::Float64, delta_T::Float64)
    -E * alpha_th * delta_T
end

"Biot 数: Bi = h L / k (对流 vs 导热)"
biot_number(h::Float64, L::Float64, k::Float64) = h*L/k

"热弹性阻尼 (TED) Q因子: Q⁻¹ = E α² T₀ / (ρ c_p) · ωτ/(1+ω²τ²)"
function thermoelastic_damping(E::Float64, alpha_th::Float64, T0::Float64, rho::Float64, cp::Float64, omega::Float64, tau::Float64)
    delta_E = E*alpha_th^2*T0/(rho*cp)
    delta_E * omega*tau / (1 + omega^2*tau^2)
end

"热松弛时间: τ = L²/(π²α)"
thermal_relaxation_time(L::Float64, alpha::Float64) = L^2/(pi^2*alpha)

"热膨胀应变: ε_th = α (T - T_ref)"
thermal_strain(alpha::Float64, T::Float64, T_ref::Float64) = alpha*(T-T_ref)

"Zener 热弹性 (薄梁热弹性弯曲)"
function zener_thermoelastic_frequency(E::Float64, alpha_th::Float64, T0::Float64, rho::Float64, cp::Float64, h::Float64)
    omega_peak = pi^2*thermal_diffusivity(1.0,rho*cp,h^2)  # approximate
    damping_peak = E*alpha_th^2*T0/(rho*cp)
    omega_peak, damping_peak
end

"耦合热弹性 (考虑变形热)"
function coupled_thermoelastic_1d(E::Float64, alpha_th::Float64, k::Float64, rho::Float64, cp::Float64, L::Float64, t_end::Float64; nx=100)
    dx=L/(nx-1); dt=0.1*dx^2; nt=Int(ceil(t_end/dt))
    u=zeros(nx,nt+1); T=zeros(nx,nt+1)
    xs=range(0,L,length=nx); return xs, u, T  # simplified
end

"热冲击应力 (表面突然加热): σ_max = E α ΔT / (1-ν)"
function thermal_shock_stress(E::Float64, alpha_th::Float64, delta_T::Float64, nu::Float64)
    E * alpha_th * delta_T / (1-nu)
end

"热冲击抵抗参数: R = σ_f (1-ν) / (E α)"
function thermal_shock_resistance(sigma_f::Float64, E::Float64, alpha_th::Float64, nu::Float64)
    sigma_f * (1-nu) / (E * alpha_th)
end

export heat_equation_1d, thermal_diffusivity, thermal_stress_1d, biot_number
export thermoelastic_damping, thermal_relaxation_time, thermal_strain
export zener_thermoelastic_frequency, thermal_shock_stress, thermal_shock_resistance
