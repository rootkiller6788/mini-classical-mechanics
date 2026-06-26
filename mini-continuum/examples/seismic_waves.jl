#!/usr/bin/env julia
include("../src/Continuum.jl"); using .Continuum

function main()
    println("="^60)
    println("  Seismic Waves — Elastic Wave Speeds & Snell")
    println("="^60)

    # Crust properties
    E, nu, rho = 70e9, 0.25, 2700.0
    lam, mu = lame_constants(E, nu)
    vp = p_wave_speed(lam, mu, rho)
    vs = s_wave_speed(mu, rho)
    vr = rayleigh_wave_speed(lam, mu, rho)

    println("\nGranite (E=$E Pa, ν=$nu, ρ=$rho kg/m³):")
    println("  P-wave: vp = $(round(vp)) m/s")
    println("  S-wave: vs = $(round(vs)) m/s")
    println("  Rayleigh: vr = $(round(vr)) m/s")
    println("  vp/vs = $(round(vp/vs, 2))")

    # Snell's law: wave propagating from crust to mantle
    rho_mantle, vp_mantle = 3300.0, 8100.0
    Z1 = acoustic_impedance(rho, vp)
    Z2 = acoustic_impedance(rho_mantle, vp_mantle)
    R = reflection_coefficient(Z1, Z2)
    T = transmission_coefficient(Z1, Z2)

    println("\nCrust→Mantle interface:")
    println("  Z_crust = $(round(Z1/1e6,2)) MRayl")
    println("  Z_mantle = $(round(Z2/1e6,2)) MRayl")
    println("  Reflection coeff: $(round(R,4))")
    println("  Transmission coeff: $(round(T,4))")

    # Critical angle
    theta_c = critical_angle(vp, vp_mantle)
    println("  Critical angle: $(round(rad2deg(theta_c),1))°")

    # 1D wave demo
    println("\n--- 1D Wave Equation (FDTD) ---")
    c = 10.0; L = 1.0; T = 0.5
    u0(x) = sin(2π*x/L)
    v0(x) = 0.0
    xs, ts, u = wave_1d_solve(c, L, T, u0, v0; nx=50, CFL=0.4)
    println("  nx=50, nt=$(size(u,2)), Δx=$(round(xs[2]-xs[1],4)), Δt=$(round(ts[2]-ts[1],6))")
    println("  u(L/2, T) ≈ $(round(u[25,end],4))")

    println("\n✅ seismic_waves.jl done.")
end
main()
