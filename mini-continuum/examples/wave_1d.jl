#!/usr/bin/env julia; include("../src/Continuum.jl"); using .Continuum
function main()
    println("="^60); println("  Continuum — Elastic Waves & Beam"); println("="^60)
    L=10.0; T=5.0; c=1.0
    u0(x)=exp(-(x-5)^2); v0(x)=0.0
    x,t,u=wave_1d_solve(c,L,T,u0,v0;nx=200)
    println("Wave grid: $(length(x))x$(length(t)), max amplitude: $(round(maximum(abs.(u)),digits=4))")
    lam,mu=lame_constants(200e9,0.3)
    println("Steel: lam=$(round(lam/1e9,digits=1))GPa, mu=$(round(mu/1e9,digits=1))GPa, E=200GPa, nu=0.3")
    vp=p_wave_speed(lam,mu,7800); vs=s_wave_speed(mu,7800)
    println("P-wave: $(round(vp))m/s, S-wave: $(round(vs))m/s")
    # Beam
    xs,w=euler_bernoulli_deflection(1e6,10,x->1000.0)
    println("Beam (L=10, EI=1e6, q=1000): max deflection ~$(round(maximum(abs.(w)),digits=4))m")
    modes=beam_natural_frequencies(1e6,10,100); println("Beam modes (first 3): $(round.(modes[1:3],digits=2)) rad/s")
    println("✅ continuum_demo done.")
end; main()
