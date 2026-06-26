#!/usr/bin/env julia
include("../src/CelestialMechanics.jl"); using .CelestialMechanics
function main()
println("="^60); println("  Hohmann Transfer — Earth→Mars"); println("="^60)

r_earth=1.0; r_mars=1.524; mu=G_SUN
hoh=hohmann_transfer(r_earth,r_mars;mu=mu)

println("\nEarth orbit: r=1.0 AU"); println("Mars orbit: r=1.524 AU")
println("\n--- Transfer Orbit ---"); println("  a_trans = $(round(hoh.a_trans,4)) AU")
println("  e_trans = $(round(hoh.e_trans,4))"); println("  Transfer time = $(round(hoh.transfer_time,4)) yr ($(round(hoh.transfer_time*365.25,1)) days)")
println("\n--- Delta-V Budget ---"); println("  Δv₁ (depart) = $(round(hoh.delta_v1,4)) AU/yr"); println("  Δv₂ (arrive) = $(round(hoh.delta_v2,4)) AU/yr")
println("  Δv_total     = $(round(hoh.delta_v_total,4)) AU/yr")
println("  In km/s: Δv₁≈$(round(hoh.delta_v1*4.74047,3)), Δv₂≈$(round(hoh.delta_v2*4.74047,3)), total≈$(round(hoh.delta_v_total*4.74047,3))")

phi=phasing_angle_hohmann(r_earth,r_mars;mu=mu)
println("\nLaunch window: phase angle φ = $(round(phi*180/pi,2))°")
println("\n✅ hohmann_transfer.jl done.")
end; main()
