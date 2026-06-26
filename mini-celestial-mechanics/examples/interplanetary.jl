#!/usr/bin/env julia
include("../src/CelestialMechanics.jl"); using .CelestialMechanics

function main()
    println("="^60)
    println("  Interplanetary Transfer — Patched Conics")
    println("="^60)

    mu_sun = G_SUN
    r_earth = 1.0; r_mars = 1.524; r_jupiter = 5.203

    println("\n--- Earth→Mars Hohmann ---")
    hoh = hohmann_transfer(r_earth, r_mars; mu=mu_sun)
    println("  Transfer a = $(round(hoh.a_trans,4)) AU")
    println("  Transfer e = $(round(hoh.e_trans,4))")
    println("  Δv total = $(round(hoh.delta_v_total,4)) AU/yr (~$(round(hoh.delta_v_total*4.74,2)) km/s)")
    println("  Time = $(round(hoh.transfer_time,3)) yr ($(round(hoh.transfer_time*365.25)) days)")

    println("\n--- Earth→Jupiter Hohmann ---")
    hoh_j = hohmann_transfer(r_earth, r_jupiter; mu=mu_sun)
    println("  Δv total = $(round(hoh_j.delta_v_total*4.74,2)) km/s")
    println("  Time = $(round(hoh_j.transfer_time,2)) yr")

    # Gravity assist at Jupiter
    println("\n--- Jupiter Gravity Assist ---")
    v_jupiter = sqrt(mu_sun / r_jupiter)  # orbital velocity
    v_inf = 5.64  # km/s typical for Jupiter flyby
    delta = gravity_assist_turn_angle(v_inf/4.74, 71492.0/149597870.7; mu=G_SUN*9.5479e-4)  # Jupiter in solar units
    dv_max = gravity_assist_delta_v(v_jupiter, v_inf/4.74, 71492.0/149597870.7; mu=G_SUN*9.5479e-4)
    println("  Jupiter orbital v = $(round(v_jupiter*4.74,1)) km/s")
    println("  Turn angle δ = $(round(rad2deg(delta),1))°")
    println("  Max Δv gain = $(round(dv_max*4.74,2)) km/s")

    # Phasing
    println("\n--- Phasing Windows ---")
    phi_em = phasing_angle_hohmann(r_earth, r_mars; mu=mu_sun)
    println("  Earth→Mars: φ = $(round(rad2deg(phi_em),1))°")
    T_syn = synodic_period_from_radii(r_earth, r_mars; mu=mu_sun)
    println("  Synodic period: $(round(T_syn,2)) yr ($(round(T_syn*365.25)) days)")

    println("\n✅ interplanetary.jl done.")
end
main()
