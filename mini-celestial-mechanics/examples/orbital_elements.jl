#!/usr/bin/env julia; include("../src/CelestialMechanics.jl"); using .CelestialMechanics
function main()
    println("="^60); println("  Celestial Mechanics"); println("="^60)
    a=1.0; e=0.0167; i=deg2rad(0.0); Omega=0.0; omega=deg2rad(102.9); M=deg2rad(100.0)
    r,v=orbital_elements_to_state(a,e,i,Omega,omega,M)
    println("Elements->State: r=$(round.(r,digits=3)), v=$(round.(v,digits=3))")
    els=state_to_orbital_elements(r,v)
    println("State->Elements: a=$(round(els.a,digits=4)), e=$(round(els.e,digits=4))")
    pts=lagrange_points(3.0e-6)  # Earth-Moon mu
    println("Earth-Moon L1 x=$(round(pts.L1[1],digits=4)) AU")
    soi=sphere_of_influence(1.0,3e-6,1.0); println("Earth SOI: $(round(soi,digits=6)) AU")
    h=hill_sphere(1.0,3e-6,1.0); println("Earth Hill: $(round(h,digits=6)) AU")
    println("✅ celestial done.")
end; main()
