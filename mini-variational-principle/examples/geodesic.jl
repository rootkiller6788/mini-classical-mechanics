#!/usr/bin/env julia
include("../src/VariationalPrinciple.jl"); using .VariationalPrinciple

function main()
    println("="^60)
    println("  Geodesics — Shortest Paths via Variational Principle")
    println("="^60)

    # Plane geodesic (straight line)
    println("\n--- Plane Geodesic ---")
    F_plane(x, y, yp) = sqrt(1 + yp^2)  # arc length
    dFdy_plane(x, y, yp) = 0.0
    dFdyp_plane(x, y, yp) = yp / sqrt(1 + yp^2 + 1e-15)

    # Straight line: y(x) = (yb-ya)/(b-a)*(x-a) + ya, slope constant
    y_plane(x) = x; yp_plane(x) = 1.0
    S_plane = evaluate_functional(F_plane, y_plane, yp_plane, 0.0, 1.0)
    println("  [0,0]→[1,1] length = $(round(S_plane,4))  (exact: √2 = $(round(sqrt(2),4)))")

    # Sphere geodesic (great circle via Beltrami)
    println("\n--- Sphere Geodesic (Great Circle) ---")
    println("  Functional: ∫ √(R²θ̇² + R²sin²θ·φ̇²) dt")
    println("  Solution: great circles")
    println("  ═ Equator: constant θ=π/2, φ varies")
    println("  ║ Meridian: constant φ, θ varies")
    R = 6371.0  # Earth radius
    dist_equator = pi * R / 2  # quarter of equator = pole→equator
    println("  Pole→Equator: $(round(dist_equator)) km (¼ circumference)")

    # Brachistochrone check
    println("\n--- Brachistochrone (Quick Check) ---")
    g = 9.81
    F_brach(x, y, yp) = sqrt((1+yp^2)/(2*g*y+1e-15))
    # Cycloid: x=a(θ-sinθ), y=a(1-cosθ)
    println("  Solution: cycloid (not solved here, see brachistochrone.jl)")

    # Catenary check
    println("\n--- Catenary ---")
    F_cat(x, y, yp) = y * sqrt(1 + yp^2)
    # cosh solution: y = a·cosh(x/a)
    a_cat = 1.0
    y_cat(x) = a_cat * cosh(x/a_cat)
    yp_cat(x) = sinh(x/a_cat)
    S_cat = evaluate_functional(F_cat, y_cat, yp_cat, -1.0, 1.0)
    println("  y=cosh(x), x∈[-1,1]: surface area = $(round(S_cat,4))")

    println("\n✅ geodesic.jl done.")
end
main()
