#!/usr/bin/env julia
include("../src/VariationalPrinciple.jl"); using .VariationalPrinciple

function main()
    println("="^60); println("  Brachistochrone — Variational Solution"); println("="^60)
    g=9.81
    # F = sqrt(1+y'^2)/sqrt(2gy), independent of x -> Beltrami identity
    F(x,y,yp) = sqrt(1+yp^2)/sqrt(2*g*abs(y)+1e-15)
    dF_dyp(x,y,yp) = yp/(sqrt(1+yp^2)*sqrt(2*g*abs(y)+1e-15))

    # At the solution (cycloid), Beltrami constant should be uniform
    x1,y1=1.0,0.5
    R,th1=2.0,1.5  # approximate cycloid params
    cycloid_y(x) = 0.2*(1-cos(0.8*x))  # rough approximation for demo
    cycloid_yp(x) = 0.16*sin(0.8*x)

    println("Beltrami constant at x=0.2: $(round(beltrami_constant(F,dF_dyp,cycloid_y,cycloid_yp,0.2),digits=6))")
    println("Beltrami constant at x=0.8: $(round(beltrami_constant(F,dF_dyp,cycloid_y,cycloid_yp,0.8),digits=6))")
    println("✅ brachistochrone.jl done.")
end; main()
