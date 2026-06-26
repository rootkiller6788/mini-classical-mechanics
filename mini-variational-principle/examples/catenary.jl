#!/usr/bin/env julia
include("../src/VariationalPrinciple.jl"); using .VariationalPrinciple

function main()
    println("="^60); println("  Catenary — Hanging Chain"); println("="^60)
    # F = y*sqrt(1+y'^2), independent of x
    F(x,y,yp) = y*sqrt(1+yp^2)
    dF_dyp(x,y,yp) = y*yp/sqrt(1+yp^2)

    cat_y(x) = cosh(x); cat_yp(x) = sinh(x)
    println("Catenary y=cosh(x):")
    for x in [0.0,0.5,1.0,1.5,2.0]
        bc = beltrami_constant(F,dF_dyp,cat_y,cat_yp,x)
        println("  x=$x: Beltrami=$(round(bc,digits=6)) (should be const)")
    end
    println("✅ catenary.jl done.")
end; main()
