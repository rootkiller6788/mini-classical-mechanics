#!/usr/bin/env julia
include("../src/Continuum.jl"); using .Continuum

function main()
    println("="^60)
    println("  Euler Buckling & Column Stability")
    println("="^60)

    E, L = 200e9, 2.0  # steel, 2m
    b, h = 0.02, 0.04   # rectangular cross section
    I = b*h^3/12
    A = b*h
    r = radius_of_gyration(I, A)
    slenderness = slenderness_ratio(1.0, L, I, A)

    println("\nSteel column: E=$E Pa, L=$L m, $(b*1000)×$(h*1000) mm")
    println("  I = $(round(I*1e8,4))e-8 m⁴")
    println("  r = $(round(r*1000,2)) mm")
    println("  λ = L/r = $(round(slenderness,1))")

    for (name, K) in BUCKLING_K
        Pcr = euler_buckling_load(E*I, L, 1, K)
        sigma_cr = Pcr / A
        println("  $name (K=$K): Pcr = $(round(Pcr/1000,2)) kN, σcr = $(round(sigma_cr/1e6,1)) MPa")
    end

    println("\n--- Critical stress vs slenderness ---")
    for lam in [50.0, 100.0, 150.0, 200.0]
        sigma_cr = critical_stress(E, lam)
        println("  λ=$lam: σcr = $(round(sigma_cr/1e6,1)) MPa")
    end

    # Euler-Bernoulli beam
    println("\n--- Simply Supported Beam (uniform load) ---")
    EI = E*I; Lb = 1.0; q = 1000.0
    w_max = beam_max_deflection_uniform(EI, Lb, q)
    println("  q=$q N/m, w_max = $(round(w_max*1000,3)) mm")

    # Cantilever
    w_cant = cantilever_deflection(EI, Lb, 500.0, Lb)
    println("  Cantilever, P=500N at tip: w = $(round(w_cant*1000,3)) mm")

    println("\n✅ buckling.jl done.")
end
main()
