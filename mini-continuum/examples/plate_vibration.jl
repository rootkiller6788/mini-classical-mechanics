#!/usr/bin/env julia
# plate_vibration.jl — 板振动与屈曲示例
include("../src/Continuum.jl")
using .Continuum

function main()
    println("="^60)
    println("  Plate Vibration & Buckling")
    println("="^60)

    # Steel square plate: a=b=0.5m, h=2mm
    E, nu, rho = 200e9, 0.3, 7800.0
    a, b, h = 0.5, 0.5, 0.002
    D = plate_flexural_rigidity(E, nu, h)
    println("\nSteel plate: $(a*1000)x$(b*1000) mm, h=$(h*1000) mm")
    println("  Flexural rigidity D = $(round(D,1)) Nm")

    # Natural frequencies
    freqs = plate_natural_frequencies(D, rho, h, a, b, 3, 3)
    println("\nNatural frequencies (Hz, simply supported):")
    for (i, w) in enumerate(freqs[1:6])
        println("  mode $i: f = $(round(w/(2*pi),1)) Hz")
    end

    # Boundary condition effects
    println("\nFundamental frequency coefficient lambda^2:")
    for (bc, lam) in BOUNDARY_LAMBDA
        f1 = plate_frequency_coefficient(D, rho, h, a, b, lam)/(2*pi)
        println("  $bc: lambda^2=$(round(lam,2)), f1=$(round(f1,1)) Hz")
    end

    # Circular plate
    radius = 0.3
    D_circ = plate_flexural_rigidity(E, nu, h)
    q = 1000.0  # N/m2
    w_max_ss = circular_plate_max_deflection(D_circ, radius, q, nu)
    w_max_cl = clamped_circular_plate_max(D_circ, radius, q)
    println("\nCircular plate (R=$radius m, q=$q Pa):")
    println("  Simply supported: w_max = $(round(w_max_ss*1000,4)) mm")
    println("  Clamped:         w_max = $(round(w_max_cl*1000,4)) mm")

    # Plate buckling
    P_cr = circular_plate_buckling(D_circ, radius)
    println("  Buckling load: P_cr = $(round(P_cr/1000,1)) kN")

    # Membrane vibration
    T, rho_s = 1000.0, 0.1  # N/m, kg/m2
    c = membrane_wave_speed(T, rho_s)
    mem_freqs = membrane_frequencies(c, a, b, 3, 3)
    println("\nMembrane (T=$T N/m, rho_s=$rho_s):")
    println("  Wave speed c = $(round(c,1)) m/s")
    for (i, w) in enumerate(mem_freqs[1:4])
        println("  mode $i: f = $(round(w/(2*pi),1)) Hz")
    end

    # Shell
    R_shell, L_shell = 0.2, 1.0
    f_shell = cylindrical_shell_frequency(E, nu, rho, R_shell, h, L_shell, 1, 2)/(2*pi)
    println("\nCylindrical shell (R=$R_shell m, L=$L_shell m):")
    println("  Fundamental freq (m=1,n=2): f = $(round(f_shell,1)) Hz")

    println("\nDone.")
end
main()
