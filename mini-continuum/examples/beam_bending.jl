#!/usr/bin/env julia
# beam_bending.jl — Euler-Bernoulli 与 Timoshenko 梁弯曲示例
include("../src/Continuum.jl")
using .Continuum

function main()
    println("="^60)
    println("  Beam Bending — Euler-Bernoulli vs Timoshenko")
    println("="^60)

    # Steel rectangular beam: b=20mm, h=40mm, L=1m
    E, nu, rho = 200e9, 0.3, 7800.0
    b, h, L = 0.02, 0.04, 1.0
    I = b*h^3/12; A = b*h
    EI = E*I
    G = E/(2*(1+nu))
    kappa = 5/6  # 矩形截面剪切系数
    GAKs = kappa*G*A

    println("\nBeam: $(b*1000)x$(h*1000) mm, L=$L m, Steel")
    println("  EI = $(round(EI/1000,1)) Nm2, GAKs = $(round(GAKs/1e6,1)) MN")

    # 均布载荷
    q = 1000.0  # N/m
    w_eb = beam_max_deflection_uniform(EI, L, q)
    println("\nUniform load q=$q N/m:")
    println("  Euler-Bernoulli: w_max = $(round(w_eb*1000,4)) mm")

    # 集中力
    P = 500.0  # N at midspan
    w_eb_point = beam_max_deflection_point(EI, L, P)
    w_tot, w_bend, w_shear = timoshenko_deflection_point(EI, L, P, GAKs)
    println("\nPoint load P=$P N at midspan:")
    println("  Euler-Bernoulli:  w_max = $(round(w_eb_point*1000,4)) mm")
    println("  Timoshenko total:  w_max = $(round(w_tot*1000,4)) mm")
    println("    bending: $(round(w_bend*1000,4)) mm, shear: $(round(w_shear*1000,4)) mm")
    alpha = shear_deformation_ratio(EI, GAKs, L)
    println("  Shear ratio alpha = $(round(alpha,4))")

    # FDM求解
    xs, w_fdm = euler_bernoulli_deflection(EI, L, x -> q)
    println("\nFDM solution ($(length(xs)) nodes):")
    println("  w_max(FDM) = $(round(maximum(abs.(w_fdm))*1000,4)) mm")

    # 悬臂梁
    w_cant_tip = cantilever_deflection(EI, L, P, L)
    w_cant_uniform = cantilever_uniform(EI, L, q, L)
    println("\nCantilever (tip load P=$P N):")
    println("  w_tip = $(round(w_cant_tip*1000,4)) mm")
    println("Cantilever (uniform q=$q N/m):")
    println("  w_tip = $(round(w_cant_uniform*1000,4)) mm")

    # 梁振动
    rhoA = rho*A
    omegas = beam_natural_frequencies(EI, rhoA, L, 5)
    println("\nNatural frequencies (Hz):")
    for (i, w) in enumerate(omegas)
        println("  mode $i: f = $(round(w/(2*pi),1)) Hz")
    end
    f_cant = cantilever_fundamental_frequency(EI, rhoA, L)/(2*pi)
    println("  cantilever fundamental: f = $(round(f_cant,1)) Hz")

    println("\nDone.")
end
main()
