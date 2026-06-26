#!/usr/bin/env julia
include("../src/Chaos.jl"); using .Chaos

function main()
    println("="^60)
    println("  Lyapunov Spectrum — Lorenz Attractor")
    println("="^60)

    # Lorenz maximum Lyapunov via two-particle method
    x0 = [1.0, 1.0, 1.0]
    lam1 = lyapunov_exponent(s -> lorenz(s), x0, 50.0, 0.005; d0=1e-6)
    println("\nLorenz (σ=10, ρ=28, β=8/3):")
    println("  λ₁ = $(round(lam1, 5)) (positive → chaotic)")

    # Logistic map λ(r)
    println("\nLogistic λ(r) scan:")
    for r in [2.5, 3.0, 3.2, 3.5, 3.57, 3.7, 3.83, 3.9, 4.0]
        lam = logistic_lyapunov(r)
        sig = lam > 0.01 ? "+" : lam < -0.01 ? "-" : "0"
        println("  λ($r) = $(round(lam,5)) $sig")
    end

    # Hénon map Lyapunov pair
    println("\nHénon map (a=1.4, b=0.3):")
    l1, l2 = henon_lyapunovs()
    println("  λ₁ = $(round(l1,5))")
    println("  λ₂ = $(round(l2,5))")
    println("  λ₁+λ₂ = $(round(l1+l2,5)) ≈ log|b| = $(round(log(0.3),5))")

    # Rössler λ
    println("\nRössler (a=0.2, b=0.2, c=5.7):")
    lam_ros = lyapunov_exponent(s -> rossler(s), [1.0, 1.0, 1.0], 50.0, 0.005)
    println("  λ₁ = $(round(lam_ros, 5))")

    println("\n✅ lyapunov_spectrum.jl done.")
end
main()
