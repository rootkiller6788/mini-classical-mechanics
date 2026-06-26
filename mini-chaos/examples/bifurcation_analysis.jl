#!/usr/bin/env julia
# bifurcation_analysis.jl — Logistic分岔图とFeigenbaum定数の検証
include("../src/Chaos.jl"); using .Chaos

function main()
    println("="^60)
    println("  Bifurcation Analysis — Logistic Map")
    println("="^60)

    # 分岔图
    println("\n--- Logistic Bifurcation Diagram ---")
    r_range = 2.8:0.002:4.0
    rs, xs = logistic_bifurcation(r_range, n_transient=500, n_plot=80)
    println("  Generated $(length(rs)) points across r ∈ [2.8, 4.0]")

    # 找分岔点
    println("\n--- Period Doubling Cascade ---")
    r_bif = find_bifurcation_points(2.8, 5)
    println("  Bifurcation points r_n:")
    for (i, r) in enumerate(r_bif)
        println("    r_$i (period=$(2^i)): $r")
    end

    # Feigenbaum δ
    deltas = estimate_feigenbaum(r_bif)
    println("\n--- Feigenbaum δ Convergence ---")
    for (i, d) in enumerate(deltas)
        println("  δ_$i = $(round(d, digits=6))  (→ $(FEIGENBAUM_DELTA))")
    end

    # Lyapunov指数 vs r
    println("\n--- Lyapunov Exponent λ(r) ---")
    for r in [3.0, 3.5, 3.57, 3.83, 3.9, 4.0]
        lam = logistic_lyapunov(r; n_iter=5000)
        status = lam > 0.01 ? "CHAOTIC" : lam < -0.01 ? "PERIODIC" : "BIFURCATION"
        println("  r=$r: λ=$(round(lam,5)) → $status")
    end

    # 周期検出
    println("\n--- Period Detection ---")
    for r in [3.2, 3.5, 3.83]
        per = find_period_at_r(r, logistic_map)
        println("  r=$r: period = $per")
    end

    println("\n✅ bifurcation_analysis.jl done.")
end
main()
