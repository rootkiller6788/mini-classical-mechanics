#!/usr/bin/env julia; include("../src/Chaos.jl"); using .Chaos
function main()
    println("="^60); println("  Chaos Suite"); println("="^60)
    # Lorenz
    x=[1,1,1]; dt=0.005
    for _ in 1:1000; x+=dt*lorenz(x); end
    println("Lorenz (5s): x=$(round.(x,digits=3))")
    lam=lyapunov_exponent(lorenz,[1,1,1],10.0,0.005)
    println("Lorenz lambda=$(round(lam,digits=4)) ($(lam>0?"CHAOTIC":"regular"))")
    # Logistic
    println("\nLogistic map Lyapunov exponents:")
    for r in [2.0,3.2,3.5,3.57,4.0]
        lam_l=logistic_lyapunov(r); println("  r=$r: lambda=$(round(lam_l,digits=4)) ($(lam_l>0.01?"chaotic":lam_l< -0.01?"periodic":"edge"))")
    end
    # Bifurcation
    bp=find_bifurcation_points(2.8,4); println("\nBifurcation points: $(round.(bp,digits=3))")
    deltas=estimate_feigenbaum(bp); println("Feigenbaum delta estimates: $(round.(deltas,digits=3)) (theory: 4.669)")
    println("✅ chaos done.")
end; main()
