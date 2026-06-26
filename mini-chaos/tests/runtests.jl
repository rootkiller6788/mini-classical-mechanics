#!/usr/bin/env julia; include("../src/Chaos.jl"); using .Chaos
p=0;f=0; macro t(n,e) quote try if $(esc(e));global p+=1;println("  ✅ ",$(esc(n)));else;global f+=1;println("  ❌ ",$(esc(n)));end;catch e;global f+=1;println("  ❌ ",$(esc(n)));end;end;end
function main()
    println("="^60); println("  Chaos Tests"); println("="^60)
    x=logistic_map(0.5,2.0); @t("logistic r=2", 0<x<1)
    x=logistic_map(0.5,4.0); @t("logistic r=4", 0<x<1)
    y=lorenz([0,1,0]); @t("Lorenz deriv length", length(y)==3)
    lam=lyapunov_exponent(lorenz,[1,1,1],1.0,0.01); @t("Lorenz lambda>0", lam>0)
    lam_l=logistic_lyapunov(4.0); @t("Logistic r=4 chaotic", lam_l>0.3)
    m=mandelbrot_iter(0,0); @t("Mandelbrot center", m==100)
    bp=find_bifurcation_points(2.8,3); @t("Bifurcation points", length(bp)>=1)
    tot=p+f; println("\n"*"="^60); println("  $p/$tot passed"); println(f==0?"🎉":"⚠️ $f FAILED")
end; main(); exit(f==0?0:1)
