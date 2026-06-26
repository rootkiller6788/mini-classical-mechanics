#!/usr/bin/env julia
include("../src/VariationalPrinciple.jl"); using .VariationalPrinciple
p=0;f=0
macro t(n,e) quote try if $(esc(e));global p+=1;println("  ✅ ",$(esc(n)));else;global f+=1;println("  ❌ ",$(esc(n)));end;catch e;global f+=1;println("  ❌ ",$(esc(n)));end;end;end
function main()
    println("="^60); println("  Variational Principle Tests"); println("="^60)
    F(x,y,yp)=yp^2+y^2; dFdy(x,y,yp)=2y; dFdyp(x,y,yp)=2yp
    J=evaluate_functional(F,x->sin(x),x->cos(x),0,pi;n=500); @t("functional >0", J>0)
    bc=beltrami_constant(F,dFdyp,x->1.0,x->0.0,0.5); @t("Beltrami const", abs(bc-1)<1e-8)
    gd=gateaux_derivative(F,x->sin,x->cos,x->cos,x->-sin(x),0,pi); @t("Gateaux deriv finite", isfinite(gd))
    xs,fd=functional_derivative(F,dFdy,dFdyp,x->sin,x->cos,0,pi;n=100); @t("functional deriv length", length(fd)==100)
    lc=legendre_condition((x,y,yp)->2.0,0,0,0); @t("Legendre F''=2>0", lc)
    w=weierstrass_excess(F,dFdyp,0,1,2,3); @t("Weierstrass finite", isfinite(w))
    nc=natural_boundary_condition(F,dFdyp,0,1,2); @t("Natural BC", abs(nc-4)<1e-8)
    lam=isoperimetric_lagrange_multiplier(F,F,0.5,0,1,0,0); @t("Isoperimetric lambda", isfinite(lam))
    tot=p+f; println("\n"*"="^60); println("  $p/$tot passed"); println(f==0?"🎉":"⚠️ $f FAILED")
end; main(); exit(f==0?0:1)
