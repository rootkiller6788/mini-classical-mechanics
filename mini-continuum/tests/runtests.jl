#!/usr/bin/env julia
include("../src/Continuum.jl"); using .Continuum; using LinearAlgebra
p=0;f=0
macro t(n,ex) quote try if ex;global p+=1;println("  * ",n);else;global f+=1;println("  X ",n);end;catch e;global f+=1;println("  X ",n,": ",e);end;end;end
function main()
    println("="^60); println("  Continuum Tests"); println("="^60)
    eps_m=strain_tensor([0.01 0 0;0 0 0;0 0 0]); eps=StrainTensor(eps_m)
    @t("strain_tensor",eps_m[1,1]>0)
    @t("StrainTensor>0",eps.xx>0)
    s=hookes_law_isotropic(eps,1.0,1.0); @t("Hooke>0",s.xx>0)
    si=hookes_law_inverse(s,1.0,1.0); @t("Hooke inv",abs(si.xx-eps.xx)<1e-8)
    lam,mu=lame_constants(200e9,0.3); @t("Lame>0",lam>0&&mu>0)
    vp=p_wave_speed(1.0,1.0,1.0); vs=s_wave_speed(1.0,1.0); @t("P>S",vp>vs)
    @t("vonMises>=0",von_mises(s)>=0)
    @t("Tresca>=0",tresca(s)>=0)
    sig=StressTensor(100e6,50e6,0,25e6,0,0); @t("vonMises_val",abs(von_mises(sig)-86.602e6)/86.602e6<0.01)
    @t("principal_stresses",length(principal_stresses(sig))==3)
    @t("stress_invariants",stress_invariants(sig).I1≈150e6)
    xs,w=euler_bernoulli_deflection(1e6,10.0,x->1000.0)
    @t("beam_out",w isa Vector)
    @t("beam_len",length(w)==200)
    @t("Euler buckling",euler_buckling_load(1e6,10.0,1,1.0)>0)
    @t("plate rigidity",plate_flexural_rigidity(200e9,0.3,0.01)>0)
    @t("Rayleigh wave",rayleigh_wave_speed(1.0,1.0,1.0)>0)
    @t("Fatigue Basquin",basquin_sn_curve(1000,500e6,-0.1)>0)
    tot=p+f;println("
"*"="^60);println("  ",p,"/",tot);println(f==0?"ALL PASSED":string(f," FAILED"))
end;main();exit(f==0?0:1)