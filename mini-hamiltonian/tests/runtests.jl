#!/usr/bin/env julia
# tests/runtests.jl — mini-hamiltonian tests
include("../src/Hamiltonian.jl")
using .Hamiltonian

passed=0; failed=0
macro test(n,e)
    quote try
        if $(esc(e)); global passed+=1; println("  ✅ ",$(esc(n)))
        else; global failed+=1; println("  ❌ ",$(esc(n)))
        end
    catch e; global failed+=1; println("  ❌ ",$(esc(n))," — ",e); end
    end
end

function main()
    println("="^60)
    println("  mini-hamiltonian Test Suite")
    println("="^60)

    # Phase space
    println("\n--- Phase Space ---")
    pp = PhasePoint([1.0,0.0], [0.0,1.0])
    @test("PhasePoint creation", pp.n==2)
    traj=PhaseTrajectory(); record!(traj,0.0,[1.0],[0.0],0.5)
    @test("record", length(traj.ts)==1 && traj.H_vals[1]==0.5)

    # Hamilton's equations: SHO
    println("\n--- Hamilton Equations ---")
    H(q,p)=0.5*(p[1]^2+q[1]^2)
    sys = analytic_hamiltonian_system(H, q->[q[1]], p->[p[1]], 1)
    f=hamiltons_ode(sys)
    y=rk4_step(f,0.0,[1.0,0.0],0.01)
    @test("SHO: q decreases", y[1]<1.0 && y[2]<0.0)

    # Symplectic Euler vs Verlet
    q0=[1.0]; p0=[0.0]; dt=0.01
    q_se,p_se=symplectic_euler_step(sys,q0,p0,dt)
    q_sv,p_sv=stormer_verlet_step(sys,q0,p0,dt)
    @test("both methods advance", q_se[1]!=q0[1] && q_sv[1]!=q0[1])

    # Energy conservation
    traj_sv=solve_hamiltonian(sys,q0,p0,2π,0.01;method=:verlet)
    Edrift=max(abs.(traj_sv.H_vals.-traj_sv.H_vals[1]))
    @test("Verlet energy conserved (SHO)", Edrift<1e-4)

    # Poisson brackets
    println("\n--- Poisson Brackets ---")
    pb1=poisson_bracket((q,p)->q[1],(q,p)->p[1],[1.0],[2.0])
    @test("{q1,p1}=1", abs(pb1-1.0)<0.01)
    pb2=poisson_bracket((q,p)->q[1],(q,p)->q[1],[1.0],[2.0])
    @test("{q1,q1}=0", abs(pb2)<1e-10)
    errs=fundamental_poisson_brackets(2,[1.0,2.0],[3.0,4.0])
    @test("fundamental brackets", max(maximum(errs[1]),maximum(errs[2]),maximum(errs[3]))<1e-4)

    # Jacobi identity
    f(q,p)=q[1]; g(q,p)=p[1]; h(q,p)=q[1]^2
    jac=verify_jacobi_identity(f,g,h,[1.0],[2.0])
    @test("Jacobi identity", abs(jac)<1e-6)

    # Canonical transforms
    println("\n--- Canonical Transforms ---")
    ct_id=identity_transform(2)
    ok,err=verify_canonical(ct_id,[1.0,2.0],[3.0,4.0])
    @test("identity is canonical", ok)
    ct_scale=scale_transform(2.0,2)
    ok_s,_=verify_canonical(ct_scale,[1.0,2.0],[3.0,4.0])
    @test("scale(2) is canonical (μν=1)", ok_s)

    # Liouville
    println("\n--- Liouville ---")
    ensemble=generate_ensemble([1.0],[0.0],30)
    @test("generate ensemble", length(ensemble)==30)
    V=ensemble_volume_proxy(ensemble)
    @test("volume proxy >0", V>0)

    # Action-angle
    println("\n--- Action-Angle ---")
    aa=HarmonicOscillatorActionAngle(1.0,1.0)
    J,th=qp_to_action_angle_ho(aa,1.0,0.0)
    @test("HO action J=E/ω", abs(J-0.5)<1e-10)
    E=energy_from_action_ho(aa,0.5)
    @test("HO E=ωJ", abs(E-0.5)<1e-10)
    q,p=action_to_qp_ho(aa,0.5,π/2)
    @test("HO q from J,θ=π/2", abs(q-1.0)<1e-10 && abs(p)<1e-10)

    println("\n"*"="^60)
    tot=passed+failed
    println("  Results: $passed / $tot passed")
    println(failed==0 ? "  🎉 ALL TESTS PASSED" : "  ⚠️  $failed FAILED")
    return failed==0
end
exit(main()?0:1)
