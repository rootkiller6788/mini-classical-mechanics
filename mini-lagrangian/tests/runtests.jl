#!/usr/bin/env julia
# tests/runtests.jl — mini-lagrangian 测试

include("../src/Lagrangian.jl")
using .Lagrangian
using LinearAlgebra

passed = 0; failed = 0

macro test(name, expr)
    quote
        try
            if $(esc(expr))
                global passed += 1; println("  ✅ ", $(esc(name)))
            else
                global failed += 1; println("  ❌ ", $(esc(name)), " — false")
            end
        catch e
            global failed += 1; println("  ❌ ", $(esc(name)), " — ", e)
        end
    end
end

function main()
    println("="^60)
    println("  mini-lagrangian Test Suite")
    println("="^60)

    # ---- Generalized Coords ----
    println("\n--- Generalized Coordinates ---")
    gc = GeneralizedCoords([1.0, 2.0, 3.0])
    @test("GeneralizedCoords creation", gc.n == 3 && gc.values == [1.0,2.0,3.0])
    gv = GeneralizedVelocities([0.5, -1.0])
    @test("GeneralizedVelocities creation", gv.n == 2)

    # Polar coord transform
    x, y = polar_to_cartesian(1.0, π/2)
    @test("polar_to_cartesian: (1,π/2)->(0,1)", abs(x) < 1e-10 && abs(y-1) < 1e-10)
    J = polar_jacobian(2.0, π/4)
    @test("polar_jacobian size", size(J) == (2,2))
    @test("polar_jacobian ∂x/∂r", abs(J[1,1] - cos(π/4)) < 1e-10)

    # Spherical
    x,y,z = spherical_to_cartesian(1.0, π/2, 0.0)
    @test("spherical: (1,π/2,0)->(1,0,0)", abs(x-1) < 1e-10 && abs(y) < 1e-10 && abs(z) < 1e-10)

    # Mass matrix
    masses = [1.0, 1.0]
    J2 = [1.0 0.0; 0.0 1.0; 0.0 0.0; 0.0 0.0]  # 2 particles in 2D, no coupling
    M = mass_matrix(masses, J2)
    @test("mass matrix diag", abs(M[1,1]-1) < 1e-10 && abs(M[1,2]) < 1e-10)

    # ---- Euler-Lagrange Solver ----
    println("\n--- EL Solver ---")
    # SHO: L = 0.5*m*xdot^2 - 0.5*k*x^2
    m,k = 1.0, 1.0
    L_sho(q,qdot) = 0.5*m*qdot[1]^2 - 0.5*k*q[1]^2
    grad_q(q,qdot) = [-k*q[1]]
    grad_qdot(q,qdot) = [m*qdot[1]]
    M_func(q) = [m]
    sys = EulerLagrangeSystem(1, L_sho, grad_q, grad_qdot, M_func)
    f = el_to_first_order(sys)

    y0 = [1.0, 0.0]; dt = 0.01
    y = rk4_step(f, 0.0, y0, dt)
    @test("EL: SHO x decreasing from 1", y[1] < 1.0)

    # ---- Noether ----
    println("\n--- Noether's Theorem ---")
    E0 = energy_from_lagrangian(L_sho, [1.0], [0.0], grad_qdot)
    @test("Energy initial = 0.5", abs(E0 - 0.5) < 1e-10)

    # Time translation → energy conservation
    energies = Float64[]
    y = copy(y0)
    for _ in 1:100
        global y = rk4_step(f, 0.0, y, 0.01)
        push!(energies, energy_from_lagrangian(L_sho, [y[1]], [y[2]], grad_qdot))
    end
    E_drift = maximum(abs.(energies .- E0) ./ abs(E0))
    @test("SHO energy conservation (RK4)", E_drift < 1e-4)

    # Translation symmetry Q
    Q_trans = translation_symmetry_Q(2, 1)
    @test("translation Q for 2 particles", length(Q_trans)==6 && Q_trans[1]==1 && Q_trans[4]==1)

    # Rotation symmetry Q_z
    Q_rot = rotation_symmetry_Q_z([[1.0,2.0,3.0], [4.0,5.0,6.0]])
    @test("rotation Q_z: particle 1", abs(Q_rot[1]+2.0)<1e-10 && abs(Q_rot[2]-1.0)<1e-10)

    # ---- Action ----
    println("\n--- Action Principle ---")
    function q_true(t); return [cos(t)]; end
    function qdot_true(t); return [-sin(t)]; end
    function q_varied(t); return [cos(t) + 0.1*sin(3t)]; end
    function qdot_varied(t); return [-sin(t) + 0.3*cos(3t)]; end
    S_true, S_var, diff = verify_least_action(L_sho, q_true, qdot_true, q_varied, qdot_varied, (0.0, 2π))
    @test("Action: true path < varied path", S_true < S_var)

    # ---- Constraints ----
    println("\n--- Constraints ---")
    c = spherical_constraint(1.0)
    @test("spherical constraint: on surface", abs(c.func([1.0,0.0,0.0], 0.0)) < 1e-10)
    @test("spherical constraint: at 2R", abs(c.func([2.0,0.0,0.0], 0.0) - 3.0) < 1e-10)
    @test("holonomic flag", c.is_holonomic)

    rod = rigid_rod_constraint(3.0)
    f_val = rod.func([0.0,0.0,0.0, 3.0,0.0,0.0], 0.0)
    @test("rigid rod: |r1-r2|=3", abs(f_val) < 1e-10)

    # ---- Legendre Transform ----
    println("\n--- Legendre Transform ---")
    M_func(q) = [1.0 0.0; 0.0 1.0]
    U_func(q) = 0.5*(q[1]^2 + q[2]^2)
    grad_U(q) = [q[1], q[2]]
    lt = standard_legendre(M_func, U_func, grad_U, 2)
    q_test = [1.0, 0.0]; qdot_test = [2.0, 3.0]
    p = qdot_to_momentum(lt, q_test, qdot_test)
    @test("qdot→p in identity M", norm(p - [2.0,3.0]) < 1e-10)
    qdot_back = momentum_to_qdot(lt, q_test, p)
    @test("p→qdot roundtrip", norm(qdot_back - qdot_test) < 1e-10)
    H_val = hamiltonian(lt, q_test, p)
    T = 0.5*(4+9); U = 0.5*(1+0)
    @test("H=T+U", abs(H_val - (T+U)) < 1e-10)
    # Hamilton's equations test
    f_ham = hamiltons_equations(lt)
    y0 = [1.0, 0.0, 0.0, 1.0]; dt=0.01
    y1 = rk4_step(f_ham, 0.0, y0, dt)
    @test("Hamilton: q1 decreases from 1", y1[1] < 1.0)

    # Poisson bracket {q1, p2} = 0
    pb = poisson_bracket((q,p)->q[1], (q,p)->p[2], [1.0,2.0], [3.0,4.0])
    @test("{q1,p2}=0", abs(pb) < 1e-6)
    pb2 = poisson_bracket((q,p)->q[1], (q,p)->p[1], [1.0,2.0], [3.0,4.0])
    @test("{q1,p1}=1", abs(pb2 - 1.0) < 0.01)

    # ---- Small Oscillations ----
    println("\n--- Small Oscillations ---")
    U_sho(q) = 0.5*(q[1]^2 + q[2]^2 + (q[1]-q[2])^2)  # coupled
    M_sho(q) = [1.0 0.0; 0.0 1.0]
    grad_U_sho(q) = numerical_gradient(U_sho, q)
    hess_U_sho(q) = numerical_hessian(U_sho, q)
    q_eq = find_equilibrium(grad_U_sho, hess_U_sho, [0.5, -0.3])
    @test("equilibrium at origin", norm(q_eq) < 1e-6)
    K_sho = hess_U_sho(q_eq)
    M_sho_val = M_sho(q_eq)
    sys_sho = solve_normal_modes(M_sho_val, K_sho, q_eq)
    @test("2 normal modes found", length(sys_sho.frequencies) == 2)
    @test("frequencies positive", all(sys_sho.frequencies .> 0))
    # M-normalization: a'M a = 1
    for a in 1:2
        nrm = dot(sys_sho.modes[:,a], M_sho_val * sys_sho.modes[:,a])
        @test("mode $a M-normalized", abs(nrm - 1.0) < 1e-10)
    end

    # ---- Summary ----
    println("\n" * "="^60)
    total = passed + failed
    println("  Results: $passed / $total passed")
    println(failed == 0 ? "  🎉 ALL TESTS PASSED" : "  ⚠️  $failed FAILED")
    println("="^60)
    return failed == 0
end

exit(main() ? 0 : 1)
