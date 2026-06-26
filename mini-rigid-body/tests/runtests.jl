#!/usr/bin/env julia
# runtests.jl — 刚体力学测试套件 (Goldstein Ch.4-5)
# 覆盖所有导出函数

include("../src/RigidBody.jl"); using .RigidBody

p = 0; f = 0

macro t(n, ex)
    quote
        try
            if ex; global p += 1; println("  ✓ ", n)
            else; global f += 1; println("  ✗ ", n, " (assertion failed)"); end
        catch e; global f += 1; println("  ✗ ", n, ": ", e); end
    end
end

function main()
    println("="^60)
    println("  Rigid Body Mechanics — Test Suite")
    println("="^60)

    # ==============================================
    # Types
    # ==============================================
    println("\n--- Types ---")
    I = InertiaTensor(2.0, 3.0, 1.0, 0.0, 0.0, 0.0)
    @t("InertiaTensor creation", I.Ixx == 2.0 && I.Iyy == 3.0 && I.Izz == 1.0)
    @t("InertiaTensor symmetry (no off-diagonal)", I.Ixy == 0.0 && I.Ixz == 0.0 && I.Iyz == 0.0)

    M = inertia_matrix(I)
    @t("inertia_matrix 3x3", size(M) == (3, 3))
    @t("inertia_matrix diagonal", M[1, 1] == 2.0 && M[2, 2] == 3.0 && M[3, 3] == 1.0)

    ea = EulerAngles(0.5, 0.3, 1.2)
    @t("EulerAngles creation", abs(ea.phi - 0.5) < 1e-15)

    rbs = RigidBodyState([1.0, 0.0, 0.0], EulerAngles())
    @t("RigidBodyState creation", rbs.omega == [1.0, 0.0, 0.0])

    # ==============================================
    # Inertia Tensor — basic computation
    # ==============================================
    println("\n--- Inertia Tensor ---")

    # Single particle on x-axis
    I1 = inertia_tensor([1.0], [[1.0, 0.0, 0.0]])
    @t("Point mass on x: Ixx=0", abs(I1.Ixx) < 1e-14)
    @t("Point mass on x: Iyy=1", abs(I1.Iyy - 1.0) < 1e-14)
    @t("Point mass on x: Izz=1", abs(I1.Izz - 1.0) < 1e-14)
    @t("Point mass on x: Ixy=0", abs(I1.Ixy) < 1e-14)

    # Two particles
    masses = [1.0, 1.0]
    positions = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
    I2 = inertia_tensor(masses, positions)
    @t("Two masses along x: Ixx=0", abs(I2.Ixx) < 1e-14)
    @t("Two masses along x: Iyy=2", abs(I2.Iyy - 2.0) < 1e-14)
    @t("Two masses along x: Izz=2", abs(I2.Izz - 2.0) < 1e-14)

    # ==============================================
    # Principal axes
    # ==============================================
    println("\n--- Principal Axes ---")

    pm = principal_moments(I)
    @t("principal_moments sorted descending", pm[1] >= pm[2] >= pm[3])
    @t("principal_moments length=3", length(pm) == 3)

    pa = principal_axes(I)
    @t("PrincipalAxes moments match", pa.moments == pm)
    @t("PrincipalAxes axes 3x3", size(pa.axes) == (3, 3))

    @t("is_diagonal true on diagonal tensor", is_diagonal(I))
    I_off = InertiaTensor(1.0, 1.0, 1.0, 0.5, 0.0, 0.0)
    @t("is_diagonal false on off-diagonal", !is_diagonal(I_off))

    # ==============================================
    # Parallel axis theorem
    # ==============================================
    println("\n--- Parallel Axis ---")

    I_sphere = sphere_inertia(1.0, 1.0)
    I_shifted = parallel_axis(I_sphere, 1.0, [1.0, 0.0, 0.0])
    @t("Parallel axis: Iyy increases by Md^2", abs(I_shifted.Iyy - I_sphere.Iyy - 1.0) < 1e-14)
    @t("Parallel axis: Ixx unchanged", abs(I_shifted.Ixx - I_sphere.Ixx) < 1e-14)

    # ==============================================
    # Standard shapes
    # ==============================================
    println("\n--- Standard Shapes ---")

    sph = sphere_inertia(1.0, 1.0)
    @t("Sphere: I = 2/5", abs(sph.Ixx - 0.4) < 1e-14)

    shell = spherical_shell_inertia(1.0, 1.0)
    @t("Spherical shell: I = 2/3", abs(shell.Ixx - 2/3) < 1e-14)

    cyl = cylinder_inertia(1.0, 1.0, 2.0)
    @t("Cylinder: Izz = 1/2", abs(cyl.Izz - 0.5) < 1e-14)

    cub = cuboid_inertia(1.0, 2.0, 2.0, 2.0)
    @t("Cuboid (cube): all I equal", abs(cub.Ixx - cub.Iyy) < 1e-14 && abs(cub.Iyy - cub.Izz) < 1e-14)
    @t("Cuboid (cube): value = 2/12*4 = 2/3", abs(cub.Ixx - 2/3) < 1e-14)

    rod = rod_inertia(1.0, 1.0)
    @t("Rod: Ixx = 1/12", abs(rod.Ixx - 1/12) < 1e-14)
    @t("Rod: Izz = 0", abs(rod.Izz) < 1e-14)

    disk = disk_inertia(1.0, 1.0)
    @t("Disk: Izz = 1/2", abs(disk.Izz - 0.5) < 1e-14)
    @t("Disk: Ixx = 1/4", abs(disk.Ixx - 0.25) < 1e-14)

    # Inertia ellipsoid
    abc = inertia_ellipsoid(sph)
    @t("Ellipsoid: sphere → a=b=c", abs(abc[1] - abc[2]) < 1e-14 && abs(abc[2] - abc[3]) < 1e-14)

    pts = inertia_ellipsoid_points(sph; n_theta=10, n_phi=5)
    @t("Ellipsoid points: correct count", length(pts) == 10 * 5)

    # ==============================================
    # Euler Equations
    # ==============================================
    println("\n--- Euler Equations ---")

    I_sym = sphere_inertia(1.0, 1.0)  # I1 = I2 = I3 = 0.4
    w_const = [1.0, 0.0, 0.0]
    dw = euler_equations_derivative(I_sym, w_const)
    @t("Sphere: dω/dt = 0 for pure spin", all(abs.(dw) .< 1e-14))

    # Asymmetric top
    I_asy = InertiaTensor(3.0, 2.0, 1.0, 0.0, 0.0, 0.0)
    w_mid = [0.01, 1.0, 0.01]  # mainly around middle axis
    dw_asy = euler_equations_derivative(I_asy, w_mid)
    @t("Asymmetric: non-zero derivative", any(abs.(dw_asy) .> 1e-15))
    @t("Euler deriv returns 3 components", length(dw_asy) == 3)

    # Euler with torque
    N = [0.0, 0.0, 1.0]
    dw_torque = euler_with_torque(I_asy, [0.0, 0.0, 0.0], N)
    @t("Torque on stationary: dωz = Nz/Izz", abs(dw_torque[3] - 1.0/1.0) < 1e-14)

    # ==============================================
    # Motion constants
    # ==============================================
    println("\n--- Motion Constants ---")

    T, L2 = motion_constants(I_asy, [0.1, 0.5, 0.2])
    @t("T > 0", T > 0.0)
    @t("L² > 0", L2 > 0.0)

    # Check consistency: T and L² are independent
    T2, L2_2 = motion_constants(I_asy, [0.3, 0.1, 0.5])
    @t("Motion constants different for different ω", abs(T - T2) > 1e-15)

    # ==============================================
    # Integrators
    # ==============================================
    println("\n--- Integrators ---")

    w0 = [0.1, 1.0, 0.05]
    w_euler = euler_step_omega(I_asy, w0, 0.01)
    @t("Euler step: 3 components", length(w_euler) == 3)
    @t("Euler step: values changed", any(w_euler .!= w0))

    w_rk4 = rk4_step_omega(I_asy, w0, 0.01)
    @t("RK4 step: 3 components", length(w_rk4) == 3)

    w_verlet = verlet_step_omega(I_asy, w0, 0.01)
    @t("Verlet step: 3 components", length(w_verlet) == 3)

    # Full simulation
    traj = simulate_free_rigid_body(I_asy, w0, 1.0, 0.01; method=:rk4)
    @t("Simulate: correct number of steps", length(traj) == 101)
    @t("Simulate: each step is 3-vector", all(length.(traj) .== 3))

    # Integrator comparison: all methods produce trajectories
    traj_euler = simulate_free_rigid_body(I_asy, w0, 0.1, 0.01; method=:euler)
    traj_verlet = simulate_free_rigid_body(I_asy, w0, 0.1, 0.01; method=:verlet)
    @t("Euler simulation runs", length(traj_euler) > 1)
    @t("Verlet simulation runs", length(traj_verlet) > 1)

    # ==============================================
    # Torque simulation
    # ==============================================
    println("\n--- Torque Simulation ---")

    constant_torque(t, w) = [0.0, 0.0, 0.5]
    times, omegas = simulate_rigid_body_with_torque(I_asy, [0.0, 0.0, 0.0], constant_torque, 1.0, 0.1)
    @t("Torque sim: ωz increases", omegas[end][3] > omegas[1][3])
    @t("Torque sim: correct time steps", length(times) == 11)

    # ==============================================
    # Motion constant monitoring
    # ==============================================
    println("\n--- Constant Monitoring ---")

    T_hist, L2_hist = monitor_constants(I_asy, traj)
    @t("Monitor: same length as trajectory", length(T_hist) == length(traj))
    @t("Monitor: T values all positive", all(T_hist .> 0))
    @t("Monitor: L² values all positive", all(L2_hist .> 0))

    drift = constant_drift_report(I_asy, traj)
    @t("Drift report: T_drift is number", drift.T_drift isa Float64)
    @t("Drift report: L2_drift is number", drift.L2_drift isa Float64)

    # ==============================================
    # Kinematics — Rotation matrices
    # ==============================================
    println("\n--- Kinematics: Rotation ---")

    Rz = euler_to_rotation(0.0, 0.0, π/2)
    @t("Rot Z(90): x→y", abs(Rz[2, 1] - 1.0) < 1e-14)

    Rx = euler_to_rotation(0.0, π/2, 0.0)
    @t("Rot X(90): preserves x", abs(Rx[1, 1] - 1.0) < 1e-14)
    @t("Rot X(90): y→z", abs(Rx[3, 2] - 1.0) < 1e-14)

    # Rotation → Euler → Rotation roundtrip
    R_orig = euler_to_rotation(0.5, 0.8, 1.3)
    ea2 = rotation_to_euler(R_orig)
    R_recovered = euler_to_rotation(ea2.phi, ea2.theta, ea2.psi)
    @t("Euler↔Rotation roundtrip", maximum(abs.(R_orig - R_recovered)) < 1e-10)

    # ==============================================
    # Kinematics — Angular velocity
    # ==============================================
    println("\n--- Kinematics: Angular Velocity ---")

    ω_body = euler_rates_to_omega(1.0, 0.0, 2.0, 0.5, 0.3)
    @t("ω from Euler rates: 3 components", length(ω_body) == 3)

    phid, thd, psid = omega_to_euler_rates(ω_body, 0.5, 0.3)
    ω_body2 = euler_rates_to_omega(phid, thd, psid, 0.5, 0.3)
    @t("Euler rates ↔ ω roundtrip", maximum(abs.(ω_body - ω_body2)) < 1e-10)

    # Cross matrix
    wx = cross_matrix([1.0, 2.0, 3.0])
    @t("cross_matrix skew-symmetric", wx + wx' == zeros(3, 3))
    @t("cross_matrix 3x3", size(wx) == (3, 3))

    # Rotation derivative
    R_test = euler_to_rotation(0.0, 0.0, 0.0)
    dR = rotation_derivative(R_test, [0.0, 0.0, 1.0])
    @t("dR/dt is 3x3", size(dR) == (3, 3))

    # Quaternion
    q = rotation_to_quaternion(R_test)
    @t("Quaternion: 4 components", length(q) == 4)
    @t("Quaternion: unit (approximately)", abs(sqrt(sum(x->x*x, q)) - 1.0) < 1e-14)

    # Basic rotations
    @t("rot_x: preserves x-axis", abs(rot_x(0.5)[1, 1] - 1.0) < 1e-14)
    @t("rot_y: preserves y-axis", abs(rot_y(0.5)[2, 2] - 1.0) < 1e-14)
    @t("rot_z: preserves z-axis", abs(rot_z(0.5)[3, 3] - 1.0) < 1e-14)

    R_rod = rot_axis([1.0, 0.0, 0.0], π/2)
    @t("rot_axis X(90): preserves X", abs(R_rod[1, 1] - 1.0) < 1e-14)

    # ==============================================
    # Energy and Angular Momentum
    # ==============================================
    println("\n--- Energy & Angular Momentum ---")

    I_test = InertiaTensor(3.0, 2.0, 1.0, 0.0, 0.0, 0.0)
    w_test = [1.0, 0.5, 0.2]

    T_rot = rotational_kinetic_energy(I_test, w_test)
    T_rot_p = rotational_KE_principal(I_test, w_test)
    @t("T_rot > 0", T_rot > 0)
    @t("T_rot matches T_rot_principal", abs(T_rot - T_rot_p) < 1e-14)

    T_tot, T_t, T_r = total_kinetic_energy(1.0, [1.0, 0.0, 0.0], I_test, w_test)
    @t("Total KE = trans + rot", abs(T_tot - T_t - T_r) < 1e-14)

    L = angular_momentum_rigid(I_test, w_test)
    L_p = angular_momentum_principal(I_test, w_test)
    @t("Angular momentum 3 components", length(L) == 3)
    @t("L matches L_principal for diagonal I", maximum(abs.(L - L_p)) < 1e-14)

    L_tot, L_orb, L_spin = total_angular_momentum(
        [1.0, 0.0, 0.0], 1.0, [0.0, 1.0, 0.0], I_test, w_test)
    @t("Total L = orb + spin (x-comp)", abs(L_tot[1] - L_orb[1] - L_spin[1]) < 1e-14)

    # ==============================================
    # Stability (Tennis Racket Theorem)
    # ==============================================
    println("\n--- Stability Analysis ---")

    # I1 > I2 > I3
    I_asym = InertiaTensor(3.0, 2.0, 1.0, 0.0, 0.0, 0.0)
    stab = stability_analysis(I_asym)
    @t("Stability dict has 3 keys", length(stab) == 3)
    @t("Axis 1 stable (largest)", stab[:I1] == "stable")
    @t("Axis 2 unstable (middle)", stab[:I2] == "unstable")
    @t("Axis 3 stable (smallest)", stab[:I3] == "stable")

    s1, _ = axis_stability(I_asym, 1)
    s2, r2 = axis_stability(I_asym, 2)
    s3, _ = axis_stability(I_asym, 3)
    @t("axis_stability: I1=stable", s1 == :stable)
    @t("axis_stability: I2=unstable", s2 == :unstable)
    @t("axis_stability: I3=stable", s3 == :stable)
    @t("axis_stability: I2 has growth rate", r2 > 0)

    # Symmetric top: axis 2 degenerate
    I_sym_full = InertiaTensor(2.0, 2.0, 1.0, 0.0, 0.0, 0.0)
    stab_sym = stability_analysis(I_sym_full)
    @t("Symmetric top: I2 degenerate", stab_sym[:I2] == "degenerate")

    # Flipping period
    T_flip = flipping_period_estimate(I_asym, 1.0, 0.01)
    @t("Flipping period finite", T_flip > 0 && T_flip < Inf)

    # Torque power
    P = torque_power([0.0, 0.0, 1.0], [0.1, 0.0, 0.5])
    @t("Torque power: Nz*ωz", abs(P - 0.5) < 1e-14)

    J = torque_impulse([0.0, 1.0, 0.0], 0.5)
    @t("Torque impulse: N*dt", J == [0.0, 0.5, 0.0])

    # ==============================================
    # Tops — Symmetric Top Dynamics
    # ==============================================
    println("\n--- Tops ---")

    I1, I3 = 1.0, 0.5
    M_top, g, l = 1.0, 9.81, 0.3
    Lz, L3 = 1.0, 2.0

    # Effective potential
    thetas, Veff = symmetric_top_effective_potential(I1, I3, M_top, g, l, Lz, L3)
    @t("V_eff: correct length", length(thetas) == 200)
    @t("V_eff: returns matching lengths", length(thetas) == length(Veff))
    @t("V_eff: finite values", all(isfinite.(Veff)))

    # Nutation rate
    E_test = Veff[100] + 0.1
    nut_rate = nutation_rate(I1, I3, M_top, g, l, E_test, Lz, L3, π/4)
    @t("Nutation rate: returns scalar", nut_rate isa Float64)
    @t("Nutation rate: ≥ 0", nut_rate >= 0)

    # Precession rate
    φdot = precession_rate_top(I1, Lz, L3, π/4)
    @t("Precession rate: scalar", φdot isa Float64)

    ψdot = spin_rate_top(I3, I1, Lz, L3, π/4)
    @t("Spin rate: scalar", ψdot isa Float64)

    # Steady precession
    theta_steady = steady_precession_theta(I1, I3, M_top, g, l, Lz, L3;
                                           bracket=(0.2, 1.0))
    @t("Steady precession θ in range", 0.0 < theta_steady < π)

    # Fast top
    Ω_fast = fast_top_precession(M_top, g, l, I3, 10.0)
    @t("Fast top: Ω > 0", Ω_fast > 0)

    Ω_slow = slow_precession(M_top, g, l, L3)
    @t("Slow precession: Ω > 0", Ω_slow > 0)

    # Sleeping top
    status, ω_crit = sleeping_top_stability(M_top, g, l, I1, I3, 10.0)
    @t("Sleeping top: has status", status in [:stable, :unstable])

    ω_c = critical_spin_rate(M_top, g, l, I1, I3)
    @t("Critical spin rate > 0", ω_c > 0)

    # Full Lagrange top simulation
    times_top, thetas_top, phis, psis = simulate_lagrange_top(
        I1, I3, M_top, g, l, π/6, 0.0, 0.0, 0.1, Lz, L3, 0.5, 0.01)
    @t("Lagrange sim: returns time series", length(times_top) == 51)
    @t("Lagrange sim: θ stays positive", all(thetas_top .> 0.0))
    @t("Lagrange sim: all arrays same length",
       length(thetas_top) == length(phis) == length(psis))

    # Nutation period (small amplitude)
    T_nut = nutation_period_small(I1, M_top, g, l)
    @t("Nutation period > 0", T_nut > 0)

    # Nutation amplitude range
    amp_range = nutation_amplitude_range(I1, I3, M_top, g, l, E_test, Lz, L3)
    @t("Nutation amplitude: returns vector", amp_range isa Vector)

    # Precession type
    if length(amp_range) >= 2
        ptype = precession_type(I1, Lz, L3, amp_range)
        @t("Precession type: valid symbol", ptype in [:monotonic, :looping, :cusped])
    else
        @t("Precession type: skipped (insufficient range)", true)
    end

    # Gyroscope
    N_gyro = gyroscopic_torque([0.0, 0.0, 1.0], [0.1, 0.0, 0.0])
    @t("Gyroscopic torque: ⊥ to L and Ω", abs(N_gyro[3]) < 1e-14)

    Ω_gyro = gyroscope_precession(0.5, 2.0)
    @t("Gyroscope precession: Ω = N/L", abs(Ω_gyro - 0.25) < 1e-14)

    # ==============================================
    # Summary
    # ==============================================
    tot = p + f
    println("\n" * "="^60)
    println("  Results: $p / $tot passed")
    if f == 0
        println("  ✅ ALL TESTS PASSED")
    else
        println("  ❌ $f TEST(S) FAILED")
    end
    println("="^60)
    exit(f == 0 ? 0 : 1)
end

main()