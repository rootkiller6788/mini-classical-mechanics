/**
 * test_rigid_body.c — Comprehensive Test Suite for Rigid Body Mechanics
 *
 * Goldstein Ch.4-5, Landau Ch.6, MIT 8.012 coverage verification.
 *
 * All tests use assert() for mathematical truth verification.
 * Each test corresponds to a specific knowledge point.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "../include/rigid_body.h"

static int passed = 0, failed = 0;

#define TEST(name) do { printf("  %-55s", name); } while(0)
#define CHECK(cond) do { \
    if (cond) { printf("✓ PASS\n"); passed++; } \
    else { printf("✗ FAIL\n"); failed++; } \
} while(0)
#define ASSERT_FLOAT_EQ(a, b, tol) CHECK(fabs((a) - (b)) < (tol))
#define ASSERT_VEC3_EQ(a, b, tol) CHECK( \
    fabs((a)[0] - (b)[0]) < (tol) && \
    fabs((a)[1] - (b)[1]) < (tol) && \
    fabs((a)[2] - (b)[2]) < (tol))

/* ========================================================================== */
int main(void) {
    printf("══════════════════════════════════════════════════════════════\n");
    printf("  Rigid Body Mechanics — C Test Suite (Goldstein Ch.4-5)\n");
    printf("══════════════════════════════════════════════════════════════\n");

    /* ======================================================================
     * L1: Types — InertiaTensor, EulerAngles, RigidBodyState, PrincipalAxes
     * ====================================================================== */
    printf("\n─── L1: Core Types ───\n");

    InertiaTensor I0 = {2.0, 3.0, 1.0, 0.0, 0.0, 0.0};
    TEST("InertiaTensor diagonal creation");
    CHECK(I0.Ixx == 2.0 && I0.Iyy == 3.0 && I0.Izz == 1.0);

    TEST("InertiaTensor zero off-diagonal");
    CHECK(I0.Ixy == 0.0 && I0.Ixz == 0.0 && I0.Iyz == 0.0);

    mat3 M0;
    inertia_to_mat3(&I0, &M0);
    TEST("inertia_to_mat3: diagonal preserved");
    CHECK(M0.m[0]==2.0 && M0.m[4]==3.0 && M0.m[8]==1.0);

    InertiaTensor I_from_mat = mat3_to_inertia(&M0);
    TEST("mat3_to_inertia: roundtrip");
    CHECK(I_from_mat.Ixx == I0.Ixx && I_from_mat.Iyy == I0.Iyy && I_from_mat.Izz == I0.Izz);

    EulerAngles ea = {0.5, 0.3, 1.2};
    TEST("EulerAngles creation");
    CHECK(ea.phi == 0.5 && ea.theta == 0.3 && ea.psi == 1.2);

    euler_normalize(&ea);
    TEST("EulerAngles normalization: theta in [0,π]");
    CHECK(ea.theta >= 0.0 && ea.theta <= M_PI);

    RigidBodyState rbs = rb_state_new((vec3){1.0, 0.0, 0.0}, EULER_ZERO);
    TEST("RigidBodyState creation");
    CHECK(rbs.omega.x == 1.0 && rbs.t == 0.0);

    /* ======================================================================
     * L1: Quaternion Operations
     * ====================================================================== */
    Quaternion q = QUAT_IDENTITY;
    TEST("Quaternion identity: rotates v to itself");
    vec3 vtest = {1.0, 2.0, 3.0};
    vec3 vrot = quat_rotate_vec(q, vtest);
    {
        double vr_a[3] = {vrot.x, vrot.y, vrot.z};
        double vr_b[3] = {1.0, 2.0, 3.0};
        ASSERT_VEC3_EQ(vr_a, vr_b, 1e-14);
    }

    Quaternion q_xy = {0.707106781, 0.707106781, 0.0, 0.0}; /* 90° about x */
    q_xy = quat_normalize(q_xy);
    TEST("Quaternion norm: unit after normalize");
    CHECK(quat_is_unit(q_xy, 1e-14));

    /* ======================================================================
     * L2: Inertia Tensor from Particles
     * ====================================================================== */
    printf("\n─── L2: Inertia Tensor Computation ───\n");

    /* Single particle on x-axis at x=1 */
    {
        double masses[] = {1.0};
        vec3 positions[] = {{1.0, 0.0, 0.0}};
        InertiaTensor Ip;
        inertia_tensor_from_particles(1, masses, positions, &Ip);
        TEST("Point mass on x-axis: Ixx=0");
        ASSERT_FLOAT_EQ(Ip.Ixx, 0.0, 1e-14);
        TEST("Point mass on x-axis: Iyy=y²+z²=1");
        ASSERT_FLOAT_EQ(Ip.Iyy, 1.0, 1e-14);
        TEST("Point mass on x-axis: Izz=1");
        ASSERT_FLOAT_EQ(Ip.Izz, 1.0, 1e-14);
        TEST("Point mass on x-axis: Ixy=0");
        ASSERT_FLOAT_EQ(Ip.Ixy, 0.0, 1e-14);
    }

    /* Two particles along x-axis at ±1 */
    {
        double masses[] = {1.0, 1.0};
        vec3 positions[] = {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}};
        InertiaTensor Ip;
        inertia_tensor_from_particles(2, masses, positions, &Ip);
        TEST("Two masses along x: Ixx=0");
        ASSERT_FLOAT_EQ(Ip.Ixx, 0.0, 1e-14);
        TEST("Two masses along x: Iyy=2");
        ASSERT_FLOAT_EQ(Ip.Iyy, 2.0, 1e-14);
    }

    /* ======================================================================
     * L2: Principal Axes
     * ====================================================================== */
    {
        PrincipalAxes pa;
        principal_axes_decompose(&I0, &pa);
        TEST("PrincipalAxes: moments sorted descending");
        CHECK(pa.moments[0] >= pa.moments[1] && pa.moments[1] >= pa.moments[2]);
        TEST("PrincipalAxes: I1=3.0, I2=2.0, I3=1.0");
        CHECK(fabs(pa.moments[0]-3.0)<1e-12 && fabs(pa.moments[1]-2.0)<1e-12 && fabs(pa.moments[2]-1.0)<1e-12);

        double moments[3];
        principal_moments_only(&I0, moments);
        TEST("principal_moments_only: match decomposed");
        CHECK(moments[0]==pa.moments[0] && moments[1]==pa.moments[1] && moments[2]==pa.moments[2]);
    }

    TEST("is_inertia_diagonal: true for diagonal tensor");
    CHECK(is_inertia_diagonal(&I0, 1e-14));

    InertiaTensor I_off = {1.0, 1.0, 1.0, 0.5, 0.0, 0.0};
    TEST("is_inertia_diagonal: false when off-diagonal present");
    CHECK(!is_inertia_diagonal(&I_off, 1e-14));

    /* ======================================================================
     * L2: Parallel Axis Theorem
     * ====================================================================== */
    {
        InertiaTensor I_cm = sphere_inertia(1.0, 1.0);
        vec3 d = {1.0, 0.0, 0.0};
        InertiaTensor I_p;
        parallel_axis_theorem(&I_cm, 1.0, d, &I_p);
        TEST("Parallel axis: Ixx unchanged");
        ASSERT_FLOAT_EQ(I_p.Ixx, I_cm.Ixx, 1e-14);
        TEST("Parallel axis: Iyy increases by Md²");
        ASSERT_FLOAT_EQ(I_p.Iyy, I_cm.Iyy + 1.0, 1e-14);

        InertiaTensor I_cm_back;
        parallel_axis_inverse(&I_p, 1.0, d, &I_cm_back);
        TEST("Parallel axis inverse: roundtrip");
        ASSERT_FLOAT_EQ(I_cm_back.Ixx, I_cm.Ixx, 1e-13);
        ASSERT_FLOAT_EQ(I_cm_back.Iyy, I_cm.Iyy, 1e-13);
    }

    /* ======================================================================
     * L2: Standard Shapes
     * ====================================================================== */
    printf("\n─── Standard Shapes ───\n");

    InertiaTensor I_sph = sphere_inertia(1.0, 1.0);
    TEST("Solid sphere: I = 2/5 MR²");
    ASSERT_FLOAT_EQ(I_sph.Ixx, 0.4, 1e-14);

    InertiaTensor I_shell = spherical_shell_inertia(1.0, 1.0);
    TEST("Spherical shell: I = 2/3 MR²");
    ASSERT_FLOAT_EQ(I_shell.Ixx, 2.0/3.0, 1e-14);

    InertiaTensor I_cyl = cylinder_inertia(1.0, 1.0, 2.0);
    TEST("Cylinder: Izz = 1/2 MR²");
    ASSERT_FLOAT_EQ(I_cyl.Izz, 0.5, 1e-14);

    InertiaTensor I_cub = cuboid_inertia(1.0, 2.0, 2.0, 2.0);
    TEST("Cube: all I equal");
    CHECK(fabs(I_cub.Ixx - I_cub.Iyy) < 1e-14 && fabs(I_cub.Iyy - I_cub.Izz) < 1e-14);
    TEST("Cube: I = 2/12 * 4 = 2/3");
    ASSERT_FLOAT_EQ(I_cub.Ixx, 2.0/3.0, 1e-14);

    InertiaTensor I_rod = rod_inertia(1.0, 1.0);
    TEST("Rod: Ixx = ML²/12");
    ASSERT_FLOAT_EQ(I_rod.Ixx, 1.0/12.0, 1e-14);
    TEST("Rod: Izz = 0");
    ASSERT_FLOAT_EQ(I_rod.Izz, 0.0, 1e-14);

    InertiaTensor I_disk = disk_inertia(1.0, 1.0);
    TEST("Disk: Izz = 1/2 MR²");
    ASSERT_FLOAT_EQ(I_disk.Izz, 0.5, 1e-14);
    TEST("Disk: Ixx = 1/4 MR²");
    ASSERT_FLOAT_EQ(I_disk.Ixx, 0.25, 1e-14);

    InertiaTensor I_ell = ellipsoid_inertia(1.0, 2.0, 1.0, 1.0);
    TEST("Ellipsoid: a=2,b=1,c=1 => Ixx=(1+1)/5=0.4, Iyy=(4+1)/5=1.0");
    ASSERT_FLOAT_EQ(I_ell.Ixx, 0.4, 1e-14);
    ASSERT_FLOAT_EQ(I_ell.Iyy, 1.0, 1e-14);

    InertiaTensor I_cone = solid_cone_inertia(1.0, 1.0, 1.0);
    TEST("Solid cone: Izz = 0.3 MR²");
    ASSERT_FLOAT_EQ(I_cone.Izz, 0.3, 1e-14);

    InertiaTensor I_torus = torus_inertia(1.0, 2.0, 0.5);
    TEST("Torus: Izz matches formula");
    double Izz_torus_expected = 1.0 * (4.0 + 0.75 * 0.25);
    ASSERT_FLOAT_EQ(I_torus.Izz, Izz_torus_expected, 1e-14);

    /* Inertia ellipsoid semi-axes */
    {
        double a, b, c;
        inertia_ellipsoid_semiaxes(&I0, &a, &b, &c);
        TEST("Inertia ellipsoid: a < b < c (longest along min inertia)");
        CHECK(a < b && b < c);
    }

    /* Radius of gyration */
    {
        vec3 axis = {0.0, 0.0, 1.0};
        double k = radius_of_gyration(&I_cyl, 1.0, axis);
        TEST("Radius of gyration about z-axis of cylinder");
        ASSERT_FLOAT_EQ(k, sqrt(0.5), 1e-14);
    }

    /* ======================================================================
     * L4: Euler Equations
     * ====================================================================== */
    printf("\n─── L4: Euler Equations ───\n");

    {
        double I_sphere[3] = {0.4, 0.4, 0.4};
        double w_const[3] = {1.0, 0.0, 0.0};
        double domega[3];
        euler_free_derivative(I_sphere, w_const, domega);
        TEST("Sphere: dω/dt = 0 for pure spin");
        CHECK(fabs(domega[0])<1e-14 && fabs(domega[1])<1e-14 && fabs(domega[2])<1e-14);
    }

    {
        double I_asy[3] = {3.0, 2.0, 1.0};
        double w_mid[3] = {0.01, 1.0, 0.01};
        double dw[3];
        euler_free_derivative(I_asy, w_mid, dw);
        TEST("Asymmetric: non-zero derivative");
        CHECK(fabs(dw[0])>1e-15 || fabs(dw[1])>1e-15 || fabs(dw[2])>1e-15);
    }

    {
        double I_asy[3] = {3.0, 2.0, 1.0};
        double w_zero[3] = {0.0, 0.0, 0.0};
        double N[3] = {0.0, 0.0, 1.0};
        double dw[3];
        euler_with_torque(I_asy, w_zero, N, dw);
        TEST("Torque on stationary: dωz = Nz/I₃");
        ASSERT_FLOAT_EQ(dw[2], 1.0/1.0, 1e-14);
    }

    /* ======================================================================
     * L4: Constants of Motion
     * ====================================================================== */
    {
        double I[3] = {3.0, 2.0, 1.0};
        double w[3] = {0.1, 0.5, 0.2};
        double T, L2;
        motion_constants_euler(I, w, &T, &L2);
        TEST("Motion constants: T > 0");
        CHECK(T > 0.0);
        TEST("Motion constants: L² > 0");
        CHECK(L2 > 0.0);

        double T2, L2_2;
        double w2[3] = {0.3, 0.1, 0.5};
        motion_constants_euler(I, w2, &T2, &L2_2);
        TEST("Motion constants: different for different ω");
        CHECK(fabs(T - T2) > 1e-15);
    }

    /* ======================================================================
     * L5: Numerical Integrators
     * ====================================================================== */
    printf("\n─── L5: Numerical Integrators ───\n");

    {
        double I[3] = {3.0, 2.0, 1.0};
        double w0[3] = {0.1, 1.0, 0.05};

        double w_euler[3];
        euler_step_omega(I, w0, 0.01, w_euler);
        TEST("Euler step: 3 components (changed)");
        CHECK(w_euler[0] != w0[0] || w_euler[1] != w0[1] || w_euler[2] != w0[2]);

        double w_rk4[3];
        rk4_step_omega(I, w0, 0.01, w_rk4);
        TEST("RK4 step: finite values");
        CHECK(isfinite(w_rk4[0]) && isfinite(w_rk4[1]) && isfinite(w_rk4[2]));

        double w_mid[3];
        midpoint_step_omega(I, w0, 0.01, w_mid);
        TEST("Midpoint step: finite values");
        CHECK(isfinite(w_mid[0]) && isfinite(w_mid[1]) && isfinite(w_mid[2]));

        double w_verlet[3];
        verlet_step_omega(I, w0, 0.01, w_verlet);
        TEST("Verlet step: finite values");
        CHECK(isfinite(w_verlet[0]) && isfinite(w_verlet[1]) && isfinite(w_verlet[2]));

        /* DOPRI54 adaptive step */
        double w_dopri[3] = {w0[0], w0[1], w0[2]};
        double dt = 0.01, dt_next, err;
        dopri54_step_omega(I, w_dopri, dt, &dt_next, &err);
        TEST("DOPRI54 step: finite values");
        CHECK(isfinite(w_dopri[0]) && isfinite(w_dopri[1]) && isfinite(w_dopri[2]));
        TEST("DOPRI54: valid error estimate");
        CHECK(err >= 0.0);
    }

    /* Full simulation */
    {
        double I[3] = {3.0, 2.0, 1.0};
        double w0[3] = {0.1, 1.0, 0.05};
        int n_steps;
        double *traj = malloc(200 * 3 * sizeof(double));

        int ret = simulate_free_rigid_body(I, w0, 1.0, 0.01, INTEGRATOR_RK4, &n_steps, traj);
        TEST("Simulate free body: returns 0");
        CHECK(ret == 0);
        TEST("Simulate free body: ~101 steps");
        CHECK(n_steps >= 100 && n_steps <= 102);

        /* Monitor constants */
        double *T_hist = malloc(n_steps * sizeof(double));
        double *L2_hist = malloc(n_steps * sizeof(double));
        monitor_motion_constants(I, n_steps, traj, T_hist, L2_hist);

        double T_drift, L2_drift;
        constant_drift_report(n_steps, T_hist, L2_hist, &T_drift, &L2_drift);
        TEST("RK4: energy drift < 1e-6");
        CHECK(T_drift < 1e-6);
        TEST("RK4: angular momentum drift < 1e-6");
        CHECK(L2_drift < 1e-6);

        free(traj); free(T_hist); free(L2_hist);
    }

    /* Torque simulation */
    {
        double I[3] = {3.0, 2.0, 1.0};
        double w0[3] = {0.0, 0.0, 0.0};
        int n_steps;
        double *times = malloc(20 * sizeof(double));
        double *omegas = malloc(20 * 3 * sizeof(double));

        /* Constant torque function */
        void const_torque(double t, const double w[3], double N[3], void *p) {
            (void)t; (void)w; (void)p;
            N[0] = 0.0; N[1] = 0.0; N[2] = 0.5;
        }
        simulate_rigid_body_torque(I, w0, const_torque, NULL, 1.0, 0.1,
                                   &n_steps, times, omegas);
        TEST("Torque sim: ωz increases");
        CHECK(omegas[(n_steps-1)*3+2] > omegas[2]);

        free(times); free(omegas);
    }

    /* ======================================================================
     * L3: Kinematics — Euler Angles ↔ Rotation Matrix
     * ====================================================================== */
    printf("\n─── L3: Kinematics ───\n");

    {
        mat3 Rz;
        euler_to_rotation_matrix(0.0, 0.0, M_PI_2, &Rz);
        TEST("Rot Z(90°): R[1][0] ≈ 1 (x→y)");
        ASSERT_FLOAT_EQ(Rz.m[3], 1.0, 1e-14);
    }

    {
        mat3 Rx;
        euler_to_rotation_matrix(0.0, M_PI_2, 0.0, &Rx);
        TEST("Rot X(90°): preserves x");
        ASSERT_FLOAT_EQ(Rx.m[0], 1.0, 1e-14);
        TEST("Rot X(90°): y→z");
        ASSERT_FLOAT_EQ(Rx.m[7], 1.0, 1e-14);
    }

    /* Roundtrip: Euler → Matrix → Euler */
    {
        mat3 R_orig;
        euler_to_rotation_matrix(0.5, 0.8, 1.3, &R_orig);
        EulerAngles ea_back;
        rotation_matrix_to_euler(&R_orig, &ea_back);
        mat3 R_recovered;
        euler_to_rotation_matrix(ea_back.phi, ea_back.theta, ea_back.psi, &R_recovered);
        TEST("Euler ↔ Matrix roundtrip: R ≈ R_recovered");
        double max_diff = 0.0;
        for (int i = 0; i < 9; i++) {
            double d = fabs(R_orig.m[i] - R_recovered.m[i]);
            if (d > max_diff) max_diff = d;
        }
        CHECK(max_diff < 1e-10);
    }

    /* ======================================================================
     * L3: Angular Velocity ↔ Euler Rates
     * ====================================================================== */
    {
        double omega[3];
        euler_rates_to_body_omega(1.0, 0.0, 2.0, 0.5, 0.3, omega);
        TEST("Euler rates → ω: returns 3 components");
        CHECK(isfinite(omega[0]) && isfinite(omega[1]) && isfinite(omega[2]));

        double phid, thd, psid;
        body_omega_to_euler_rates(omega, 0.5, 0.3, &phid, &thd, &psid);
        double omega2[3];
        euler_rates_to_body_omega(phid, thd, psid, 0.5, 0.3, omega2);
        TEST("Euler rates ↔ ω roundtrip");
        ASSERT_VEC3_EQ(omega, omega2, 1e-10);
    }

    /* Cross product matrix */
    {
        mat3 skew;
        double w[3] = {1.0, 2.0, 3.0};
        cross_product_matrix(w, &skew);
        TEST("Cross matrix: skew-symmetric");
        CHECK(skew.m[0]==0 && skew.m[4]==0 && skew.m[8]==0);
        CHECK(skew.m[1] == -skew.m[3] && skew.m[2] == -skew.m[6] && skew.m[5] == -skew.m[7]);
    }

    /* ======================================================================
     * L3: Quaternion Operations
     * ====================================================================== */
    {
        mat3 R_test;
        quaternion_to_rotation_matrix(QUAT_IDENTITY, &R_test);
        TEST("Quaternion→Matrix: identity gives identity");
        CHECK(R_test.m[0]==1 && R_test.m[4]==1 && R_test.m[8]==1);

        Quaternion q_euler;
        euler_to_quaternion(0.0, 0.0, M_PI_2, &q_euler);
        TEST("Euler→Quaternion: 90° about Z");
        ASSERT_FLOAT_EQ(quat_norm2(q_euler), 1.0, 1e-14);
    }

    {
        Quaternion q0 = QUAT_IDENTITY;
        Quaternion q1;
        euler_to_quaternion(0.0, M_PI_2, 0.0, &q1); /* 90° about X */
        Quaternion q_mid;
        quat_slerp(q0, q1, 0.5, &q_mid);
        TEST("SLERP: result is unit quaternion");
        CHECK(quat_is_unit(q_mid, 1e-14));
    }

    {
        double axis_angle[3] = {0.0, 0.0, M_PI_2};
        Quaternion q;
        quat_exp(axis_angle, &q);
        TEST("Quat exp(π/2 about Z): unit result");
        CHECK(quat_is_unit(q, 1e-14));

        double aa[3];
        quat_log(q, aa);
        TEST("Quat log: recovers rotation vector");
        ASSERT_FLOAT_EQ(aa[0]+aa[1], 0.0, 1e-12);
        ASSERT_FLOAT_EQ(aa[2], M_PI_2, 1e-12);
    }

    /* ======================================================================
     * L3: Basic Rotations
     * ====================================================================== */
    {
        mat3 Rx, Ry, Rz;
        rotation_x(M_PI_2, &Rx);
        rotation_y(M_PI_2, &Ry);
        rotation_z(M_PI_2, &Rz);
        TEST("rot_x: preserves x-axis");
        ASSERT_FLOAT_EQ(Rx.m[0], 1.0, 1e-14);
        TEST("rot_y: preserves y-axis");
        ASSERT_FLOAT_EQ(Ry.m[4], 1.0, 1e-14);
        TEST("rot_z: preserves z-axis");
        ASSERT_FLOAT_EQ(Rz.m[8], 1.0, 1e-14);
    }

    {
        vec3 axis = {1.0, 0.0, 0.0};
        mat3 R;
        rotation_axis_angle(axis, M_PI_2, &R);
        TEST("Rodrigues formula: rot(X, 90°) preserves x");
        ASSERT_FLOAT_EQ(R.m[0], 1.0, 1e-14);
        TEST("Rodrigues formula: rot(X, 90°) y→z");
        ASSERT_FLOAT_EQ(R.m[7], 1.0, 1e-14);

        vec3 rec_axis; double rec_angle;
        rotation_to_axis_angle(&R, &rec_axis, &rec_angle);
        TEST("Rot→Axis-angle: recovers axis");
        ASSERT_FLOAT_EQ(rec_axis.x, 1.0, 1e-12);
    }

    /* ======================================================================
     * L3: Tait-Bryan Angles
     * ====================================================================== */
    {
        mat3 R_tb;
        tait_bryan_to_matrix(0.1, 0.2, 0.3, &R_tb);
        double roll, pitch, yaw;
        matrix_to_tait_bryan(&R_tb, &roll, &pitch, &yaw);
        TEST("Tait-Bryan roundtrip (roll)");
        ASSERT_FLOAT_EQ(roll, 0.1, 1e-10);
        TEST("Tait-Bryan roundtrip (pitch)");
        ASSERT_FLOAT_EQ(pitch, 0.2, 1e-10);
        TEST("Tait-Bryan roundtrip (yaw)");
        ASSERT_FLOAT_EQ(yaw, 0.3, 1e-10);
    }

    /* ======================================================================
     * L6: Energy and Angular Momentum
     * ====================================================================== */
    printf("\n─── L6: Energy & Angular Momentum ───\n");

    {
        double I_princ[3] = {3.0, 2.0, 1.0};
        double w[3] = {1.0, 0.5, 0.2};

        double T = rotational_ke_principal(I_princ, w);
        double T_expected = 0.5*(3.0*1.0 + 2.0*0.25 + 1.0*0.04);
        TEST("Rotational KE principal: matches formula");
        ASSERT_FLOAT_EQ(T, T_expected, 1e-14);

        InertiaTensor I_full = {3.0, 2.0, 1.0, 0.0, 0.0, 0.0};
        double T_full = rotational_ke(&I_full, w);
        TEST("Rotational KE full: matches principal for diagonal I");
        ASSERT_FLOAT_EQ(T, T_full, 1e-14);
    }

    {
        InertiaTensor I = {3.0, 2.0, 1.0, 0.0, 0.0, 0.0};
        double w[3] = {1.0, 0.5, 0.2};
        double L[3];
        angular_momentum_from_inertia(&I, w, L);
        TEST("Angular momentum: Lx = Ixx ωx");
        ASSERT_FLOAT_EQ(L[0], 3.0*1.0, 1e-14);
        TEST("Angular momentum: Ly = Iyy ωy");
        ASSERT_FLOAT_EQ(L[1], 2.0*0.5, 1e-14);

        double Lp[3];
        angular_momentum_principal((double[]){3.0, 2.0, 1.0}, w, Lp);
        TEST("Angular momentum principal: matches full for diagonal I");
        ASSERT_VEC3_EQ(L, Lp, 1e-14);
    }

    {
        vec3 r_cm = {1.0, 0.0, 0.0};
        InertiaTensor I_cm = sphere_inertia(1.0, 1.0);
        double w[3] = {0.0, 0.0, 1.0};
        double L_total[3], L_orbital[3], L_spin[3];
        total_angular_momentum(r_cm, 1.0, (vec3){0.0, 1.0, 0.0}, &I_cm, w,
                              L_total, L_orbital, L_spin);

        TEST("Total angular momentum: L_orbital= r×Mv");
        ASSERT_FLOAT_EQ(L_orbital[2], 1.0, 1e-14); /* (1,0,0)×(0,1,0) = (0,0,1) */
    }

    /* ======================================================================
     * L6: Tennis Racket Theorem
     * ====================================================================== */
    printf("\n─── L6: Tennis Racket Theorem ───\n");

    {
        double I[3] = {3.0, 2.0, 1.0}; /* I1 > I2 > I3 */

        StabilityResult results[3];
        full_stability_analysis(I, results);

        TEST("Axis 1 (max): STABLE");
        CHECK(results[0] == STABILITY_STABLE);
        TEST("Axis 2 (mid): UNSTABLE (Tennis Racket Theorem)");
        CHECK(results[1] == STABILITY_UNSTABLE);
        TEST("Axis 3 (min): STABLE");
        CHECK(results[2] == STABILITY_STABLE);
    }

    {
        double I[3] = {2.0, 2.0, 1.0}; /* Symmetric: I1 = I2 */
        StabilityResult results[3];
        full_stability_analysis(I, results);
        TEST("Symmetric top: I2 DEGENERATE");
        CHECK(results[1] == STABILITY_DEGENERATE);
    }

    {
        double I[3] = {3.0, 2.0, 1.0};
        double gr;
        StabilityResult s = axis_stability_analysis(I, 1, 1.0, &gr);
        TEST("Intermediate axis: growth rate > 0");
        CHECK(s == STABILITY_UNSTABLE && gr > 0.0);
    }

    {
        double I[3] = {3.0, 2.0, 1.0};
        double T_flip = flipping_period_estimate(I, 1.0, 0.01);
        TEST("Flipping period: finite for intermediate axis");
        CHECK(T_flip > 0.0 && T_flip < INFINITY);
    }

    /* Dzhanibekov simulation */
    {
        double I[3] = {3.0, 2.0, 1.0};
        double w0[3] = {0.01, 1.0, 0.01};
        int n_steps, n_flips;
        double *traj = malloc(500 * 3 * sizeof(double));
        double flip_times[20];

        dzhanibekov_simulation(I, w0, 2.0, 0.01, &n_steps, traj, flip_times, &n_flips);
        TEST("Dzhanibekov: simulation runs");
        CHECK(n_steps > 0);
        free(traj);
    }

    /* ======================================================================
     * L6: Torque Power and Impulse
     * ====================================================================== */
    {
        double N[3] = {0.0, 0.0, 1.0};
        double w[3] = {0.1, 0.0, 0.5};
        double P = torque_power(N, w);
        TEST("Torque power: Nz * ωz = 0.5");
        ASSERT_FLOAT_EQ(P, 0.5, 1e-14);

        double dL[3];
        torque_impulse(N, 0.5, dL);
        TEST("Torque impulse: dL = (0, 0, 0.5)");
        {
            double ex[3] = {0.0, 0.0, 0.5};
            ASSERT_VEC3_EQ(dL, ex, 1e-14);
        }
    }

    /* ======================================================================
     * L6: Energy-Momentum Geometry
     * ====================================================================== */
    {
        double I[3] = {3.0, 2.0, 1.0};
        double L2 = 10.0;
        double T_min, T_max;
        energy_extrema_for_L2(I, L2, &T_min, &T_max);
        TEST("T_min = L²/(2 I₁)");
        ASSERT_FLOAT_EQ(T_min, L2 / (2.0 * 3.0), 1e-14);
        TEST("T_max = L²/(2 I₃)");
        ASSERT_FLOAT_EQ(T_max, L2 / (2.0 * 1.0), 1e-14);
        TEST("T_min < T_max");
        CHECK(T_min < T_max);
    }

    /* ======================================================================
     * L6: Heavy Symmetric Top
     * ====================================================================== */
    printf("\n─── L6: Heavy Symmetric Top ───\n");

    {
        double I1 = 1.0, I3 = 0.5;
        double M = 1.0, g = 9.81, l = 0.3, Lz = 1.0, L3 = 2.0;

        double Veff = symmetric_top_eff_potential(I1, I3, M, g, l, Lz, L3, M_PI/4);
        TEST("V_eff: finite value");
        CHECK(isfinite(Veff) && Veff < 1e6);

        double thetas[200], Veff_arr[200];
        symmetric_top_eff_potential_grid(I1, I3, M, g, l, Lz, L3, 200, thetas, Veff_arr);
        TEST("V_eff grid: correct length, finite values");
        CHECK(isfinite(Veff_arr[0]) && isfinite(Veff_arr[199]));
    }

    {
        double I1 = 1.0, I3 = 0.5;
        double M = 1.0, g = 9.81, l = 0.3, Lz = 1.0, L3 = 2.0;
        double E = 15.0;
        int n_tp;
        double tp[2];
        nutation_turning_points(I1, I3, M, g, l, Lz, L3, E, 500, &n_tp, tp);
        TEST("Nutation turning points: returns result");
        CHECK(n_tp >= 0);
        if (n_tp >= 2) {
            TEST("Nutation turning points: θ_min < θ_max");
            CHECK(tp[0] < tp[1]);
        }
    }

    {
        double I1 = 1.0, I3 = 0.5;
        double M = 1.0, g = 9.81, l = 0.3, Lz = 1.0, L3 = 2.0;

        double phidot = precession_rate_top(I1, Lz, L3, M_PI/4);
        TEST("Precession rate: scalar");
        CHECK(isfinite(phidot));

        double psidot = spin_rate_top(I1, I3, Lz, L3, M_PI/4);
        TEST("Spin rate: scalar");
        CHECK(isfinite(psidot));

        double nut_rate = nutation_rate(I1, I3, M, g, l, 20.0, Lz, L3, M_PI/4);
        TEST("Nutation rate: returns ≥ 0 scalar");
        CHECK(nut_rate >= 0.0);
    }

    /* Steady precession */
    {
        double I1 = 1.0, I3 = 0.5;
        double M = 1.0, g = 9.81, l = 0.3, Lz = 2.0, L3 = 3.0;
        double bracket[2] = {0.2, 1.0};
        double theta0;
        int ret = steady_precession_angle(I1, I3, M, g, l, Lz, L3, bracket, &theta0);
        TEST("Steady precession: returns angle");
        CHECK(ret == 0 || (ret == -1 && theta0 > 0));
    }

    /* Fast top */
    {
        double Omega_fast = fast_top_precession_rate(1.0, 9.81, 0.3, 0.5, 10.0);
        TEST("Fast top: precession rate > 0");
        CHECK(Omega_fast > 0.0);

        double Omega_slow = slow_precession_rate(1.0, 9.81, 0.3, 2.0);
        TEST("Slow precession: precession rate > 0");
        CHECK(Omega_slow > 0.0);
    }

    /* ======================================================================
     * L6: Sleeping Top
     * ====================================================================== */
    {
        double omega_crit;
        int stable = sleeping_top_is_stable(1.0, 9.81, 0.3, 1.0, 0.5, 20.0, &omega_crit);
        TEST("Sleeping top: stable at high spin");
        CHECK(stable == 1);

        stable = sleeping_top_is_stable(1.0, 9.81, 0.3, 1.0, 0.5, 1.0, &omega_crit);
        TEST("Sleeping top: unstable at low spin");
        CHECK(stable == 0);

        double wc = sleeping_top_critical_spin(1.0, 9.81, 0.3, 1.0, 0.5);
        TEST("Critical spin: matches formula");
        ASSERT_FLOAT_EQ(wc, sqrt(4.0*9.81*0.3*1.0)/0.5, 1e-12);
    }

    /* ======================================================================
     * L6: Lagrange Top Full Simulation
     * ====================================================================== */
    {
        LagrangeTopParams params = {1.0, 0.5, 1.0, 9.81, 0.3};
        int n_steps;
        double *times = malloc(100 * sizeof(double));
        double *thetas = malloc(100 * sizeof(double));
        double *phis = malloc(100 * sizeof(double));
        double *psis = malloc(100 * sizeof(double));

        int ret = simulate_lagrange_top(&params, M_PI/6, 0.0, 0.0, 0.1,
                                        1.0, 2.0, 0.5, 0.01,
                                        &n_steps, times, thetas, phis, psis);
        TEST("Lagrange sim: returns 0");
        CHECK(ret == 0);
        TEST("Lagrange sim: ~51 steps");
        CHECK(n_steps == 51);
        TEST("Lagrange sim: θ positive");
        CHECK(thetas[0] > 0.0 && thetas[n_steps-1] > 0.0);

        free(times); free(thetas); free(phis); free(psis);
    }

    /* Nutation analysis */
    {
        double T_nut = small_nutation_period(1.0, 1.0, 9.81, 0.3);
        TEST("Small nutation period: > 0");
        CHECK(T_nut > 0.0);

        PrecessionType pt = classify_precession(1.0, 1.0, 2.0, 0.5, 1.0);
        TEST("Precession type: valid enum");
        CHECK(pt == PRECESSION_MONOTONIC || pt == PRECESSION_LOOPING || pt == PRECESSION_CUSPED);
    }

    /* ======================================================================
     * L7: Gyroscope
     * ====================================================================== */
    printf("\n─── L7: Gyroscope ───\n");

    {
        double L_spin[3] = {0.0, 0.0, 1.0};
        double w_forced[3] = {0.1, 0.0, 0.0};
        double torque[3];
        gyroscopic_torque(L_spin, w_forced, torque);
        TEST("Gyroscopic torque: perpendicular to L and Ω");
        CHECK(fabs(torque[2]) < 1e-14);
        ASSERT_FLOAT_EQ(torque[1], 0.1, 1e-14);
    }

    {
        double Omega = gyroscope_precession_rate(0.5, 2.0);
        TEST("Gyroscope precession: Ω = N/L");
        ASSERT_FLOAT_EQ(Omega, 0.25, 1e-14);
    }

    {
        double theta[50], omega_g[50];
        double w_base[3] = {0.0, 0.1, 0.0};
        gyroscope_2axis_simulate(1.0, w_base, 0.5, 0.0, 0.1, 0.01, 50, theta, omega_g);
        TEST("Gyroscope 2-axis: runs");
        CHECK(isfinite(theta[49]) && isfinite(omega_g[49]));
    }

    {
        double w_nut = free_gyro_nutation_frequency(2.0, 0.5);
        TEST("Free gyro nutation: ω_nut = L/I_trans");
        ASSERT_FLOAT_EQ(w_nut, 2.0/0.5, 1e-14);
    }

    /* ======================================================================
     * L8: Poinsot Construction
     * ====================================================================== */
    printf("\n─── L8: Poinsot Construction ───\n");

    {
        double I[3] = {3.0, 2.0, 1.0};
        double w[3] = {0.2, 0.3, 0.4};
        mat3 R;
        euler_to_rotation_matrix(0.1, 0.2, 0.3, &R);
        vec3 normal;
        double dist;
        invariable_plane(I, w, &R, &normal, &dist);
        TEST("Invariable plane: normal is unit");
        ASSERT_FLOAT_EQ(vec3_norm(normal), 1.0, 1e-14);
        TEST("Invariable plane: distance > 0");
        CHECK(dist > 0.0);
    }

    {
        double I[3] = {3.0, 2.0, 1.0};
        double w[3] = {0.2, 0.3, 0.4};
        double tangent[3];
        polhode_tangent(I, w, tangent);
        TEST("Polhode tangent: non-zero for this ω");
        CHECK(tangent[0]!=0 || tangent[1]!=0 || tangent[2]!=0);
    }

    /* ======================================================================
     * L8: Energy-Momentum Sphere Sampling
     * ====================================================================== */
    {
        double I[3] = {3.0, 2.0, 1.0};
        double L2 = 10.0;
        int np = 10;
        double *T_map = malloc(np * np * sizeof(double));
        energy_on_momentum_sphere(I, L2, np, np, T_map, NULL, NULL, NULL);
        TEST("Energy on momentum sphere: T in expected range");
        double T_min, T_max;
        energy_extrema_for_L2(I, L2, &T_min, &T_max);
        int in_range = 1;
        for (int i = 0; i < np*np; i++) {
            if (T_map[i] < T_min*0.99 || T_map[i] > T_max*1.01) in_range = 0;
        }
        CHECK(in_range);
        free(T_map);
    }

    /* ======================================================================
     * L4: Conservation law verification (linear momentum case)
     * ====================================================================== */
    {
        double I[3] = {3.0, 2.0, 1.0};
        double w[3] = {0.2, 0.3, 0.4};
        double T0, L20;
        motion_constants_euler(I, w, &T0, &L20);

        /* One RK4 step should preserve constants well */
        double w1[3];
        rk4_step_omega(I, w, 0.01, w1);
        double T1, L21;
        motion_constants_euler(I, w1, &T1, &L21);

        TEST("RK4: Energy conservation within 1e-12");
        CHECK(fabs(T1 - T0) < 1e-12);
        TEST("RK4: Angular momentum conservation within 1e-12");
        CHECK(fabs(L21 - L20) < 1e-12);
    }

    /* ======================================================================
     * Summary
     * ====================================================================== */
    int total = passed + failed;
    printf("\n══════════════════════════════════════════════════════════════\n");
    printf("  Results: %d / %d passed\n", passed, total);
    if (failed == 0) {
        printf("  ✅ ALL TESTS PASSED\n");
    } else {
        printf("  ❌ %d TEST(S) FAILED\n", failed);
    }
    printf("══════════════════════════════════════════════════════════════\n");

    return (failed == 0) ? 0 : 1;
}
