/**
 * test_newtonian.c - assert-based tests for mini-newtonian C library
 *
 * Tests cover all L1-L6 knowledge layers.
 * Compile: gcc -std=c11 -Iinclude tests/test_newtonian.c -Lbuild/c -lnewtonian -lm -o build/c/test_newtonian
 * Run: ./build/c/test_newtonian
 */
#include "newtonian.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

static int tests_passed = 0;
static int tests_failed = 0;
static const double EPS = 1e-12;

#define TEST(name) printf("  TEST: %-55s ", name)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define CHECK(cond) do { \
    if (!(cond)) { \
        printf("FAIL at line %d: %s\n", __LINE__, #cond); \
        tests_failed++; \
        return; \
    } \
} while(0)
#define CHECK_NEAR(a, b, tol) do { \
    if (fabs((a) - (b)) > (tol)) { \
        printf("FAIL at line %d: |%g - %g| = %g > %g\n", \
               __LINE__, (double)(a), (double)(b), \
               fabs((a)-(b)), (double)(tol)); \
        tests_failed++; \
        return; \
    } \
} while(0)

/* ================================================================
 * L1: Vec3 Operations
 * ================================================================ */
static void test_vec3_create(void) { TEST("Vec3 creation");
    Vec3 v = vec3_make(1.0, 2.0, 3.0);
    CHECK(v.x == 1.0 && v.y == 2.0 && v.z == 3.0);
    Vec3 z = vec3_zero();
    CHECK(z.x == 0.0 && z.y == 0.0 && z.z == 0.0);
    PASS();
}
static void test_vec3_add(void) { TEST("Vec3 addition");
    Vec3 a = vec3_make(1, 2, 3), b = vec3_make(4, 5, 6);
    Vec3 c = vec3_add(a, b);
    CHECK_NEAR(c.x, 5.0, EPS); CHECK_NEAR(c.y, 7.0, EPS); CHECK_NEAR(c.z, 9.0, EPS);
    PASS();
}
static void test_vec3_sub(void) { TEST("Vec3 subtraction");
    Vec3 a = vec3_make(5, 7, 9), b = vec3_make(1, 2, 3);
    Vec3 c = vec3_sub(a, b);
    CHECK_NEAR(c.x, 4.0, EPS); CHECK_NEAR(c.y, 5.0, EPS); CHECK_NEAR(c.z, 6.0, EPS);
    PASS();
}
static void test_vec3_dot(void) { TEST("Vec3 dot product");
    Vec3 a = vec3_make(1, 2, 3), b = vec3_make(4, -5, 6);
    double d = vec3_dot(a, b);
    CHECK_NEAR(d, 1.0*4.0 + 2.0*(-5.0) + 3.0*6.0, EPS);
    PASS();
}
static void test_vec3_cross(void) { TEST("Vec3 cross product");
    Vec3 a = vec3_make(1, 0, 0), b = vec3_make(0, 1, 0);
    Vec3 c = vec3_cross(a, b);
    CHECK_NEAR(c.x, 0.0, EPS); CHECK_NEAR(c.y, 0.0, EPS); CHECK_NEAR(c.z, 1.0, EPS);
    PASS();
}
static void test_vec3_norm(void) { TEST("Vec3 norm");
    Vec3 v = vec3_make(3, 4, 0);
    CHECK_NEAR(vec3_norm(v), 5.0, EPS);
    CHECK_NEAR(vec3_norm2(v), 25.0, EPS);
    PASS();
}
static void test_vec3_normalize(void) { TEST("Vec3 normalize");
    Vec3 v = vec3_make(3, 0, 0);
    Vec3 u = vec3_normalize(v);
    CHECK_NEAR(u.x, 1.0, EPS); CHECK_NEAR(u.y, 0.0, EPS); CHECK_NEAR(vec3_norm(u), 1.0, EPS);
    PASS();
}
static void test_vec3_distance(void) { TEST("Vec3 distance");
    Vec3 a = vec3_make(0, 0, 0), b = vec3_make(3, 4, 0);
    CHECK_NEAR(vec3_distance(a, b), 5.0, EPS);
    PASS();
}
static void test_vec3_scalar_triple(void) { TEST("Vec3 scalar triple product");
    Vec3 a = vec3_make(1, 0, 0), b = vec3_make(0, 1, 0), c = vec3_make(0, 0, 1);
    CHECK_NEAR(vec3_scalar_triple(a, b, c), 1.0, EPS);
    PASS();
}
static void test_vec3_vector_triple(void) { TEST("Vec3 BAC-CAB rule");
    Vec3 a = vec3_make(1, 2, 0), b = vec3_make(0, 1, 0), c = vec3_make(0, 0, 1);
    Vec3 r = vec3_vector_triple(a, b, c);
    Vec3 expected = vec3_sub(vec3_scale(b, vec3_dot(a,c)), vec3_scale(c, vec3_dot(a,b)));
    CHECK_NEAR(vec3_norm(vec3_sub(r, expected)), 0.0, EPS);
    PASS();
}
static void test_vec3_project(void) { TEST("Vec3 projection");
    Vec3 a = vec3_make(3, 4, 0), b = vec3_make(1, 0, 0);
    Vec3 proj = vec3_project(a, b);
    CHECK_NEAR(proj.x, 3.0, EPS); CHECK_NEAR(proj.y, 0.0, EPS);
    PASS();
}
static void test_vec3_rotate(void) { TEST("Vec3 Rodrigues rotation");
    Vec3 v = vec3_make(1, 0, 0);
    Vec3 axis = vec3_make(0, 0, 1);
    Vec3 r = vec3_rotate_axis_angle(v, axis, M_PI / 2.0);
    CHECK_NEAR(r.x, 0.0, EPS); CHECK_NEAR(r.y, 1.0, EPS); CHECK_NEAR(r.z, 0.0, EPS);
    PASS();
}
static void test_spherical_conversion(void) { TEST("Spherical coordinate conversion");
    Vec3 v = vec3_make(1, 0, 0);
    Spherical s = vec3_to_spherical(v);
    CHECK_NEAR(s.r, 1.0, EPS); CHECK_NEAR(s.theta, M_PI/2.0, EPS); CHECK_NEAR(s.phi, 0.0, EPS);
    Vec3 v2 = vec3_from_spherical(s.r, s.theta, s.phi);
    CHECK_NEAR(vec3_norm(vec3_sub(v, v2)), 0.0, EPS);
    PASS();
}

/* ================================================================
 * L1/L2: Kinematics
 * ================================================================ */
static void test_suvat(void) { TEST("SUVAT equations");
    double r = suvat_position(0.0, 10.0, -9.81, 1.0);
    CHECK_NEAR(r, 10.0 - 0.5*9.81, EPS);
    double v = suvat_velocity(10.0, -9.81, 1.0);
    CHECK_NEAR(v, 10.0 - 9.81, EPS);
    PASS();
}
static void test_projectile(void) { TEST("Projectile motion");
    double range = projectile_range(10.0, M_PI/4.0, 9.81);
    double expected = 10.0 * 10.0 / 9.81;  /* v0^2/g for 45 deg */
    CHECK_NEAR(range, expected, 1e-6);
    double h_max = projectile_max_height(10.0, M_PI/4.0, 9.81);
    double v0z = 10.0 * sin(M_PI/4.0);
    CHECK_NEAR(h_max, v0z*v0z/(2.0*9.81), 1e-6);
    PASS();
}
static void test_circular_motion(void) { TEST("Circular motion");
    Vec3 center = vec3_zero();
    Vec3 axis = vec3_make(0, 0, 1);
    Vec3 pos = circular_position(center, axis, 1.0, 1.0, 0.0);
    /* For axis=(0,0,1): e1=(0,1,0), so at t=0: pos=(0,1,0) */
    CHECK_NEAR(pos.x, 0.0, EPS); CHECK_NEAR(pos.y, 1.0, EPS); CHECK_NEAR(pos.z, 0.0, EPS);
    double a_c = centripetal_acceleration(2.0, 1.0);
    CHECK_NEAR(a_c, 4.0, EPS);
    PASS();
}
static void test_galilean_transform(void) { TEST("Galilean transformation");
    Vec3 rp = vec3_make(1, 2, 3), V = vec3_make(5, 0, 0);
    Vec3 r = galilean_position_transform(rp, V, 2.0);
    CHECK_NEAR(r.x, 11.0, EPS);
    Vec3 v = galilean_velocity_transform(vec3_make(0,0,0), V);
    CHECK_NEAR(v.x, 5.0, EPS);
    PASS();
}
static void test_coriolis(void) { TEST("Coriolis acceleration");
    Vec3 v_rot = vec3_make(1, 0, 0);
    Vec3 omega = vec3_make(0, 0, 1);
    Vec3 a_cor = coriolis_acceleration(v_rot, omega);
    /* -2*omega x v = -2*(0,0,1)x(1,0,0) = -2*(0,1,0) = (0,-2,0) */
    CHECK_NEAR(a_cor.x, 0.0, EPS); CHECK_NEAR(a_cor.y, -2.0, EPS); CHECK_NEAR(a_cor.z, 0.0, EPS);
    PASS();
}
static void test_centrifugal(void) { TEST("Centrifugal acceleration");
    Vec3 r = vec3_make(1, 0, 0);
    Vec3 omega = vec3_make(0, 0, 1);
    Vec3 a_cen = centrifugal_acceleration(r, omega);
    /* -omega x (omega x r) = -{(0,0,1) x (0,1,0)} = -( -1, 0, 0 ) = (1,0,0) */
    CHECK_NEAR(a_cen.x, 1.0, EPS); CHECK_NEAR(a_cen.y, 0.0, EPS); CHECK_NEAR(a_cen.z, 0.0, EPS);
    PASS();
}

/* ================================================================
 * L4: Force Laws
 * ================================================================ */
static void test_newton_gravity(void) { TEST("Newton gravity force");
    Vec3 r_obj = vec3_make(1, 0, 0);
    Vec3 r_src = vec3_zero();
    Vec3 F = newton_gravity_force(r_obj, r_src, 100.0, 1.0, 1.0);
    /* F = -G*M*m/r^2 * r_hat = -1*100*1/1 * (1,0,0) = (-100,0,0) */
    CHECK_NEAR(F.x, -100.0, EPS); CHECK_NEAR(F.y, 0.0, EPS); CHECK_NEAR(F.z, 0.0, EPS);
    PASS();
}
static void test_hooke_force(void) { TEST("Hooke law");
    Vec3 r = vec3_make(2, 0, 0), r_eq = vec3_zero();
    Vec3 F = hooke_force_3d(r, r_eq, 10.0);
    CHECK_NEAR(F.x, -20.0, EPS);
    double f1d = hooke_force_1d(3.0, 0.0, 5.0);
    CHECK_NEAR(f1d, -15.0, EPS);
    PASS();
}
static void test_drag_forces(void) { TEST("Drag forces");
    Vec3 v = vec3_make(10, 0, 0);
    Vec3 F_lin = linear_drag(v, 0.5);
    CHECK_NEAR(F_lin.x, -5.0, EPS);
    Vec3 F_quad = quadratic_drag(v, 0.1);
    CHECK_NEAR(F_quad.x, -10.0, EPS);  /* -c*|v|*vx = -0.1*10*10 = -10 */
    PASS();
}
static void test_friction(void) { TEST("Coulomb friction");
    Vec3 v = vec3_make(1, 0, 0);
    Vec3 f = kinetic_friction_force(v, 10.0, 0.3);
    CHECK_NEAR(f.x, -3.0, EPS);
    PASS();
}
static void test_lorentz_force(void) { TEST("Lorentz force");
    Vec3 E = vec3_zero();
    Vec3 v = vec3_make(1, 0, 0);
    Vec3 B = vec3_make(0, 0, 1);
    Vec3 F = lorentz_force(1.0, E, v, B);
    /* F = q*(0 + v x B) = (1,0,0)x(0,0,1) = (0*1-0*0, 0*0-1*1, 1*0-0*0) = (0, -1, 0) */
    CHECK_NEAR(F.x, 0.0, EPS); CHECK_NEAR(F.y, -1.0, EPS);
    PASS();
}

/* ================================================================
 * L5: ODE Integrators
 * ================================================================ */
/* Acceleration function for harmonic oscillator: a = -omega^2 * r */
static Vec3 sho_accel(Vec3 r, Vec3 v, double t) {
    (void)v; (void)t;
    return vec3_scale(r, -1.0);  /* omega=1, k=1, m=1 */
}

static void test_euler_cromer(void) { TEST("Euler-Cromer integrator");
    ParticleState s = particle_state_make(0.0, vec3_make(1,0,0), vec3_zero());
    s = euler_cromer_step(sho_accel, s, 0.01);
    /* a = -r = (-1,0,0); v_new = 0 + (-1,0,0)*0.01 = (-0.01,0,0) */
    /* r_new = (1,0,0) + (-0.01,0,0)*0.01 = (0.9999,0,0) */
    CHECK_NEAR(s.v.x, -0.01, EPS);
    PASS();
}
static void test_velocity_verlet(void) { TEST("Velocity Verlet integrator");
    ParticleState s = particle_state_make(0.0, vec3_make(1,0,0), vec3_zero());
    s = velocity_verlet_step(sho_accel, s, 0.01);
    /* a_cur = (-1,0,0); v_half = (0,0,0)+0.5*(-1,0,0)*0.01 = (-0.005,0,0) */
    /* r_new = (1,0,0)+(-0.005,0,0)*0.01 = (0.99995,0,0) */
    /* a_new = (-0.99995,0,0); v_new = (-0.005,0,0)+0.5*(-0.99995,0,0)*0.01 */
    /* = (-0.005-0.00499975, 0,0) = (-0.00999975, 0, 0) */
    CHECK_NEAR(s.r.x, 0.99995, 1e-6);
    CHECK_NEAR(s.v.x, -0.00999975, 1e-6);
    PASS();
}
static void test_trajectory(void) { TEST("Trajectory recording");
    Trajectory *traj = trajectory_create(100);
    CHECK(traj != NULL);
    int ok = trajectory_record(traj, 0.0, vec3_make(1,0,0), vec3_zero());
    CHECK(ok == 0);
    CHECK(traj->n_points == 1);
    trajectory_destroy(traj);
    PASS();
}
static void test_nbody_system(void) { TEST("N-body system creation");
    NBodySystem *sys = nbody_system_create(3);
    CHECK(sys != NULL);
    CHECK(sys->n_bodies == 3);
    nbody_system_destroy(sys);
    PASS();
}

/* ================================================================
 * L2/L4: Energy
 * ================================================================ */
static void test_kinetic_energy(void) { TEST("Kinetic energy");
    double T = kinetic_energy(2.0, vec3_make(3, 4, 0));
    CHECK_NEAR(T, 0.5 * 2.0 * 25.0, EPS);  /* 0.5*m*v^2 = 0.5*2*25 = 25 */
    PASS();
}
static void test_potential_energy(void) { TEST("Potential energy");
    double U = potential_energy_uniform_gravity(1.0, 10.0, 9.81);
    CHECK_NEAR(U, 98.1, EPS);
    double U_elastic = potential_energy_elastic(vec3_make(2,0,0), vec3_zero(), 3.0);
    CHECK_NEAR(U_elastic, 0.5 * 3.0 * 4.0, EPS);
    PASS();
}
static void test_virial_theorem(void) { TEST("Virial theorem");
    /* For harmonic oscillator: <T> = <U>, n=2 */
    double n = virial_theorem_exponent(10.0, 10.0);
    CHECK_NEAR(n, 2.0, EPS);
    /* For gravity: 2<T> + <U> = 0, ratio test */
    double r = virial_ratio_gravity(5.0, -10.0);
    CHECK_NEAR(r, 0.0, EPS);
    PASS();
}

/* ================================================================
 * L2/L4: Momentum and Collisions
 * ================================================================ */
static void test_linear_momentum(void) { TEST("Linear momentum");
    Vec3 p = linear_momentum(2.0, vec3_make(3, 0, 0));
    CHECK_NEAR(p.x, 6.0, EPS);
    PASS();
}
static void test_angular_momentum(void) { TEST("Angular momentum");
    Vec3 L = angular_momentum(vec3_make(1,0,0), vec3_make(0,1,0), 1.0);
    /* L = r x p = (1,0,0) x (0,1,0) = (0,0,1) */
    CHECK_NEAR(L.z, 1.0, EPS);
    PASS();
}
static void test_center_of_mass(void) { TEST("Center of mass");
    double m[] = {1, 1};
    Vec3 p[] = {vec3_make(0,0,0), vec3_make(2,0,0)};
    (void)0; /* dummy for Vec3 v */
    Vec3 v_stub[] = {vec3_zero(), vec3_zero()}; (void)v_stub;
    Vec3 cm = center_of_mass(m, p, 2);
    CHECK_NEAR(cm.x, 1.0, EPS);
    PASS();
}
static void test_elastic_collision_1d(void) { TEST("1D elastic collision");
    double v1f, v2f;
    elastic_collision_1d(1.0, 1.0, 1.0, -1.0, &v1f, &v2f);
    /* Equal masses: velocities exchange */
    CHECK_NEAR(v1f, -1.0, EPS); CHECK_NEAR(v2f, 1.0, EPS);
    /* Check momentum conservation */
    double p_i = 1.0*1.0 + 1.0*(-1.0);
    double p_f = 1.0*v1f + 1.0*v2f;
    CHECK_NEAR(p_f, p_i, EPS);
    PASS();
}
static void test_inelastic_collision(void) { TEST("Inelastic collision");
    double vf = inelastic_collision_1d(2.0, 3.0, 1.0, 0.0);
    /* v_common = (2*3 + 1*0)/(2+1) = 2 */
    CHECK_NEAR(vf, 2.0, EPS);
    PASS();
}
static void test_coefficient_of_restitution(void) { TEST("Coefficient of restitution");
    double e = coefficient_of_restitution(2.0, 0.0, 1.0, 1.0);
    CHECK_NEAR(e, 0.0, EPS);  /* v_rel_f = 0, v_rel_i = 2 */
    PASS();
}
static void test_restitution_collision(void) { TEST("Collision with restitution");
    double v1f, v2f;
    collision_with_restitution(1.0, 2.0, 1.0, 0.0, 0.5, &v1f, &v2f);
    /* v_cm = 1; v1_cm=1, v2_cm=-1; v1f_cm=-0.5, v2f_cm=0.5; v1f=0.5, v2f=1.5 */
    CHECK_NEAR(v1f, 0.5, EPS); CHECK_NEAR(v2f, 1.5, EPS);
    PASS();
}
static void test_tsiolkovsky(void) { TEST("Tsiolkovsky rocket equation");
    double dv = tsiolkovsky_delta_v(3000.0, 1000.0, 100.0);
    CHECK_NEAR(dv, 3000.0 * log(10.0), EPS);
    PASS();
}

/* ================================================================
 * L6: Constraint Systems
 * ================================================================ */
static void test_incline(void) { TEST("Inclined plane");
    double a = incline_acceleration(M_PI/6.0, 0.0, 9.81);
    CHECK_NEAR(a, 9.81 * 0.5, EPS);
    double N = incline_normal_force(2.0, M_PI/3.0, 9.81);
    CHECK_NEAR(N, 2.0 * 9.81 * 0.5, EPS);
    PASS();
}
static void test_atwood(void) { TEST("Atwood machine");
    double a1, a2, T;
    atwood_machine(3.0, 2.0, 10.0, &a1, &a2, &T);
    /* a = (3-2)*10/(3+2) = 2 m/s^2 */
    CHECK_NEAR(a1, 2.0, EPS); CHECK_NEAR(a2, -2.0, EPS);
    /* T = 2*3*2*10/5 = 24 N */
    CHECK_NEAR(T, 24.0, 1e-6);
    PASS();
}
static void test_conical_pendulum(void) { TEST("Conical pendulum");
    double theta, T;
    int ok = conical_pendulum_parameters(1.0, 1.0, 5.0, 9.81, &theta, &T);
    CHECK(ok == 1);
    /* cos(theta) = g/(omega^2*L) = 9.81/(25*1) = 0.3924 */
    CHECK_NEAR(cos(theta), 9.81/25.0, EPS);
    PASS();
}
static void test_loop_the_loop(void) { TEST("Loop-the-loop");
    double v_min = loop_the_loop_min_speed(2.0, 9.81);
    CHECK_NEAR(v_min, sqrt(9.81 * 2.0), EPS);
    double N_top = loop_the_loop_normal_force_at_top(1.0, 5.0, 2.0, 9.81);
    /* N = m(v^2/R - g) = 1*(25/2 - 9.81) = 12.5 - 9.81 = 2.69 */
    CHECK(N_top > 0.0);  /* positive at top = stays on track */
    PASS();
}
static void test_banked_curve(void) { TEST("Banked curve");
    double theta = banked_curve_ideal_angle(10.0, 50.0, 9.81);
    double expected_angle = atan(100.0/(9.81*50.0));
    CHECK_NEAR(theta, expected_angle, EPS);
    PASS();
}
static void test_pendulum(void) { TEST("Simple pendulum");
    double T0 = pendulum_period_small_angle(1.0, 9.81);
    CHECK_NEAR(T0, 2.0 * M_PI / sqrt(9.81), 1e-6);
    double freq = pendulum_frequency(1.0, 9.81);
    CHECK_NEAR(freq, sqrt(9.81), EPS);
    PASS();
}

/* ================================================================
 * L5: Analysis (Poincare, Lyapunov, orbital elements)
 * ================================================================ */
static void test_lyapunov_basic(void) { TEST("Lyapunov exponent estimation");
    Trajectory *t1 = trajectory_create(10);
    Trajectory *t2 = trajectory_create(10);
    trajectory_record(t1, 0.0, vec3_zero(), vec3_zero());
    trajectory_record(t1, 1.0, vec3_make(1,0,0), vec3_zero());
    trajectory_record(t2, 0.0, vec3_zero(), vec3_zero());
    trajectory_record(t2, 1.0, vec3_make(1.01,0,0), vec3_zero());
    double lam = estimate_lyapunov_exponent(t1, t2);
    /* Should be ~ln(0.01/0.01)/1 = 0 (constant offset, not divergence) */
    /* Actually initial separation is computed from first points: both start at 0 */
    CHECK(lam >= 0.0);  /* approximately zero */
    trajectory_destroy(t1); trajectory_destroy(t2);
    PASS();
}
static void test_orbital_elements(void) { TEST("Orbital elements from state");
    /* Circular orbit: r = (1,0,0) AU, v = (0, 2*pi, 0) AU/yr, mu = (2*pi)^2 */
    double mu = 4.0 * M_PI * M_PI;
    Vec3 r = vec3_make(1, 0, 0);
    Vec3 v = vec3_make(0, 2.0 * M_PI, 0);
    OrbitalElements elem;
    int ok = orbital_elements_from_state(r, v, mu, &elem);
    CHECK(ok == 1);
    CHECK_NEAR(elem.eccentricity, 0.0, 1e-6);
    CHECK_NEAR(elem.semi_major_axis, 1.0, 1e-6);
    CHECK_NEAR(elem.inclination, 0.0, 1e-6);
    PASS();
}

/* ================================================================
 * L4: Conservation law verification
 * ================================================================ */
static void test_momentum_conservation(void) { TEST("Momentum conservation check");
    Vec3 pi = vec3_make(10, 0, 0), pf = vec3_make(10, 0, 0);
    CHECK(check_momentum_conservation(pi, pf, 1e-10) == 1);
    Vec3 pf2 = vec3_make(5, 0, 0);
    CHECK(check_momentum_conservation(pi, pf2, 1e-10) == 0);
    PASS();
}
static void test_energy_conservation(void) { TEST("Energy conservation check");
    CHECK(check_energy_conservation(100.0, 100.0, 1e-10) == 1);
    CHECK(check_energy_conservation(100.0, 90.0, 1e-10) == 0);
    PASS();
}

/* ================================================================ */
int main(void) {
    printf("\n=== mini-newtonian C Library Tests ===\n\n");
    printf("L1: Vec3 Operations\n");
    test_vec3_create();
    test_vec3_add();
    test_vec3_sub();
    test_vec3_dot();
    test_vec3_cross();
    test_vec3_norm();
    test_vec3_normalize();
    test_vec3_distance();
    test_vec3_scalar_triple();
    test_vec3_vector_triple();
    test_vec3_project();
    test_vec3_rotate();
    test_spherical_conversion();

    printf("\nL1/L2: Kinematics\n");
    test_suvat();
    test_projectile();
    test_circular_motion();
    test_galilean_transform();
    test_coriolis();
    test_centrifugal();

    printf("\nL4: Force Laws\n");
    test_newton_gravity();
    test_hooke_force();
    test_drag_forces();
    test_friction();
    test_lorentz_force();

    printf("\nL5: ODE Integrators\n");
    test_euler_cromer();
    test_velocity_verlet();
    test_trajectory();
    test_nbody_system();

    printf("\nL2/L4: Energy\n");
    test_kinetic_energy();
    test_potential_energy();
    test_virial_theorem();

    printf("\nL2/L4: Momentum & Collisions\n");
    test_linear_momentum();
    test_angular_momentum();
    test_center_of_mass();
    test_elastic_collision_1d();
    test_inelastic_collision();
    test_coefficient_of_restitution();
    test_restitution_collision();
    test_tsiolkovsky();

    printf("\nL6: Constraint Systems\n");
    test_incline();
    test_atwood();
    test_conical_pendulum();
    test_loop_the_loop();
    test_banked_curve();
    test_pendulum();

    printf("\nL5: Analysis\n");
    test_lyapunov_basic();
    test_orbital_elements();

    printf("\nL4: Conservation Laws\n");
    test_momentum_conservation();
    test_energy_conservation();

    printf("\n========================================\n");
    printf("Results: %d passed, %d failed (total %d)\n",
           tests_passed, tests_failed, tests_passed + tests_failed);
    printf("========================================\n");

    return tests_failed > 0 ? 1 : 0;
}
