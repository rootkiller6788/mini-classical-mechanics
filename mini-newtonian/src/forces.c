/**
 * forces.c - Newtonian force laws implementation
 * Each function implements one independent physics knowledge point.
 * References: Newton (1687), Hooke (1678), Stokes (1851), Coulomb (1785), Lorentz (1895)
 * Course: MIT 8.012 Ch.3-4, Goldstein Ch.1-2
 */
#include "forces.h"
#include <math.h>

/* ============ L4: Universal Gravitation (Newton 1687) ============ */
Vec3 newton_gravity_force(Vec3 r_obj, Vec3 r_source, double M, double m, double G) {
    Vec3 dr = vec3_sub(r_obj, r_source);
    double dist = vec3_norm(dr);
    if (dist < 1e-10) return vec3_zero();
    Vec3 r_hat = vec3_normalize(dr);
    return vec3_scale(r_hat, -G * M * m / (dist * dist));
}
Vec3 gravity_acceleration_field(Vec3 r_obj, Vec3 r_source, double M, double G) {
    Vec3 dr = vec3_sub(r_obj, r_source);
    double dist = vec3_norm(dr);
    if (dist < 1e-10) return vec3_zero();
    return vec3_scale(vec3_normalize(dr), -G * M / (dist * dist));
}
Vec3 gravity_nbody_force(Vec3 r_obj, double m, Vec3 sources_pos[],
                         double sources_mass[], int n_sources, double G) {
    Vec3 f_total = vec3_zero();
    for (int i = 0; i < n_sources; i++) {
        Vec3 fi = newton_gravity_force(r_obj, sources_pos[i], sources_mass[i], m, G);
        f_total = vec3_add(f_total, fi);
    }
    return f_total;
}
Vec3 uniform_gravity_force(double m, Vec3 direction, double g) {
    return vec3_scale(vec3_normalize(direction), m * g);
}
Vec3 gravity_potential_gradient(Vec3 r_obj, Vec3 r_source, double M, double G) {
    return gravity_acceleration_field(r_obj, r_source, M, G);
}

/* ============ L4: Hooke's Law (1660) ============ */
Vec3 hooke_force_3d(Vec3 r, Vec3 r_eq, double k) {
    return vec3_scale(vec3_sub(r, r_eq), -k);
}
double hooke_force_1d(double x, double x_eq, double k) {
    return -k * (x - x_eq);
}
double hooke_anharmonic_force_1d(double x, double x_eq, double k, double beta) {
    double dx = x - x_eq;
    return -k * dx - beta * dx * dx * dx;
}

/* ============ L4: Drag forces (Stokes 1851, Rayleigh) ============ */
Vec3 linear_drag(Vec3 v, double b) {
    return vec3_scale(v, -b);
}
Vec3 quadratic_drag(Vec3 v, double c) {
    double speed = vec3_norm(v);
    return vec3_scale(v, -c * speed);
}
Vec3 combined_drag(Vec3 v, double b, double c) {
    double speed = vec3_norm(v);
    return vec3_scale(v, -(b + c * speed));
}
double drag_reynolds_threshold(double rho_fluid, double eta,
                               double characteristic_length) {
    (void)rho_fluid; (void)eta; (void)characteristic_length;
    return 2000.0;
}

/* ============ L4: Coulomb/Amontons Friction ============ */
Vec3 kinetic_friction_force(Vec3 velocity, double normal_mag, double mu_k) {
    double speed = vec3_norm(velocity);
    if (speed < 1e-10) return vec3_zero();
    return vec3_scale(velocity, -mu_k * normal_mag / speed);
}
double static_friction_max(double normal_mag, double mu_s) {
    return mu_s * normal_mag;
}
double static_friction_threshold(double normal_mag, double mu_s,
                                 double applied_tangential) {
    double f_max = mu_s * normal_mag;
    if (applied_tangential <= f_max) return applied_tangential;
    return f_max;
}

/* ============ L4: Lorentz Force (1895) ============ */
Vec3 lorentz_force(double q, Vec3 E, Vec3 v, Vec3 B) {
    Vec3 v_cross_B = vec3_cross(v, B);
    return vec3_scale(vec3_add(E, v_cross_B), q);
}
double cyclotron_frequency(double q, double B_mag, double m) {
    return fabs(q) * B_mag / m;
}
double cyclotron_radius(double v_perp, double q, double B_mag, double m) {
    return m * v_perp / (fabs(q) * B_mag);
}
Vec3 magnetic_mirror_force(double mu_magnetic_moment, Vec3 B_gradient) {
    return vec3_scale(B_gradient, -mu_magnetic_moment);
}

/* ============ L2: Composite force functions ============ */
Vec3 damped_harmonic_force(Vec3 r, Vec3 r_eq, Vec3 v, double k, double b) {
    return vec3_add(hooke_force_3d(r, r_eq, k), linear_drag(v, b));
}
Vec3 driven_damped_force(Vec3 r, Vec3 r_eq, Vec3 v, double k, double b,
                         double F0, double omega_drive, double t) {
    Vec3 restoring = hooke_force_3d(r, r_eq, k);
    Vec3 damping = linear_drag(v, b);
    Vec3 driving = vec3_make(F0 * cos(omega_drive * t), 0.0, 0.0);
    return vec3_add(restoring, vec3_add(damping, driving));
}
Vec3 van_der_pol_force(Vec3 r, Vec3 r_eq, Vec3 v, double k, double mu) {
    Vec3 dx = vec3_sub(r, r_eq);
    Vec3 restoring = vec3_scale(dx, -k);
    double x2 = vec3_norm2(dx);
    Vec3 nonlinear = vec3_scale(v, mu * (1.0 - x2));
    return vec3_add(restoring, nonlinear);
}
Vec3 duffing_force_3d(Vec3 r, Vec3 r_eq, double k, double alpha, double beta) {
    Vec3 dx = vec3_sub(r, r_eq);
    double x = dx.x;
    double f_x = -k * x - alpha * x * x - beta * x * x * x;
    return vec3_make(f_x, 0.0, 0.0);
}

/* ============ L4/L7: Gravitational multipole expansion ============ */
Vec3 gravity_monopole_term(Vec3 r, double M, double G) {
    double dist = vec3_norm(r);
    if (dist < 1e-15) return vec3_zero();
    return vec3_scale(vec3_normalize(r), -G * M / (dist * dist));
}
Vec3 gravity_dipole_term(Vec3 r, Vec3 dipole_moment, double G) {
    double dist = vec3_norm(r);
    if (dist < 1e-15) return vec3_zero();
    Vec3 r_hat = vec3_normalize(r);
    double d_dot_r = vec3_dot(dipole_moment, r_hat);
    Vec3 term1 = vec3_scale(r_hat, 3.0 * d_dot_r);
    return vec3_scale(vec3_sub(term1, dipole_moment), G / (dist * dist * dist));
}
Vec3 gravity_quadrupole_term(Vec3 r, double Q_tensor[3][3], double G) {
    double dist = vec3_norm(r);
    if (dist < 1e-15) return vec3_zero();
    double Q_eff = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            Q_eff += Q_tensor[i][j] * Q_tensor[i][j];
    Q_eff = sqrt(Q_eff);
    double factor = -G * Q_eff / (dist * dist * dist * dist);
    return vec3_scale(vec3_normalize(r), factor);
}
