/**
 * energy.c - Work, energy, and conservation law implementations
 * Each function implements an independent knowledge point.
 * Course: MIT 8.012 Ch.5, Goldstein Ch.2, Feynman Vol.I Ch.13-14
 */
#include "energy.h"
#include <math.h>

/* L1: Kinetic Energy T = 1/2 * m * |v|^2 */
double kinetic_energy(double m, Vec3 v) {
    return 0.5 * m * vec3_norm2(v);
}
double kinetic_energy_scalar(double m, double speed) {
    return 0.5 * m * speed * speed;
}
double total_kinetic_energy(double masses[], Vec3 velocities[], int n) {
    double T = 0.0;
    for (int i = 0; i < n; i++) T += kinetic_energy(masses[i], velocities[i]);
    return T;
}

/* L1: Potential Energy */
double potential_energy_uniform_gravity(double m, double height, double g) {
    return m * g * height;
}
double potential_energy_newtonian(Vec3 r_obj, Vec3 r_source, double M, double m, double G) {
    double dist = vec3_norm(vec3_sub(r_obj, r_source));
    if (dist < 1e-10) return -1e300;
    return -G * M * m / dist;
}
double total_potential_newtonian(double masses[], Vec3 positions[], int n, double G) {
    double U = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            U += potential_energy_newtonian(positions[i], positions[j], masses[j], masses[i], G);
    return U;
}
double potential_energy_elastic(Vec3 r, Vec3 r_eq, double k) {
    return 0.5 * k * vec3_norm2(vec3_sub(r, r_eq));
}
double potential_energy_elastic_1d(double x, double x_eq, double k) {
    double dx = x - x_eq;
    return 0.5 * k * dx * dx;
}

/* L2: Effective potential for central force problems
 * U_eff(r) = U(r) + L^2/(2*m*r^2)
 * The centrifugal barrier L^2/(2mr^2) acts as repulsive potential.
 * Theorem: Angular momentum conservation reduces 3D central force
 * problem to 1D radial motion in U_eff.
 * Course: MIT 8.012 Ch.9, Goldstein Ch.3.3
 */
double effective_potential(double r, double U_func(double r),
                           double angular_momentum, double m) {
    if (r < 1e-15 || m < 1e-15) return 1e300;
    return U_func(r) + (angular_momentum * angular_momentum) / (2.0 * m * r * r);
}

/* L2: Work
 * W = integral_C F.dr
 * For constant force: W = F.dr
 * For variable force: numerical midpoint rule
 * Power: P = dW/dt = F.v
 */
double work_constant_force(Vec3 F, Vec3 dr) {
    return vec3_dot(F, dr);
}
double work_along_path(Vec3 (*F_func)(Vec3 position), Vec3 path[], int n_points) {
    double W = 0.0;
    for (int i = 1; i < n_points; i++) {
        Vec3 dr = vec3_sub(path[i], path[i-1]);
        Vec3 midpoint = vec3_scale(vec3_add(path[i], path[i-1]), 0.5);
        W += vec3_dot(F_func(midpoint), dr);
    }
    return W;
}
double power_instantaneous(Vec3 F, Vec3 v) {
    return vec3_dot(F, v);
}

/* L4: Energy Conservation
 * For conservative system: dE/dt = 0 where E = T + U.
 * Energy drift = (E(t) - E0) / |E0| measures integrator quality.
 * Symplectic integrators: bounded oscillation. Non-symplectic: monotonic drift.
 */
void compute_energy_series(double m, const Vec3 positions[], const Vec3 velocities[],
                           int n, double (*U_func)(Vec3),
                           double energies[], double *E0) {
    for (int i = 0; i < n; i++) {
        double T = kinetic_energy(m, velocities[i]);
        energies[i] = T + U_func(positions[i]);
    }
    *E0 = energies[0];
}
void energy_drift_series(const double energies[], int n, double E0, double drift[]) {
    double absE0 = fabs(E0);
    if (absE0 < 1e-300) absE0 = 1.0;
    for (int i = 0; i < n; i++) drift[i] = (energies[i] - E0) / absE0;
}

/* L4: Conservative force check via finite-difference curl
 * F(r) is conservative iff curl(F) = 0 everywhere.
 * Numerical test: evaluate curl at point using central differences.
 */
int is_force_conservative(Vec3 (*F_func)(Vec3 position), Vec3 eval_point, double h) {
    Vec3 p = eval_point;
    Vec3 F_xp = F_func(vec3_make(p.x + h, p.y, p.z));
    Vec3 F_xm = F_func(vec3_make(p.x - h, p.y, p.z));
    Vec3 F_yp = F_func(vec3_make(p.x, p.y + h, p.z));
    Vec3 F_ym = F_func(vec3_make(p.x, p.y - h, p.z));
    Vec3 F_zp = F_func(vec3_make(p.x, p.y, p.z + h));
    Vec3 F_zm = F_func(vec3_make(p.x, p.y, p.z - h));
    double curl_x = ((F_zp.y - F_zm.y) - (F_yp.z - F_ym.z)) / (2.0 * h);
    double curl_y = ((F_xp.z - F_xm.z) - (F_zp.x - F_zm.x)) / (2.0 * h);
    double curl_z = ((F_yp.x - F_ym.x) - (F_xp.y - F_xm.y)) / (2.0 * h);
    double curl_norm = sqrt(curl_x*curl_x + curl_y*curl_y + curl_z*curl_z);
    return (curl_norm < 1e-6) ? 1 : 0;
}

/* L4: Mechanical Energy Balance
 * W_ext + W_nc = Delta(E) where W_nc is work by non-conservative forces.
 */
double external_work_balance(double E_initial, double E_final, double W_nc) {
    return (E_final - E_initial) - W_nc;
}

/* L4: Virial Theorem (Clausius 1870)
 * For bound systems: <T> = -1/2 <sum F_i.r_i>
 * For power-law potentials U ~ r^n: 2<T> = n<U>
 * Gravity (n=-1): <T> = -<U>/2 ; Harmonic (n=2): <T> = <U> (equipartition)
 * Course: Goldstein Ch.3.4, Pathria Statistical Mechanics Ch.3
 */
double virial_ratio_gravity(double kinetic, double potential) {
    if (fabs(potential) < 1e-15) return 0.0;
    return (kinetic / (-0.5 * potential)) - 1.0;
}
double virial_theorem_exponent(double kinetic, double potential) {
    if (fabs(potential) < 1e-15) return 0.0;
    return 2.0 * kinetic / potential;
}
