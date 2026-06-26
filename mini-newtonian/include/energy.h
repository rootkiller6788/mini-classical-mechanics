/**
 * energy.h - Work, energy, and conservation laws (L2 Concepts / L4 Laws)
 *
 * Core concepts:
 *   - Work-Energy Theorem: W = Delta(K)
 *   - Conservative forces: F = -grad(U), W_path_independent
 *   - Mechanical energy: E = T + U
 *   - Energy conservation in isolated systems
 *
 * Course: MIT 8.012 Ch.5, Goldstein Ch.2, Feynman Vol.I Ch.13
 */
#ifndef ENERGY_H
#define ENERGY_H

#include "vec3.h"

/* ==== L1: Kinetic Energy ==== */
/* T = 1/2 * m * |v|^2  (translational kinetic energy) */
double kinetic_energy(double m, Vec3 v);
double kinetic_energy_scalar(double m, double speed);  /* T = 0.5*m*speed^2 */
double total_kinetic_energy(double masses[], Vec3 velocities[], int n);

/* ==== L1: Potential Energy ==== */
/* Gravitational potential (uniform field): U = m*g*h */
double potential_energy_uniform_gravity(double m, double height, double g);

/* Gravitational potential (Newtonian): U = -G*M*m/r */
double potential_energy_newtonian(Vec3 r_obj, Vec3 r_source, double M, double m, double G);

/* N-body total gravitational potential (each pair counted once) */
double total_potential_newtonian(double masses[], Vec3 positions[], int n, double G);

/* Elastic potential: U = 1/2 * k * |r - r_eq|^2 */
double potential_energy_elastic(Vec3 r, Vec3 r_eq, double k);
double potential_energy_elastic_1d(double x, double x_eq, double k);

/* Effective potential for central force problems: U_eff = U(r) + L^2/(2mr^2) */
double effective_potential(double r, double U_func(double r),
                           double angular_momentum, double m);

/* ==== L2: Work ==== */
/* Work by constant force: W = F . dr */
double work_constant_force(Vec3 F, Vec3 dr);

/* Work along a path (midpoint rule numerical integration) */
double work_along_path(Vec3 (*F_func)(Vec3 position),
                       Vec3 path[], int n_points);

/* Power: P = F . v (instantaneous) */
double power_instantaneous(Vec3 F, Vec3 v);

/* ==== L4: Energy Conservation ==== */
/* Compute total mechanical energy E = T + U for a trajectory */
void compute_energy_series(double m, const Vec3 positions[], const Vec3 velocities[],
                           int n, double (*U_func)(Vec3),
                           double energies[], double *E0);

/* Relative energy drift: (E(t) - E0) / |E0| */
void energy_drift_series(const double energies[], int n, double E0,
                         double drift[]);

/* Check if a force field is conservative (finite difference curl test) */
int is_force_conservative(Vec3 (*F_func)(Vec3 position),
                          Vec3 eval_point, double h);

/* ==== L4: Mechanical Energy Balance ==== */
/* W_ext = Delta(E)  (Work done by external forces = change in mechanical energy) */
double external_work_balance(double E_initial, double E_final,
                             double W_nonconservative);

/* ==== L4: Virial Theorem (Clausius, 1870) ==== */
/* For bound systems: <T>_time = -1/2 * <F.r>_time */
/* For gravity/power-law potentials: <T> = -(n/2)*<U> where U ~ r^n */
double virial_ratio_gravity(double kinetic, double potential);
double virial_theorem_exponent(double kinetic, double potential);

#endif /* ENERGY_H */
