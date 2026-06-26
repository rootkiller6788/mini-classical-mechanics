/**
 * forces.h — Newtonian force laws (L1 Definitions / L4 Fundamental Laws)
 *
 * All forces in Newtonian mechanics. Each force law is an independent
 * knowledge point with a specific physical origin.
 *
 * Course mapping:
 *   MIT 8.012 §3.1-4.4 Newton's Laws and Forces
 *   Goldstein §1.3-1.5 Mechanics of a particle
 *   Feynman Lectures Vol.I §12 Force
 *
 * Core theorem: Newton's Second Law F = ma
 *   (implemented via acceleration functions in integrators)
 */

#ifndef FORCES_H
#define FORCES_H

#include "vec3.h"

/* ===== L4: Universal Gravitation (Newton, 1687) ===== */
/* F_12 = -G * M1 * M2 / r^2 * r_hat
 * Knowledge points:
 *  - Inverse-square law
 *  - Superposition principle for N bodies
 *  - Uniform gravity as near-surface approximation
 */
Vec3 newton_gravity_force(Vec3 r_obj, Vec3 r_source, double M, double m, double G);
Vec3 gravity_acceleration_field(Vec3 r_obj, Vec3 r_source, double M, double G);
Vec3 gravity_nbody_force(Vec3 r_obj, double m, Vec3 sources_pos[],
                         double sources_mass[], int n_sources, double G);
Vec3 uniform_gravity_force(double m, Vec3 direction, double g);
Vec3 gravity_potential_gradient(Vec3 r_obj, Vec3 r_source, double M, double G);

/* ===== L4: Hooke's Law (Hooke, 1660) ===== */
/* F = -k * (x - x_eq)
 * Knowledge points: linear elasticity, equilibrium position
 */
Vec3 hooke_force_3d(Vec3 r, Vec3 r_eq, double k);
double hooke_force_1d(double x, double x_eq, double k);
double hooke_anharmonic_force_1d(double x, double x_eq, double k, double beta);

/* ===== L4: Drag forces (Stokes 1851, Rayleigh) ===== */
/* Knowledge points:
 *  - Linear drag (Stokes): F = -b*v (low Reynolds number, laminar)
 *  - Quadratic drag: F = -c*|v|*v (high Reynolds number, turbulent)
 *  - Combined model with Reynolds-dependent transition
 */
Vec3 linear_drag(Vec3 v, double b);
Vec3 quadratic_drag(Vec3 v, double c);
Vec3 combined_drag(Vec3 v, double b, double c);
double drag_reynolds_threshold(double rho_fluid, double eta,
                               double characteristic_length);

/* ===== L4: Friction (Amontons/Coulomb, 1699/1785) ===== */
/* Knowledge points:
 *  - Kinetic (sliding) friction: f_k = mu_k * N
 *  - Static friction: f_s <= mu_s * N
 *  - Angle of repose: theta_c = arctan(mu_s)
 */
Vec3 kinetic_friction_force(Vec3 velocity, double normal_mag, double mu_k);
double static_friction_max(double normal_mag, double mu_s);
double static_friction_threshold(double normal_mag, double mu_s,
                                 double applied_tangential);

/* ===== L4: Lorentz force (Lorentz, 1895) ===== */
/* F = q*(E + v x B)
 * Bridge to electromagnetism module (mini-electromagnetism)
 */
Vec3 lorentz_force(double q, Vec3 E, Vec3 v, Vec3 B);
double cyclotron_frequency(double q, double B_mag, double m);
double cyclotron_radius(double v_perp, double q, double B_mag, double m);
Vec3 magnetic_mirror_force(double mu_magnetic_moment, Vec3 B_gradient);

/* ===== L2: Composite force functions ===== */
/* Knowledge points: superposition of independent force laws
 * Each function demonstrates composition of building-block forces
 */
Vec3 damped_harmonic_force(Vec3 r, Vec3 r_eq, Vec3 v, double k, double b);
Vec3 driven_damped_force(Vec3 r, Vec3 r_eq, Vec3 v, double k, double b,
                         double F0, double omega_drive, double t);
Vec3 van_der_pol_force(Vec3 r, Vec3 r_eq, Vec3 v, double k, double mu);
Vec3 duffing_force_3d(Vec3 r, Vec3 r_eq, double k, double alpha, double beta);

/* ===== L4/L7: Gravitational multipole expansion ===== */
/* Knowledge point: far-field expansion for distributed mass
 * monopole (M/r) + dipole (CM offset) + quadrupole correction
 */
Vec3 gravity_monopole_term(Vec3 r, double M, double G);
Vec3 gravity_dipole_term(Vec3 r, Vec3 dipole_moment, double G);
Vec3 gravity_quadrupole_term(Vec3 r, double Q_tensor[3][3], double G);

#endif /* FORCES_H */
