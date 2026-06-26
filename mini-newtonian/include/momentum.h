/**
 * momentum.h - Linear/angular momentum, impulse, collisions (L2/L4)
 *
 * Core theorems:
 *   - Conservation of linear momentum: Sigma(p_i) = const (isolated system)
 *   - Conservation of angular momentum: Sigma(L_i) = const (central forces)
 *   - Impulse-momentum: J = Delta(p)
 *   - Center of mass theorem: F_ext = M*a_cm
 *
 * Course: MIT 8.012 Ch.6, Goldstein Ch.2-3
 */
#ifndef MOMENTUM_H
#define MOMENTUM_H

#include "vec3.h"

/* ==== L1: Single-particle momentum ==== */
Vec3 linear_momentum(double m, Vec3 v);
Vec3 angular_momentum(Vec3 r, Vec3 v, double m);       /* L = r x p */
Vec3 angular_momentum_about(Vec3 r, Vec3 v, double m, Vec3 origin);

/* ==== L2: System momentum ==== */
Vec3 total_linear_momentum(double masses[], Vec3 velocities[], int n);
Vec3 total_angular_momentum(double masses[], Vec3 positions[],
                            Vec3 velocities[], int n);
Vec3 center_of_mass(double masses[], Vec3 positions[], int n);
Vec3 center_of_mass_velocity(double masses[], Vec3 velocities[], int n);
double reduced_mass(double m1, double m2);

/* ==== L4: Impulse ==== */
Vec3 impulse_constant_force(Vec3 F, double dt);
Vec3 impulse_from_momentum_change(double m, Vec3 v_initial, Vec3 v_final);
Vec3 average_force_from_impulse(Vec3 impulse, double dt);
int verify_impulse_momentum(double m, Vec3 v_before, Vec3 v_after,
                            Vec3 impulse, double tol);

/* ==== L6: 1D Collisions ==== */
/* Elastic: both momentum and kinetic energy conserved */
void elastic_collision_1d(double m1, double v1, double m2, double v2,
                          double *v1f, double *v2f);

/* Perfectly inelastic: stick together, momentum conserved */
double inelastic_collision_1d(double m1, double v1, double m2, double v2);

/* Coefficient of restitution: e = |v2f - v1f| / |v2i - v1i| */
double coefficient_of_restitution(double v1i, double v2i,
                                  double v1f, double v2f);

/* Collision with arbitrary restitution coefficient e */
void collision_with_restitution(double m1, double v1, double m2, double v2,
                                double e, double *v1f, double *v2f);

/* Energy loss in partially inelastic collision */
double collision_energy_loss(double m1, double v1i, double m2, double v2i,
                             double v1f, double v2f);

/* ==== L6: 3D Collisions ==== */
/* Elastic collision of smooth spheres in 3D */
void elastic_collision_3d(double m1, Vec3 v1, double m2, Vec3 v2,
                          Vec3 collision_normal,
                          Vec3 *v1f, Vec3 *v2f);

/* ==== L7: Rocket equation (Tsiolkovsky, 1903) ==== */
double tsiolkovsky_delta_v(double ve, double m0, double mf);
double mass_ratio_for_delta_v(double delta_v, double ve);
double rocket_mass_flow_rate(double thrust, double ve);
double rocket_final_mass(double m0, double burn_time, double mass_rate);

/* ==== L4: Conservation law verification ==== */
int check_momentum_conservation(Vec3 p_initial, Vec3 p_final, double tol);
int check_angular_momentum_conservation(Vec3 L_initial, Vec3 L_final, double tol);
int check_energy_conservation(double E_initial, double E_final, double tol);

#endif /* MOMENTUM_H */
