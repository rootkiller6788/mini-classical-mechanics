/**
 * constraints.h - Constrained mechanical systems (L2 Concepts / L6 Canonical)
 *
 * Constraints reduce degrees of freedom through geometric relations.
 * Treated via constraint forces (normal forces, tension).
 *
 * Course: MIT 8.012 Ch.7, Kleppner & Kolenkow Ch.3, Goldstein Ch.1.3
 *
 * Key concepts:
 *   - Normal force: N perpendicular to constraint surface
 *   - Tension: constant magnitude along ideal massless string
 *   - Holonomic constraints: f(r1, r2, ..., t) = 0
 *   - D'Alembert's principle: virtual work of constraint forces = 0
 */
#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include "vec3.h"

/* ==== L6: Inclined plane ==== */
double incline_acceleration(double theta, double mu, double g);
double incline_normal_force(double m, double theta, double g);
double angle_of_repose(double mu_s);
double incline_time_to_bottom(double length, double theta, double mu, double g);
double incline_velocity_at_bottom(double length, double theta, double mu, double g);

/* ==== L6: Pulley systems ==== */
void atwood_machine(double m1, double m2, double g,
                    double *a1, double *a2, double *tension);
double pulley_mechanical_advantage(int n_segments, double load_mass, double g);
double pulley_system_acceleration(double m_load, double F_pull,
                                  int n_segments, double g);
/* Generalized Atwood with friction on one side */
void atwood_with_friction(double m1, double m2, double theta,
                          double mu, double g, double *a, double *T);

/* ==== L6: Tension analysis ==== */
Vec3 rope_tension_force(Vec3 attach_point, Vec3 anchor_point, double tension_magnitude);
double tension_in_accelerating_rope(double attached_mass, double acceleration,
                                    double other_force_magnitude, double g);

/* ==== L6: Conical pendulum ==== */
int conical_pendulum_parameters(double m, double L, double omega, double g,
                                double *theta, double *tension);
double conical_pendulum_critical_omega(double L, double g);
double conical_pendulum_radius_at_omega(double L, double omega, double g);

/* ==== L6: Curved surface constraints ==== */
double surface_normal_force(double m, double v, double R, double theta, double g);
int detachment_condition(double m, double v, double R, double theta, double g,
                         double *normal_force);
double loop_the_loop_min_speed(double R, double g);
double loop_the_loop_normal_force_at_top(double m, double v, double R, double g);
double loop_the_loop_normal_force_at_bottom(double m, double v, double R, double g);

/* ==== L6: Incline-pulley combination ==== */
void incline_pulley_system(double m1_on_incline, double m2_hanging,
                           double theta, double mu, double g,
                           double *acceleration, double *tension);
int incline_pulley_direction(double m1, double m2, double theta,
                             double mu_s, double g);

/* ==== L6: Banked curves ==== */
double banked_curve_ideal_angle(double v, double R, double g);
double max_safe_speed_banked(double R, double theta, double mu, double g);
double min_safe_speed_banked(double R, double theta, double mu, double g);

/* ==== L4: Effective gravity in accelerating frames ==== */
Vec3 effective_gravity_vec(Vec3 a_frame, Vec3 g_true);
double apparent_weight_in_elevator(double m, double a_elevator, double g);

/* ==== L6: Simple pendulum (small-angle approximation) ==== */
double pendulum_period_small_angle(double L, double g);
double pendulum_frequency(double L, double g);
double pendulum_period_large_angle(double L, double theta0, double g);
double pendulum_angular_acceleration(double theta, double L, double g);

/* ==== L4: D'Alembert's principle ==== */
/* Virtual work on system with holonomic constraints */
double dalembert_virtual_work(int n_particles, Vec3 applied_forces[],
                              Vec3 virtual_displacements[]);

#endif /* CONSTRAINTS_H */
