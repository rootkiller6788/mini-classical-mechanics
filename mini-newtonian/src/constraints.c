/**
 * constraints.c - Constrained mechanical systems
 * Each function implements an independent knowledge point.
 * Course: MIT 8.012 Ch.7, Kleppner & Kolenkow Ch.3, Goldstein Ch.1.3
 */
#include "constraints.h"
#include <math.h>
#include <stdlib.h>

/* ================================================================
 * L6: Inclined plane
 *
 * Knowledge point: Decomposing gravity into components parallel
 * and perpendicular to the incline. Friction opposes motion.
 *
 * a_parallel = g*(sin(theta) - mu*cos(theta))
 * N = mg*cos(theta)
 * Angle of repose: theta_c = arctan(mu_s)
 * ================================================================ */

double incline_acceleration(double theta, double mu, double g) {
    return g * (sin(theta) - mu * cos(theta));
}
double incline_normal_force(double m, double theta, double g) {
    return m * g * cos(theta);
}
double angle_of_repose(double mu_s) {
    return atan(mu_s);
}
double incline_time_to_bottom(double length, double theta, double mu, double g) {
    double a = incline_acceleration(theta, mu, g);
    if (a <= 0.0) return -1.0;  /* no motion or slides up */
    return sqrt(2.0 * length / a);
}
double incline_velocity_at_bottom(double length, double theta, double mu, double g) {
    double a = incline_acceleration(theta, mu, g);
    if (a <= 0.0) return 0.0;
    return sqrt(2.0 * a * length);
}

/* ================================================================
 * L6: Pulley systems
 *
 * Knowledge point: Ideal pulley redirects tension force.
 * Massless, frictionless string: tension constant throughout.
 *
 * Atwood machine:
 *   a = (m1 - m2)*g / (m1 + m2)  [m1 down positive]
 *   T = 2*m1*m2*g / (m1 + m2)
 *
 * Mechanical advantage: n rope segments -> force multiplied by n
 * ================================================================ */

void atwood_machine(double m1, double m2, double g,
                    double *a1, double *a2, double *tension) {
    double M = m1 + m2;
    double a = (m1 - m2) * g / M;  /* positive = m1 down */
    *a1 = a;
    *a2 = -a;
    *tension = 2.0 * m1 * m2 * g / M;
}
double pulley_mechanical_advantage(int n_segments, double load_mass, double g) {
    return load_mass * g / n_segments;  /* effort needed */
}
double pulley_system_acceleration(double m_load, double F_pull,
                                  int n_segments, double g) {
    double F_net = n_segments * F_pull - m_load * g;
    return F_net / m_load;
}
void atwood_with_friction(double m1, double m2, double theta,
                          double mu, double g, double *a, double *T) {
    /* m1 on incline (angle theta, friction mu), m2 hanging vertically */
    double M = m1 + m2;
    double f1_parallel = m1 * g * sin(theta);
    double f1_friction = mu * m1 * g * cos(theta);
    double a_val = (m2 * g - f1_parallel - f1_friction) / M;
    *a = a_val;
    *T = m2 * (g - a_val);
}

/* ================================================================
 * L6: Tension analysis
 * ================================================================ */

Vec3 rope_tension_force(Vec3 attach_point, Vec3 anchor_point, double tension_magnitude) {
    Vec3 direction = vec3_normalize(vec3_sub(anchor_point, attach_point));
    return vec3_scale(direction, tension_magnitude);
}
double tension_in_accelerating_rope(double attached_mass, double acceleration,
                                    double other_force_magnitude, double g) {
    return attached_mass * (g + acceleration) + other_force_magnitude;
}

/* ================================================================
 * L6: Conical pendulum
 *
 * Knowledge point: Circular motion with a string constraint.
 * T*cos(theta) = mg  (vertical balance)
 * T*sin(theta) = m*omega^2*L*sin(theta) = m*v^2/(L*sin(theta))
 * -> cos(theta) = g/(omega^2*L)
 *
 * Critical omega: omega_min = sqrt(g/L)
 * ================================================================ */

int conical_pendulum_parameters(double m, double L, double omega, double g,
                                double *theta, double *tension) {
    double omega2L = omega * omega * L;
    if (omega2L < g) return 0;  /* omega too small, not a cone */
    double cos_theta = g / omega2L;
    if (cos_theta > 1.0 || cos_theta <= 0.0) return 0;
    *theta = acos(cos_theta);
    *tension = m * g / cos_theta;
    return 1;
}
double conical_pendulum_critical_omega(double L, double g) {
    return sqrt(g / L);
}
double conical_pendulum_radius_at_omega(double L, double omega, double g) {
    double omega2L = omega * omega * L;
    if (omega2L < g) return 0.0;
    double cos_theta = g / omega2L;
    return L * sqrt(1.0 - cos_theta * cos_theta);
}

/* ================================================================
 * L6: Curved surface constraints
 *
 * Knowledge point: Normal force on a curved surface has two contributions:
 *   1. mg*cos(theta) from gravity component perpendicular to surface
 *   2. m*v^2/R from centripetal acceleration requirement
 *
 * Detachment condition: N <= 0 (surface can only push, not pull)
 *
 * Loop-the-loop:
 *   At top: N + mg = mv^2/R -> N = m(v^2/R - g)
 *   Minimum speed at top: v_min = sqrt(gR)
 *   At bottom: N - mg = mv^2/R -> N = m(v^2/R + g)
 * ================================================================ */

double surface_normal_force(double m, double v, double R, double theta, double g) {
    return m * (g * cos(theta) + v * v / R);
}
int detachment_condition(double m, double v, double R, double theta, double g,
                         double *normal_force) {
    *normal_force = surface_normal_force(m, v, R, theta, g);
    return (*normal_force <= 0.0) ? 1 : 0;
}
double loop_the_loop_min_speed(double R, double g) {
    return sqrt(g * R);  /* at the top */
}
double loop_the_loop_normal_force_at_top(double m, double v, double R, double g) {
    return m * (v * v / R - g);  /* can be negative = fall off */
}
double loop_the_loop_normal_force_at_bottom(double m, double v, double R, double g) {
    return m * (v * v / R + g);
}

/* ================================================================
 * L6: Incline-pulley combination
 *
 * Knowledge point: Combined constraints require solving coupled
 * equations. m1 on incline connected via pulley to hanging m2.
 *
 * Equation: m2*g - T = m2*a  (vertical)
 *           T - m1*g*sin(theta) - mu*m1*g*cos(theta) = m1*a  (incline)
 * -> a = (m2 - m1*(sin(theta) + mu*cos(theta))) * g / (m1 + m2)
 * ================================================================ */

void incline_pulley_system(double m1_on_incline, double m2_hanging,
                           double theta, double mu, double g,
                           double *acceleration, double *tension) {
    double M = m1_on_incline + m2_hanging;
    double f_incline = m1_on_incline * (sin(theta) + mu * cos(theta));
    double a = (m2_hanging - f_incline) * g / M;
    *acceleration = a;
    *tension = m2_hanging * (g - a);
}
int incline_pulley_direction(double m1, double m2, double theta,
                             double mu_s, double g) {
    double force_m2 = m2 * g;
    double force_m1_down = m1 * g * sin(theta);
    double friction_limit = mu_s * m1 * g * cos(theta);
    if (fabs(force_m2 - force_m1_down) <= friction_limit) return 0;  /* static */
    return (force_m2 > force_m1_down) ? 1 : -1;  /* 1=m1 up, -1=m1 down */
}

/* ================================================================
 * L6: Banked curves
 *
 * Knowledge point: On a frictionless banked curve, the normal
 * force provides the centripetal component.
 * tan(theta) = v^2/(gR) -> theta_ideal = arctan(v^2/(gR))
 *
 * With friction (mu):
 * v_max = sqrt(gR * (tan(theta) + mu) / (1 - mu*tan(theta)))
 * v_min = sqrt(gR * (tan(theta) - mu) / (1 + mu*tan(theta)))
 * ================================================================ */

double banked_curve_ideal_angle(double v, double R, double g) {
    double tan_theta = v * v / (g * R);
    return atan(tan_theta);
}
double max_safe_speed_banked(double R, double theta, double mu, double g) {
    double tan_t = tan(theta);
    double denom = 1.0 - mu * tan_t;
    if (denom <= 0.0) return 1e300;  /* no slip possible */
    return sqrt(g * R * (tan_t + mu) / denom);
}
double min_safe_speed_banked(double R, double theta, double mu, double g) {
    double tan_t = tan(theta);
    double denom = 1.0 + mu * tan_t;
    double num = tan_t - mu;
    if (num < 0.0) return 0.0;  /* friction holds even at rest */
    return sqrt(g * R * num / denom);
}

/* ================================================================
 * L4: Effective gravity in accelerating frames
 * g_eff = g_true - a_frame
 * In an elevator: apparent weight = m*(g + a)  (a>0 upward)
 * ================================================================ */

Vec3 effective_gravity_vec(Vec3 a_frame, Vec3 g_true) {
    return vec3_sub(g_true, a_frame);
}
double apparent_weight_in_elevator(double m, double a_elevator, double g) {
    return m * (g + a_elevator);
}

/* ================================================================
 * L6: Simple pendulum
 *
 * Knowledge point:
 * Small-angle: T = 2*pi*sqrt(L/g), omega = sqrt(g/L), theta(t) = theta0*cos(omega*t)
 * Large-angle: T = 2*pi*sqrt(L/g) * (1 + theta0^2/16 + 11*theta0^4/3072 + ...)
 * Full equation: d^2(theta)/dt^2 = -(g/L)*sin(theta)
 *
 * Course: MIT 8.012 Ch.10, Goldstein Ch.6
 * ================================================================ */

double pendulum_period_small_angle(double L, double g) {
    return 2.0 * M_PI * sqrt(L / g);
}
double pendulum_frequency(double L, double g) {
    return sqrt(g / L);
}
double pendulum_period_large_angle(double L, double theta0, double g) {
    /* Series expansion to 4th order:
     * T = T0 * (1 + theta0^2/16 + 11*theta0^4/3072 + 173*theta0^6/737280)
     */
    double T0 = pendulum_period_small_angle(L, g);
    double t2 = theta0 * theta0;
    double t4 = t2 * t2;
    double t6 = t4 * t2;
    double factor = 1.0 + t2/16.0 + 11.0*t4/3072.0 + 173.0*t6/737280.0;
    return T0 * factor;
}
double pendulum_angular_acceleration(double theta, double L, double g) {
    return -(g / L) * sin(theta);
}

/* ================================================================
 * L4: D'Alembert's principle
 *
 * Knowledge point: The virtual work of constraint forces is zero.
 * For a system in static/dynamic equilibrium with holonomic constraints:
 * sum_i (F_applied_i + F_constraint_i) . delta_r_i = 0
 * where delta_r_i are virtual displacements compatible with constraints.
 *
 * This is the foundation for Lagrangian mechanics.
 * ================================================================ */

double dalembert_virtual_work(int n_particles, Vec3 applied_forces[],
                              Vec3 virtual_displacements[]) {
    double delta_W = 0.0;
    for (int i = 0; i < n_particles; i++) {
        delta_W += vec3_dot(applied_forces[i], virtual_displacements[i]);
    }
    return delta_W;
}
