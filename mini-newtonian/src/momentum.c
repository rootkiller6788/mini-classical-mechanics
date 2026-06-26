/**
 * momentum.c - Linear/angular momentum, impulse, collisions
 * Each function implements an independent physics knowledge point.
 * Course: MIT 8.012 Ch.6, Goldstein Ch.2-3, Kleppner & Kolenkow Ch.4
 */
#include "momentum.h"
#include <math.h>

/* === L1: Single-particle momentum ===
 * p = m*v
 * L = r x p = m(r x v)
 */
Vec3 linear_momentum(double m, Vec3 v) {
    return vec3_scale(v, m);
}
Vec3 angular_momentum(Vec3 r, Vec3 v, double m) {
    return vec3_scale(vec3_cross(r, v), m);
}
Vec3 angular_momentum_about(Vec3 r, Vec3 v, double m, Vec3 origin) {
    return angular_momentum(vec3_sub(r, origin), v, m);
}

/* === L2: System momentum === */
Vec3 total_linear_momentum(double masses[], Vec3 velocities[], int n) {
    Vec3 p_total = vec3_zero();
    for (int i = 0; i < n; i++)
        p_total = vec3_add(p_total, linear_momentum(masses[i], velocities[i]));
    return p_total;
}
Vec3 total_angular_momentum(double masses[], Vec3 positions[], Vec3 velocities[], int n) {
    Vec3 L_total = vec3_zero();
    for (int i = 0; i < n; i++)
        L_total = vec3_add(L_total, angular_momentum(positions[i], velocities[i], masses[i]));
    return L_total;
}

/* Center of mass: R_cm = sum(m_i*r_i) / sum(m_i)
 * Knowledge point: CM moves as if all mass and external forces concentrated there.
 * Theorem: F_ext = M * a_cm (Center of Mass Theorem)
 */
Vec3 center_of_mass(double masses[], Vec3 positions[], int n) {
    Vec3 cm = vec3_zero();
    double M_total = 0.0;
    for (int i = 0; i < n; i++) {
        cm = vec3_add(cm, vec3_scale(positions[i], masses[i]));
        M_total += masses[i];
    }
    if (M_total < 1e-15) return vec3_zero();
    return vec3_div(cm, M_total);
}
Vec3 center_of_mass_velocity(double masses[], Vec3 velocities[], int n) {
    Vec3 v_cm = vec3_zero();
    double M_total = 0.0;
    for (int i = 0; i < n; i++) {
        v_cm = vec3_add(v_cm, vec3_scale(velocities[i], masses[i]));
        M_total += masses[i];
    }
    if (M_total < 1e-15) return vec3_zero();
    return vec3_div(v_cm, M_total);
}

/* Reduced mass: mu = m1*m2 / (m1+m2)
 * Knowledge point: Two-body problem reduces to one-body with reduced mass.
 * Course: Goldstein Ch.3.1
 */
double reduced_mass(double m1, double m2) {
    return m1 * m2 / (m1 + m2);
}

/* === L4: Impulse J = F*dt = Delta(p) ===
 * Impulse-Momentum Theorem: integral(F dt) = p_final - p_initial
 */
Vec3 impulse_constant_force(Vec3 F, double dt) {
    return vec3_scale(F, dt);
}
Vec3 impulse_from_momentum_change(double m, Vec3 v_initial, Vec3 v_final) {
    return vec3_scale(vec3_sub(v_final, v_initial), m);
}
Vec3 average_force_from_impulse(Vec3 impulse, double dt) {
    if (fabs(dt) < 1e-15) return vec3_zero();
    return vec3_div(impulse, dt);
}
int verify_impulse_momentum(double m, Vec3 v_before, Vec3 v_after,
                            Vec3 impulse, double tol) {
    Vec3 delta_p = vec3_scale(vec3_sub(v_after, v_before), m);
    double diff = vec3_norm(vec3_sub(delta_p, impulse));
    return (diff < tol) ? 1 : 0;
}

/* === L6: 1D Collisions ===
 * Elastic: both momentum and kinetic energy conserved
 *   v1f = ((m1-m2)/(m1+m2))*v1 + (2*m2/(m1+m2))*v2
 *   v2f = (2*m1/(m1+m2))*v1 + ((m2-m1)/(m1+m2))*v2
 */
void elastic_collision_1d(double m1, double v1, double m2, double v2,
                          double *v1f, double *v2f) {
    double M = m1 + m2;
    *v1f = (m1 - m2) / M * v1 + 2.0 * m2 / M * v2;
    *v2f = 2.0 * m1 / M * v1 + (m2 - m1) / M * v2;
}

/* Perfectly inelastic: stick together, v_common = (m1*v1 + m2*v2)/(m1+m2)
 * Kinetic energy is NOT conserved (maximum loss), only momentum conserved.
 */
double inelastic_collision_1d(double m1, double v1, double m2, double v2) {
    return (m1 * v1 + m2 * v2) / (m1 + m2);
}

/* Coefficient of restitution: e = |v2f - v1f| / |v2i - v1i|
 * e = 1: perfectly elastic
 * e = 0: perfectly inelastic
 * 0 < e < 1: partially elastic (real collisions)
 */
double coefficient_of_restitution(double v1i, double v2i, double v1f, double v2f) {
    double rel_i = fabs(v2i - v1i);
    double rel_f = fabs(v2f - v1f);
    if (rel_i < 1e-15) return 1.0;
    return rel_f / rel_i;
}

/* Collision with arbitrary restitution e (0 <= e <= 1)
 * Uses center-of-mass frame: v_cm' = -e * v_cm
 */
void collision_with_restitution(double m1, double v1, double m2, double v2,
                                double e, double *v1f, double *v2f) {
    double v_cm = (m1 * v1 + m2 * v2) / (m1 + m2);
    double v1_cm = v1 - v_cm;
    double v2_cm = v2 - v_cm;
    *v1f = -e * v1_cm + v_cm;
    *v2f = -e * v2_cm + v_cm;
}

/* Energy loss in partially inelastic collision
 * Delta(E) = E_final - E_initial (negative for energy loss)
 */
double collision_energy_loss(double m1, double v1i, double m2, double v2i,
                             double v1f, double v2f) {
    double E_i = 0.5 * m1 * v1i * v1i + 0.5 * m2 * v2i * v2i;
    double E_f = 0.5 * m1 * v1f * v1f + 0.5 * m2 * v2f * v2f;
    return E_f - E_i;
}

/* === L6: 3D Collisions ===
 * Smooth sphere elastic collision in 3D.
 * Normal components exchange elastically, tangential components unchanged.
 * normal: unit vector from sphere1 center to sphere2 center at contact.
 */
void elastic_collision_3d(double m1, Vec3 v1, double m2, Vec3 v2,
                          Vec3 collision_normal,
                          Vec3 *v1f, Vec3 *v2f) {
    Vec3 n = vec3_normalize(collision_normal);
    double v1n = vec3_dot(v1, n);
    double v2n = vec3_dot(v2, n);
    Vec3 v1t = vec3_sub(v1, vec3_scale(n, v1n));
    Vec3 v2t = vec3_sub(v2, vec3_scale(n, v2n));
    double v1nf, v2nf;
    elastic_collision_1d(m1, v1n, m2, v2n, &v1nf, &v2nf);
    *v1f = vec3_add(v1t, vec3_scale(n, v1nf));
    *v2f = vec3_add(v2t, vec3_scale(n, v2nf));
}

/* === L7: Rocket Equation (Tsiolkovsky 1903) ===
 * Delta_v = v_e * ln(m0/mf)
 * Knowledge points:
 *   1. Variable mass system: thrust from expelled propellant
 *   2. Logarithmic dependence: diminishing returns from more fuel
 *   3. Mass ratio = exp(Delta_v / v_e): exponential mass penalty
 * Course: MIT 8.012 Ch.6.4, Kleppner & Kolenkow Ch.4.6
 */
double tsiolkovsky_delta_v(double ve, double m0, double mf) {
    if (mf <= 0.0 || m0 <= mf) return 0.0;
    return ve * log(m0 / mf);
}
double mass_ratio_for_delta_v(double delta_v, double ve) {
    return exp(delta_v / ve);
}
double rocket_mass_flow_rate(double thrust, double ve) {
    if (fabs(ve) < 1e-15) return 0.0;
    return thrust / ve;
}
double rocket_final_mass(double m0, double burn_time, double mass_rate) {
    double mf = m0 - mass_rate * burn_time;
    return (mf > 0.0) ? mf : 0.0;
}

/* === L4: Conservation law verification === */
int check_momentum_conservation(Vec3 p_initial, Vec3 p_final, double tol) {
    return vec3_norm(vec3_sub(p_final, p_initial)) < tol ? 1 : 0;
}
int check_angular_momentum_conservation(Vec3 L_initial, Vec3 L_final, double tol) {
    return vec3_norm(vec3_sub(L_final, L_initial)) < tol ? 1 : 0;
}
int check_energy_conservation(double E_initial, double E_final, double tol) {
    return fabs(E_final - E_initial) < tol * fmax(fabs(E_initial), 1.0) ? 1 : 0;
}
