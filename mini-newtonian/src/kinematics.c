/**
 * kinematics.c - Classical kinematics implementation
 *
 * Each function implements one independent kinematics knowledge point.
 * All functions O(1), no numerical integration.
 *
 * Theorems:
 *   Galileo (1638): constant-acceleration laws, projectile parabola theorem
 *   Newton (1687): centripetal acceleration v^2/r for circular motion
 *   Coriolis (1835): rotating-frame apparent forces
 *   Frenet (1847): space curve Frenet frame
 *
 * Course mapping: MIT 8.012 Ch.1-2, Goldstein Ch.1
 */
#include "kinematics.h"
#include <math.h>

/* ===== L1: 1D constant-acceleration motion (SUVAT equations) ===== */

double suvat_position(double r0, double v0, double a, double t) {
    return r0 + v0 * t + 0.5 * a * t * t;
}
double suvat_velocity(double v0, double a, double t) {
    return v0 + a * t;
}
double suvat_velocity_from_displacement(double v0, double a, double dr) {
    double v_sq = v0 * v0 + 2.0 * a * dr;
    if (v_sq < 0.0) return 0.0;
    return sqrt(v_sq);
}

/* ===== L1: 3D constant-acceleration motion ===== */

Vec3 uniform_motion(Vec3 r0, Vec3 v, double t) {
    return vec3_add(r0, vec3_scale(v, t));
}
Vec3 accelerated_motion(Vec3 r0, Vec3 v0, Vec3 a, double t) {
    Vec3 vt = vec3_scale(v0, t);
    Vec3 at2 = vec3_scale(a, 0.5 * t * t);
    return vec3_add(r0, vec3_add(vt, at2));
}
Vec3 accelerated_velocity(Vec3 v0, Vec3 a, double t) {
    return vec3_add(v0, vec3_scale(a, t));
}

/* ===== L2: Projectile motion - Galileo projectile theorem (1638) ===== */
/* Knowledge points:
 * 1. Trajectory parabola: superposition of horizontal uniform + vertical accelerated
 * 2. Flight time: solve z(t)=0 quadratic
 * 3. Range formula: R = v0^2 sin(2*theta)/g
 * 4. Maximum height: h_max = v0z^2/(2g)
 * 5. Range from elevated launch
 * 6. Optimal launch angle (inverse problem)
 */

ProjectileState projectile_at_time(Vec3 r0, Vec3 v0, double t, Vec3 gravity) {
    ProjectileState ps;
    Vec3 gt = vec3_scale(gravity, 0.5 * t * t);
    ps.position = vec3_add(r0, vec3_add(vec3_scale(v0, t), gt));
    ps.velocity = vec3_add(v0, vec3_scale(gravity, t));
    return ps;
}
double projectile_flight_time(double v0z, double z0, double g) {
    double disc = v0z * v0z + 2.0 * g * z0;
    if (disc < 0.0) return -1.0;
    return (v0z + sqrt(disc)) / g;
}
double projectile_range(double v0, double theta, double g) {
    return v0 * cos(theta) * projectile_flight_time(v0 * sin(theta), 0.0, g);
}
double projectile_max_height(double v0, double theta, double g) {
    double v0z = v0 * sin(theta);
    return (v0z * v0z) / (2.0 * g);
}
double projectile_range_from_height(double v0, double theta, double h0, double g) {
    double v0x = v0 * cos(theta);
    double v0z = v0 * sin(theta);
    double disc = v0z * v0z + 2.0 * g * h0;
    if (disc < 0.0) return 0.0;
    return v0x * (v0z + sqrt(disc)) / g;
}
double optimal_launch_angle(double v0, double target_range, double g) {
    double sin2theta = target_range * g / (v0 * v0);
    if (fabs(sin2theta) > 1.0) return -1.0;
    return 0.5 * asin(sin2theta);
}
double projectile_trajectory_height(double x, double v0, double theta, double g) {
    double cos_t = cos(theta);
    return x * tan(theta) - (g * x * x) / (2.0 * v0 * v0 * cos_t * cos_t);
}

/* ===== L3: Uniform circular motion ===== */
/* Knowledge points:
 * 1. Parametric position on circle
 * 2. Tangential velocity v = omega * r
 * 3. Centripetal acceleration a_c = v^2/r = omega^2 * r
 * 4. Period T = 2*pi/omega
 */

Vec3 circular_position(Vec3 center, Vec3 axis_normal,
                       double radius, double omega, double t) {
    Vec3 axis = vec3_normalize(axis_normal);
    Vec3 e1, e2;
    if (fabs(axis.z) < 0.999) {
        e1 = vec3_normalize(vec3_cross(axis, vec3_make(0, 0, 1)));
    } else {
        e1 = vec3_normalize(vec3_cross(axis, vec3_make(1, 0, 0)));
    }
    e2 = vec3_cross(axis, e1);
    double angle = omega * t;
    Vec3 offset = vec3_add(
        vec3_scale(e1, radius * cos(angle)),
        vec3_scale(e2, radius * sin(angle))
    );
    return vec3_add(center, offset);
}
Vec3 circular_velocity(Vec3 axis_normal, double radius, double omega, double t) {
    Vec3 axis = vec3_normalize(axis_normal);
    Vec3 e1, e2;
    if (fabs(axis.z) < 0.999) {
        e1 = vec3_normalize(vec3_cross(axis, vec3_make(0, 0, 1)));
    } else {
        e1 = vec3_normalize(vec3_cross(axis, vec3_make(1, 0, 0)));
    }
    e2 = vec3_cross(axis, e1);
    double angle = omega * t;
    return vec3_add(
        vec3_scale(e1, -radius * omega * sin(angle)),
        vec3_scale(e2,  radius * omega * cos(angle))
    );
}
Vec3 circular_acceleration(Vec3 center, Vec3 position, double omega) {
    Vec3 radial = vec3_sub(position, center);
    return vec3_scale(radial, -omega * omega);
}
double centripetal_acceleration(double v, double r) {
    if (fabs(r) < 1e-15) return 0.0;
    return (v * v) / r;
}
double centripetal_acceleration_omega(double omega, double r) {
    return omega * omega * r;
}
double circular_period(double omega) {
    if (fabs(omega) < 1e-15) return 1e300;
    return 2.0 * M_PI / omega;
}
double circular_frequency(double period) {
    if (fabs(period) < 1e-15) return 1e300;
    return 2.0 * M_PI / period;
}

/* ===== L3: Non-uniform circular motion - acceleration decomposition ===== */
/* Knowledge point: a_total = a_centripetal + a_tangential
 *   a_c points toward center, a_t is along velocity
 */
CurvilinearAccel decompose_acceleration(Vec3 position, Vec3 velocity,
                                        Vec3 acceleration, Vec3 center) {
    CurvilinearAccel ca;
    Vec3 radial_dir = vec3_normalize(vec3_sub(position, center));
    Vec3 vel_dir = vec3_normalize(velocity);
    ca.centripetal_dir = vec3_neg(radial_dir);
    ca.tangential_dir = vel_dir;
    ca.centripetal_mag = vec3_dot(acceleration, ca.centripetal_dir);
    ca.tangential_mag  = vec3_dot(acceleration, ca.tangential_dir);
    return ca;
}

/* ===== L2: Galilean transformations ===== */
/* Knowledge point: Galileo relativity - physics laws identical in all inertial frames
 *   r = r_prime + V*t,  v = v_prime + V,  a = a_prime (invariant!)
 */
Vec3 galilean_position_transform(Vec3 r_prime, Vec3 V_frame, double t) {
    return vec3_add(r_prime, vec3_scale(V_frame, t));
}
Vec3 galilean_velocity_transform(Vec3 v_prime, Vec3 V_frame) {
    return vec3_add(v_prime, V_frame);
}
Vec3 galilean_acceleration_invariance(Vec3 a_prime) {
    return a_prime;
}

/* ===== L4: Rotating reference frames ===== */
/* Knowledge points - Non-inertial apparent forces:
 * 1. Coriolis: -2*omega x v_rot  (deflects moving objects)
 * 2. Centrifugal: -omega x (omega x r)  (outward radial)
 * 3. Euler: -d(omega)/dt x r  (angular acceleration term)
 * Theorem: a_inertial = a_rot + 2*omega x v_rot + omega x (omega x r) + d(omega)/dt x r
 * Course: MIT 8.012 Sec.7.2, Goldstein Sec.4.10
 */
Vec3 coriolis_acceleration(Vec3 v_rot, Vec3 omega) {
    return vec3_scale(vec3_cross(omega, v_rot), -2.0);
}
Vec3 centrifugal_acceleration(Vec3 r_rot, Vec3 omega) {
    Vec3 cross_prod = vec3_cross(omega, r_rot);
    return vec3_neg(vec3_cross(omega, cross_prod));
}
Vec3 euler_acceleration(Vec3 r, Vec3 omega_dot) {
    return vec3_neg(vec3_cross(omega_dot, r));
}
Vec3 rotating_frame_accel(Vec3 a_inertial, Vec3 v_rot,
                          Vec3 r_rot, Vec3 omega) {
    Vec3 cor = coriolis_acceleration(v_rot, omega);
    Vec3 cen = centrifugal_acceleration(r_rot, omega);
    return vec3_add(a_inertial, vec3_add(cor, cen));
}
Vec3 inertial_to_rotating_velocity(Vec3 v_inertial, Vec3 r, Vec3 omega) {
    return vec3_sub(v_inertial, vec3_cross(omega, r));
}

/* ===== L5: Frenet-Serret frame for space curves ===== */
/* Knowledge points (differential geometry of curves):
 *   T = v/|v|  (unit tangent)
 *   N = dT/dt / |dT/dt|  (principal normal)
 *   B = T x N  (binormal)
 *   curvature kappa = |v x a| / |v|^3
 *   radius of curvature rho = 1/kappa
 * Course: Calculus on manifolds, MIT 8.012 supplemental
 */
double trajectory_curvature(Vec3 velocity, Vec3 acceleration) {
    double v = vec3_norm(velocity);
    if (v < 1e-15) return 0.0;
    return vec3_norm(vec3_cross(velocity, acceleration)) / (v * v * v);
}
double trajectory_radius_of_curvature(Vec3 velocity, Vec3 acceleration) {
    double kappa = trajectory_curvature(velocity, acceleration);
    if (kappa < 1e-15) return 1e300;
    return 1.0 / kappa;
}
Vec3 frenet_tangent(Vec3 velocity) {
    return vec3_normalize(velocity);
}
Vec3 frenet_normal(Vec3 velocity, Vec3 acceleration) {
    Vec3 T = frenet_tangent(velocity);
    double v = vec3_norm(velocity);
    if (v < 1e-15) return vec3_zero();
    Vec3 a_proj_T = vec3_scale(T, vec3_dot(acceleration, T));
    Vec3 a_normal = vec3_sub(acceleration, a_proj_T);
    return vec3_normalize(a_normal);
}
Vec3 frenet_binormal(Vec3 velocity, Vec3 acceleration) {
    return vec3_cross(frenet_tangent(velocity),
                      frenet_normal(velocity, acceleration));
}
