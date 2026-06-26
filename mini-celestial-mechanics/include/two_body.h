/**
 * two_body.h — Two-Body Problem: vis-viva, Energy, Propagation, Time of Flight
 *
 * The Keplerian two-body problem is the foundation of celestial mechanics.
 * All formulas derive from conservation of energy and angular momentum.
 *
 * References:
 *   - Goldstein, Poole & Safko, Classical Mechanics, Sec.3.5-3.7
 *   - Murray & Dermott, Solar System Dynamics, Sec.2.3-2.5
 *   - Vallado, Fundamentals of Astrodynamics, Sec.2.3
 *
 * Knowledge: L4 Fundamental Laws (vis-viva, conservation),
 *   L2 Core Concepts (orbital energy, period, angular momentum),
 *   L5 Computational (propagation, time-of-flight)
 */

#ifndef TWO_BODY_H
#define TWO_BODY_H

#include "celestial_types.h"

/* ---- L4: Vis-Viva Equation and Energy ---- */

/** Vis-viva: v^2 = mu*(2/r - 1/a) */
double vis_viva(double r, double a, double mu);

/** Specific orbital energy: epsilon = -mu/(2a) */
double orbital_energy(double a, double mu);

/** Orbital period: T = 2*pi*sqrt(a^3/mu) (elliptic only) */
double orbital_period(double a, double mu);

/** Mean motion: n = sqrt(mu/a^3) = 2*pi/T */
double mean_motion(double a, double mu);

/* ---- L2: Orbital Geometry ---- */

/** Radial distance: r(nu) = a*(1-e^2)/(1+e*cos(nu)) */
double radial_distance_elliptic(double a, double e, double nu);

/** Radial distance for hyperbolic orbit: r = a*(e^2-1)/(1+e*cos(nu)) */
double radial_distance_hyperbolic(double a, double e, double nu);

/** Semi-latus rectum: p = a*(1-e^2) */
double semi_latus_rectum(double a, double e);

/** Periapsis distance: r_p = a*(1-e) */
double periapsis_distance(double a, double e);

/** Apoapsis distance: r_a = a*(1+e) (elliptic only, INFINITY if e>=1) */
double apoapsis_distance(double a, double e);

/** Periapsis speed from vis-viva */
double periapsis_velocity(double a, double e, double mu);

/** Apoapsis speed from vis-viva */
double apoapsis_velocity(double a, double e, double mu);

/** Circular orbit speed: v_circ = sqrt(mu/r) */
double circular_velocity(double r, double mu);

/** Escape speed: v_esc = sqrt(2*mu/r) */
double escape_velocity(double r, double mu);

/* ---- L5: Position and Velocity in Orbital Plane ---- */

/** Position in perifocal frame at true anomaly nu */
Vector3 position_in_orbital_plane(double a, double e, double nu);

/** Velocity in perifocal frame at true anomaly nu */
Vector3 velocity_in_orbital_plane(double a, double e, double nu, double mu);

/** Position in perifocal frame for hyperbolic orbit */
Vector3 position_in_orbital_plane_hyperbolic(double a, double e, double nu);

/** Velocity in perifocal frame for hyperbolic orbit */
Vector3 velocity_in_orbital_plane_hyperbolic(double a, double e, double nu, double mu);

/* ---- L5: Time of Flight and Propagation ---- */

/** Time of flight on elliptic arc from nu1 to nu2 */
double time_of_flight_elliptic(double a, double e, double nu1, double nu2, double mu);

/** Time of flight on hyperbolic arc from nu1 to nu2 */
double time_of_flight_hyperbolic(double a, double e, double nu1, double nu2, double mu);

/** Time of flight on parabolic arc (Barker) */
double time_of_flight_parabolic(double p, double nu1, double nu2, double mu);

/** Time since periapsis passage */
double time_since_periapsis(double a, double e, double nu, double mu);

/** Propagate true anomaly by delta_t */
double propagate_anomaly(double a, double e, double nu0, double dt, double mu);

/** Propagate state vector by delta_t using Kepler equation */
void propagate_state_vector(Vector3 r0, Vector3 v0, double dt, double mu,
                             Vector3 *r_out, Vector3 *v_out);

/* ---- L2: Velocity Decomposition ---- */

/** Velocity components: radial (v_r) and transverse (v_theta) */
void velocity_components(double a, double e, double nu, double mu,
                          double *v_r, double *v_theta);

/** Flight path angle gamma = atan(v_r / v_theta) */
double flight_path_angle(double a, double e, double nu);

/* ---- L3: Laplace-Runge-Lenz and Angular Momentum Vectors ---- */

/** Eccentricity vector: e = (v x h)/mu - r_hat */
Vector3 eccentricity_vector(Vector3 r, Vector3 v, double mu);

/** Semi-major axis from state vector: 1/a = 2/r - v^2/mu */
double semi_major_axis_from_state(Vector3 r, Vector3 v, double mu);

/** Eccentricity from state vector (magnitude of e_vec) */
double eccentricity_from_state(Vector3 r, Vector3 v, double mu);

/** Specific angular momentum magnitude: h = |r x v| */
double angular_momentum_magnitude(Vector3 r, Vector3 v);

/** Specific angular momentum vector: h = r x v */
Vector3 angular_momentum_vector(Vector3 r, Vector3 v);

/** Areal velocity (Kepler's Second Law): dA/dt = h/2 */
double areal_velocity(Vector3 r, Vector3 v);

/* ---- L2: Orbit Shape Parameters ---- */

/** Flattening factor: f = (r_a - r_p)/r_a */
double flattening_factor(double a, double e);

/** Semi-minor axis: b = a*sqrt(1-e^2) (elliptic) */
double semi_minor_axis(double a, double e);

/** Impact parameter for hyperbolic orbit: b_inf = |a|*sqrt(e^2-1) */
double impact_parameter_hyperbolic(double a, double e);

/** Bending angle for hyperbolic passage: delta = 2*asin(1/e) */
double bending_angle_hyperbolic(double e);

#endif /* TWO_BODY_H */