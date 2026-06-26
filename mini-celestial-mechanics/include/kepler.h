/**
 * kepler.h — Kepler Equation Solvers and Orbital Element Conversions
 *
 * The Kepler equation is the transcendental equation that relates
 * time to position in the two-body problem:
 *
 *   Elliptic:   M = E - e sin(E)          (Kepler 1609)
 *   Hyperbolic: M = e sinh(H) - H          (Lambert 1761)
 *   Parabolic:  M_p = tan(ν/2) + (1/3)tan³(ν/2)  (Barker 1757)
 *   Universal:  √μ t = a·x + ...           (Battin, universal variable)
 *
 * References:
 *   - Goldstein, Poole & Safko §3.5-3.7
 *   - Murray & Dermott §2.4-2.8
 *   - Vallado §2.2-2.5
 *   - Danby, Fundamentals of Celestial Mechanics (1988), Ch.6
 *
 * Knowledge: L1 Definitions, L4 Fundamental Laws, L5 Newton iteration
 */

#ifndef KEPLER_H
#define KEPLER_H

#include "celestial_types.h"

/* ================================================================
 * L4: Kepler Equation Solvers
 * ================================================================ */

/**
 * Solve elliptic Kepler equation: M = E - e sin(E)
 *
 * Uses the Newton-Raphson method with Danby's (1988) cubic-convergent
 * initial guess for robust convergence. Guaranteed convergence for
 * 0 ≤ e < 1 with ≤12 iterations to machine precision.
 *
 * @param M  mean anomaly [rad], may be un-wrapped (reduced mod 2π internally)
 * @param e  eccentricity, 0 ≤ e < 1
 * @param tol convergence tolerance (default: 1e-12)
 * @return   eccentric anomaly E [rad] ∈ [0, 2π)
 *
 * Ref: Danby (1988), Murray & Dermott Eq. (2.68)
 */
double solve_kepler(double M, double e, double tol);

/**
 * Solve hyperbolic Kepler equation: M_h = e sinh(H) - H
 *
 * Newton-Raphson with asymptotic initial guess.
 *
 * @param M_h  hyperbolic mean anomaly
 * @param e    eccentricity (> 1.0)
 * @param tol  convergence tolerance
 * @return     hyperbolic eccentric anomaly H
 *
 * Ref: Vallado §2.2.4
 */
double solve_kepler_hyperbolic(double M_h, double e, double tol);

/**
 * Solve Barker's equation for parabolic orbits (e = 1)
 *
 * Uses the analytic cubic formula:
 *   Let A = (3/2)M_p, B = (A + √(A²+1))^(1/3)
 *   Then tan(ν/2) = B - 1/B
 *
 * @param M_p  parabolic mean anomaly
 * @return     true anomaly ν [rad]
 *
 * Ref: Barker (1757), Murray & Dermott §2.6
 */
double solve_kepler_parabolic(double M_p);

/**
 * Universal Kepler equation solver — automatically selects the correct
 * formulation based on eccentricity.
 *
 * @param x  anomaly (M for elliptic, M_h for hyperbolic, M_p for parabolic)
 * @param e  eccentricity
 * @return   corresponding anomaly (E for e<1, H for e>1, ν for e=1)
 */
double solve_kepler_universal(double x, double e, double tol);

/**
 * Universal variable formulation of Kepler's equation (Battin's method).
 *
 * Given position r0, radial velocity σ0, time interval Δt, and α=1/a,
 * solve for the universal variable x using Newton's method.
 *
 * Used for robust propagation of any conic orbit (elliptic, parabolic,
 * hyperbolic) with a single unified formulation.
 *
 * @param r0    initial radius magnitude
 * @param sigma0 initial radial velocity (r·v)/√μ
 * @param dt    time interval √μ * Δt
 * @param alpha 1/a (reciprocal semi-major axis)
 * @param tol   convergence tolerance
 * @return      universal variable x (such that Δt = ...)
 *
 * Ref: Battin §4.4, Vallado §2.4
 */
double solve_kepler_universal_variable(double r0, double sigma0,
                                        double dt, double alpha, double tol);

/* ---- Stumpff functions (for universal variable formulation) ---- */

/** Stumpff c-function: c(z) = sum (-z)^k/(2k+2)!
 *  c(0)=1/2, c(z>0)=(1-cos(sqrt(z)))/z, c(z<0)=(cosh(sqrt(-z))-1)/(-z) */
double stumpff_c(double z);

/** Stumpff s-function: s(z) = sum (-z)^k/(2k+3)!
 *  s(0)=1/6, s(z>0)=(sqrt(z)-sin(sqrt(z)))/(z*sqrt(z)),
 *  s(z<0)=(sinh(sqrt(-z))-sqrt(-z))/((-z)*sqrt(-z)) */
double stumpff_s(double z);

/* ================================================================
 * L2: Anomaly Conversions (elliptic)
 * ================================================================ */

/**
 * Eccentric anomaly → True anomaly
 *   tan(ν/2) = √((1+e)/(1-e)) · tan(E/2)
 */
double eccentric_to_true_anomaly(double E, double e);

/**
 * True anomaly → Eccentric anomaly
 *   tan(E/2) = √((1-e)/(1+e)) · tan(ν/2)
 */
double true_to_eccentric_anomaly(double nu, double e);

/**
 * Eccentric anomaly → Mean anomaly
 *   M = E - e sin(E)   (Kepler's equation forward)
 */
double eccentric_to_mean_anomaly(double E, double e);

/**
 * Mean anomaly → Eccentric anomaly
 *   Solves M = E - e sin(E) via solve_kepler()
 */
double mean_to_eccentric_anomaly(double M, double e, double tol);

/**
 * True anomaly → Mean anomaly
 *   ν → E → M
 */
double true_to_mean_anomaly(double nu, double e);

/**
 * Mean anomaly → True anomaly
 *   M → E → ν
 */
double mean_to_true_anomaly(double M, double e, double tol);

/* ---- Hyperbolic anomaly conversions ---- */

/**
 * Hyperbolic eccentric anomaly → True anomaly
 *   tan(ν/2) = √((e+1)/(e-1)) · tanh(H/2)
 */
double hyperbolic_to_true_anomaly(double H, double e);

/**
 * True anomaly → Hyperbolic eccentric anomaly
 *   tanh(H/2) = √((e-1)/(e+1)) · tan(ν/2)
 */
double true_to_hyperbolic_anomaly(double nu, double e);

/**
 * Hyperbolic mean anomaly → Hyperbolic eccentric anomaly
 *   M_h = e sinh(H) - H
 */
double mean_to_hyperbolic_anomaly(double M_h, double e, double tol);

/* ================================================================
 * L4: Orbital Elements ↔ State Vector (Cartesian)
 * ================================================================ */

/**
 * Convert classical orbital elements to Cartesian state vector (r, v).
 *
 * Algorithm (3-1-3 Euler rotation):
 *   1. Compute position and velocity in the orbit plane
 *        r_orb = r_mag [cos ν, sin ν, 0]^T
 *        v_orb = μ/h [-sin ν, e + cos ν, 0]^T
 *   2. Rotate by R_z(Omega) · R_x(i) · R_z(omega)
 *        to transform from perifocal to ECI/ECL frame
 *
 * @param el       orbital elements
 * @param mu       gravitational parameter
 * @param r_out    output position vector (3D)
 * @param v_out    output velocity vector (3D)
 *
 * Ref: Murray & Dermott Eq. (2.122), Vallado §2.5.1
 */
void orbital_elements_to_state(const OrbitalElements *el, double mu,
                                Vector3 *r_out, Vector3 *v_out);

/**
 * Convert Cartesian state vector to classical orbital elements.
 *
 * Algorithm:
 *   1. h = r × v (specific angular momentum)
 *   2. Determine inclination i from h_z/h
 *   3. Node vector n = k × h, determine Ω from n
 *   4. Semi-major axis from vis-viva: 1/a = 2/r - v²/μ
 *   5. Eccentricity vector: e_vec = (v×h)/μ - r̂
 *   6. Argument of periapsis ω from n·e
 *   7. True anomaly ν from e·r
 *
 * All quadrant ambiguities are resolved using the sign of the
 * z-component of the relevant vectors.
 *
 * @param r    position vector
 * @param v    velocity vector
 * @param mu   gravitational parameter
 * @param el   output orbital elements
 *
 * Ref: Murray & Dermott §2.8, Vallado §2.5.2
 */
void state_to_orbital_elements(Vector3 r, Vector3 v, double mu,
                                OrbitalElements *el);

/**
 * Round-trip validation: elements → state → elements
 *
 * Returns the maximum absolute error among a, e, i, Omega, omega.
 * Should be < 1e-12 for non-degenerate orbits.
 *
 * @param el  input elements
 * @param mu  gravitational parameter
 * @return    maximum error among the five shape/orientation elements
 */
double kepler_roundtrip_error(const OrbitalElements *el, double mu);

/**
 * Compute the position and velocity magnitudes at a given true anomaly
 * for an elliptic orbit.
 *
 * r(ν) = a(1-e²)/(1+e cos ν)
 * v(ν) = √(μ(2/r - 1/a))  (vis-viva)
 *
 * @param a     semi-major axis
 * @param e     eccentricity
 * @param nu    true anomaly [rad]
 * @param mu    gravitational parameter
 * @param r_mag output radial distance
 * @param v_mag output speed
 */
void orbit_state_at_anomaly(double a, double e, double nu, double mu,
                             double *r_mag, double *v_mag);

#endif /* KEPLER_H */
