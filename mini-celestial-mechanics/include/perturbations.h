/**
 * perturbations.h — Orbital Perturbations
 *
 * Models for non-Keplerian forces that disturb the ideal two-body orbit:
 * J2 (oblateness), atmospheric drag, solar radiation pressure,
 * third-body gravitation, and general relativistic precession.
 *
 * References:
 *   - Vallado, Fundamentals of Astrodynamics, Sec.8.2-8.7, Sec.9.2
 *   - Murray & Dermott, Solar System Dynamics, Sec.10.2-10.4
 *   - Kaula, Theory of Satellite Geodesy (1966)
 *   - Goldstein Sec.7.6 (GR precession)
 *
 * Knowledge:
 *   L4 Fundamental Laws — Gauss/Lagrange planetary equations
 *   L5 Computational — combined perturbation propagation
 *   L6 Canonical Systems — J2 problem, sun-synchronous orbits
 *   L7 Applications — satellite drag, SRP modeling
 */

#ifndef PERTURBATIONS_H
#define PERTURBATIONS_H

#include "celestial_types.h"

/* ---- L6: J2 Gravitational Perturbation ---- */

/** J2 acceleration in ECI frame [km/s^2] */
Vector3 j2_acceleration(Vector3 r, double mu, double J2, double R_body);

/** J2 secular (long-period) element rates:
 *  Omega_dot (RAAN drift), omega_dot (periapsis drift), M_dot */
J2SecularRates j2_secular_rates(double a, double e, double i,
                                 double mu, double J2, double R_body);

/** Sun-synchronous orbit inclination for given (a, e):
 *  solves cos(i) = -Omega_sun_dot / (1.5*n*J2*(R/p)^2) */
double sun_sync_inclination(double a, double e, double mu, double J2, double R_body);

/** J2+J4 combined acceleration (higher-order geopotential) */
Vector3 j2_j4_acceleration(Vector3 r, double mu, double J2, double J4, double R_body);

/* ---- L7: Atmospheric Drag ---- */

/** Atmospheric drag acceleration (exponential density model)
 *  a_drag = -0.5*C_D*(A/m)*rho*v^2 * v_hat
 *  where rho = rho0*exp(-(r-r0)/H) */
Vector3 atmospheric_drag_acceleration(Vector3 r, Vector3 v,
                                       double C_D, double A_over_m,
                                       double rho0, double r0, double H);

/** Drag-induced orbital decay rate (circular orbit approximation)
 *  da/dt = -C_D*(A/m)*rho*n*a^2 */
double drag_decay_rate(double a, double C_D, double A_over_m, double rho, double mu);

/** Estimate remaining orbital lifetime from drag (simplified)
 *  Uses energy dissipation rate and assumes exponential atmosphere. */
double drag_lifetime_estimate(double a0, double e0, double C_D,
                               double A_over_m, double rho0, double r0,
                               double H, double mu);

/** Harris-Priester atmospheric density model (altitude-dependent) */
double atmospheric_density_harris_priester(double altitude_km, double F10_7);

/* ---- L7: Solar Radiation Pressure (SRP) ---- */

/** SRP acceleration with cylindrical shadow model
 *  a_srp = -C_R*(A/m)*P_sun*(1 AU)^2 * (r_sat-r_sun)/|r_sat-r_sun|^3 */
Vector3 srp_acceleration(Vector3 r_sat, Vector3 r_sun,
                          double C_R, double A_over_m, double P_sun);

/** Shadow function: 0=umbra, 0-1=penumbra, 1=sunlight */
double shadow_function(Vector3 r_sat, Vector3 r_sun, double R_body, double R_sun);

/* ---- L7: Third-Body Perturbation ---- */

/** Third-body gravitational acceleration (direct + indirect)
 *  a_3b = -mu_b*[ (r_sat-r_b)/|r_sat-r_b|^3 + r_b/|r_b|^3 ] */
Vector3 third_body_acceleration(Vector3 r_sat, Vector3 r_body, double mu_body);

/** Solar perturbation on Earth satellite */
Vector3 solar_perturbation_acceleration(Vector3 r_sat, Vector3 r_sun);

/** Lunar perturbation on Earth satellite */
Vector3 lunar_perturbation_acceleration(Vector3 r_sat, Vector3 r_moon);

/* ---- L8: General Relativistic Corrections ---- */

/** GR periapsis precession rate [rad/s]
 *  omega_dot_GR = 3*n*mu/(a*c^2*(1-e^2))
 *  For Mercury: ~43 arcsec/century */
double gr_precession_rate(double a, double e, double mu, double c);

/** GR correction to Newtonian acceleration (post-Newtonian 1PN)
 *  a_GR = -(mu/r^2)*[ (1+3*v^2/c^2)*r_hat - ... ] */
Vector3 gr_post_newtonian_acceleration(Vector3 r, Vector3 v, double mu, double c);

/** Classical (Newtonian) precession rate from perturbing mass */
double precession_rate_classical(double a, double m_perturber, double M_central);

/* ---- L5: Combined Perturbation Models ---- */

/** Propagate orbital elements with J2 secular rates only */
OrbitalElements propagate_with_j2(const OrbitalElements *el, double dt,
                                   double mu, double J2, double R_body);

/** Sum of all selected perturbation accelerations at a point */
Vector3 combined_perturbation(Vector3 r, Vector3 v, Vector3 r_sun,
                               Vector3 r_moon, double mu, double J2,
                               double R_body, double C_D, double A_over_m,
                               double rho0, double r0, double H,
                               double C_R, double P_sun, PerturbationFlags flags);

/** Gauss planetary equations: da/dt, de/dt, di/dt, ... from RSW perturbation */
void gauss_planetary_equations(const OrbitalElements *el, double mu,
                                const PerturbationRSW *accel,
                                double *a_dot, double *e_dot, double *i_dot,
                                double *Omega_dot, double *omega_dot, double *M_dot);

#endif /* PERTURBATIONS_H */