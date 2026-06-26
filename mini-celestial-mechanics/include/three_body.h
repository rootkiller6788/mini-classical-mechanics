/**
 * three_body.h — Circular Restricted Three-Body Problem (CRTBP)
 *
 * Lagrange points, Jacobi integral, zero-velocity surfaces,
 * Hill/SoI/Roche limits, and CRTBP trajectory integration.
 *
 * The CRTBP models the motion of a massless particle under the
 * gravitational influence of two massive bodies in circular orbits
 * about their common barycenter.
 *
 * References:
 *   - Murray & Dermott, Solar System Dynamics, Sec.3.8-3.13
 *   - Szebehely, Theory of Orbits: The Restricted Problem of
 *     Three Bodies (1967), Ch.1-4
 *   - Vallado, Sec.10.6
 *
 * Knowledge: L6 Canonical Systems (CRTBP, Lagrange points),
 *   L8 Advanced Topics (Hill region, Lyapunov orbits)
 */

#ifndef THREE_BODY_H
#define THREE_BODY_H

#include "celestial_types.h"

/* ---- L6: CRTBP Equations of Motion ---- */

/**
 * CRTBP equations in the rotating (synodic) frame.
 *
 * State: [x, y, z, vx, vy, vz] (position and velocity in rotating frame)
 * mu: mass ratio m2/(m1+m2), m1 at x=-mu, m2 at x=1-mu
 *
 * Equations:
 *   x'' - 2y' = dOmega/dx
 *   y'' + 2x' = dOmega/dy
 *   z'' = dOmega/dz
 *
 * where Omega is the effective potential (pseudo-potential).
 *
 * Returns derivative: [vx, vy, vz, ax, ay, az]
 *
 * Ref: Szebehely Eq.(1.2.14)
 */
void crtbp_equations(const double state[6], double mu, double deriv[6]);

/** CRTBP effective potential (pseudo-potential):
 *  Omega(x,y,z) = (x^2+y^2)/2 + (1-mu)/r1 + mu/r2 + mu*(1-mu)/2
 *  where r1 = sqrt((x+mu)^2 + y^2 + z^2) to m1
 *        r2 = sqrt((x-1+mu)^2 + y^2 + z^2) to m2 */
double crtbp_effective_potential(double x, double y, double z, double mu);

/** Jacobi constant (integral of motion):
 *  C_J = 2*Omega(x,y,z) - v^2
 *  Conserved in the CRTBP (no other forces). */
double jacobi_constant(const double state[6], double mu);

/** Compute z-coordinate of the zero-velocity surface at (x,y)
 *  for a given Jacobi constant C_J.
 *  ZVS: v^2 = 2*Omega - C_J = 0
 *  Returns sqrt(2*Omega - C_J) if >=0, else NaN. */
double zero_velocity_surface_z(double x, double y, double C_J, double mu, double z_guess);

/** Check if a point (x,y,z) is accessible for given Jacobi constant */
int point_accessible(double x, double y, double z, double C_J, double mu);

/* ---- L6: Lagrange Points ---- */

/** Compute all five Lagrange point locations for mass ratio mu */
LagrangePoints lagrange_points(double mu);

/** Compute Jacobi constants at all five Lagrange points */
void lagrange_jacobi_constants(double mu, double C_J[5]);

/** Find a collinear Lagrange point via bisection */
double find_collinear_lagrange(double mu, double x0, double x1, double tol);

/** Classify Hill region connectivity for a given C_J
 *  Returns a string describing which regions are connected. */
const char *hill_region(double C_J, double mu);

/* ---- L8: Stability and Manifolds ---- */

/** Check linear stability of a collinear Lagrange point:
 *  Returns the eigenvalues (complex conjugate pairs).
 *  Positive real part -> unstable. */
void lagrange_point_stability(double mu, int point_idx,
                               double eigenvalues[6]);

/** Compute the state transition matrix (STM) for CRTBP
 *  at a given state. STM is 6x6, stored row-major in phi[36]. */
void crtbp_stm(const double state[6], double mu,
                const double phi[36], double phi_dot[36]);

/* ---- L6: Spheres of Influence ---- */

/** Laplace sphere of influence radius:
 *  r_SOI = a * (m/M)^(2/5) */
double sphere_of_influence(double a, double m, double M);

/** Hill sphere radius:
 *  r_H = a * (m/(3*M))^(1/3) */
double hill_sphere(double a, double m, double M);

/** Roche limit for a rigid satellite:
 *  d = 1.26 * R_primary * (rho_primary/rho_satellite)^(1/3) */
double roche_limit_rigid(double R_primary, double rho_primary, double rho_satellite);

/** Roche limit for a fluid satellite:
 *  d = 2.44 * R_primary * (rho_primary/rho_satellite)^(1/3) */
double roche_limit_fluid(double R_primary, double rho_primary, double rho_satellite);

/* ---- L5: CRTBP Numerical Integration ---- */

/** Integrate CRTBP trajectory using RK4 fixed-step.
 *  state0[6]: initial state in rotating frame
 *  mu: mass ratio
 *  t_end: integration duration (in CRTBP time units)
 *  dt: fixed step size
 *  n_steps_out: output number of steps taken
 *  Returns: array of states [n_steps_out+1][6], caller must free. */
double **integrate_crtbp(const double state0[6], double mu,
                          double t_end, double dt, int *n_steps_out);

/** Free memory allocated by integrate_crtbp */
void free_crtbp_trajectory(double **traj, int n_steps);

/** Single RK4 step for CRTBP */
void crtbp_rk4_step(double state[6], double mu, double dt);

#endif /* THREE_BODY_H */