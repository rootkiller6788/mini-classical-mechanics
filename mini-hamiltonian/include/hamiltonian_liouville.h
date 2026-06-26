/**
 * @file hamiltonian_liouville.h
 * @brief Liouville theorem, phase space volume conservation, and ergodic theory.
 *
 * Reference: Goldstein Ch.9, Arnold Ch.5, Khinchin Mathematical Foundations
 *            of Statistical Mechanics
 *
 * L2 Core Concepts: Phase space density and Liouville equation
 * L4 Fundamental Laws: Liouville theorem — volume preservation under Hamiltonian flow
 * L8 Advanced Topics: Ergodicity, mixing, Poincare recurrence, KS-entropy
 *
 * University mapping:
 *   MIT 8.012 — Liouville's theorem and phase space
 *   Stanford PHYSICS 370 — Statistical mechanics foundations
 *   Oxford CMT — Ergodic theory and Hamiltonian chaos
 */

#ifndef HAMILTONIAN_LIOUVILLE_H
#define HAMILTONIAN_LIOUVILLE_H

#include "hamiltonian_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L4: Liouville theorem — core verification
 * ================================================================ */

/**
 * Generate an ensemble of phase points around (q0, p0).
 *
 * Points are drawn from a Gaussian distribution with std sigma_q, sigma_p
 * in each coordinate/momentum direction, using Box-Muller transform.
 *
 * This creates a "phase space cloud" suitable for testing Liouville's theorem:
 * evolve the cloud and track the volume proxy.
 *
 * @param n_particles  Number of particles in the ensemble
 * @param q0, p0       Center of the ensemble
 * @param sigma_q, sigma_p  Spread in coordinates and momenta
 * @param ensemble_q   Output: ensemble coordinates, n_particles × n, row-major
 * @param ensemble_p   Output: ensemble momenta, n_particles × n, row-major
 */
void generate_ensemble(size_t n, size_t n_particles,
                        const double *q0, const double *p0,
                        double sigma_q, double sigma_p,
                        double *ensemble_q, double *ensemble_p);

/**
 * Compute a proxy for phase space volume of an ensemble.
 *
 * Uses the product of RMS spreads: V_proxy = sigma_q * sigma_p
 * where sigma_q = RMS of |q - <q>| and sigma_p = RMS of |p - <p>|.
 *
 * This is a coarse proxy; for exact volume one would need the convex hull
 * or determinant of the covariance matrix.
 */
double phase_volume_proxy(size_t n, size_t n_particles,
                           const double *ensemble_q,
                           const double *ensemble_p);

/**
 * Verify Liouville's theorem: evolve an ensemble and check volume conservation.
 *
 * Theorem (Liouville, 1838):
 *   The phase flow of Hamilton's equations preserves the volume element
 *   dV = dq_1...dq_n dp_1...dp_n in phase space.
 *
 * Equivalently: the velocity field in phase space has zero divergence.
 *   div(v) = sum_i (d/dq_i(dq_i/dt) + d/dp_i(dp_i/dt))
 *          = sum_i (d^2 H/dq_i dp_i - d^2 H/dp_i dq_i) = 0
 *
 * The proof relies on equality of mixed partial derivatives (Schwarz theorem).
 *
 * @param tol  Relative volume change below which Liouville is considered verified
 * @return     Relative volume drift = (V_final - V_initial) / |V_initial|
 */
double verify_liouville(const HamiltonianSystem *sys,
                         const double *q0, const double *p0,
                         size_t n_particles, double sigma_q, double sigma_p,
                         double t_end, double dt,
                         IntegratorType integrator);

/**
 * Compute the phase space velocity divergence at a point.
 *
 * div(v) = sum_i [ d/dq_i(dq_i/dt) + d/dp_i(dp_i/dt) ]
 *        = sum_i [ d/dq_i(dH/dp_i) + d/dp_i(-dH/dq_i) ]
 *        = sum_i [ d^2H/dq_i dp_i - d^2H/dp_i dq_i ]
 *
 * For smooth H, this is identically zero by equality of mixed partials.
 * Non-zero values indicate either:
 *   (a) non-Hamiltonian flow (dissipative/external forces), or
 *   (b) numerical error in finite difference approximation.
 */
double phase_velocity_divergence(const HamiltonianSystem *sys,
                                  const double *q, const double *p,
                                  double eps);

/* ================================================================
 * Liouville equation — density evolution
 * ================================================================ */

/**
 * Evolve a phase space density rho(q, p) according to the Liouville equation.
 *
 *   drho/dt + {rho, H} = 0
 *
 * Using upwind finite differences on a uniform grid.
 * The Hamiltonian flow advects the density in phase space without
 * changing its value along trajectories.
 *
 * @param rho       Initial density grid (n_grid × n_grid), modified in-place
 * @param q_range   Grid range in q
 * @param p_range   Grid range in p
 * @param n_grid    Grid resolution per dimension
 * @param t_end, dt Integration parameters
 */
void evolve_liouville_density(LiouvilleDensity *rho,
                               const HamiltonianSystem *sys,
                               double t_end, double dt);

/**
 * Compute the Gibbs entropy of a phase space density.
 *
 *   S = -k_B ∫ rho(q,p) ln(rho(q,p)) dq dp
 *
 * For a Hamiltonian flow: dS/dt = 0 (entropy is conserved).
 * In contrast, non-Hamiltonian (dissipative) systems have dS/dt > 0.
 *
 * @param kB  Boltzmann constant (set to 1.0 for dimensionless entropy)
 */
double gibbs_entropy(const LiouvilleDensity *rho, double kB);

/**
 * Compute the fine-grained entropy vs. coarse-grained entropy.
 *
 * Fine-grained entropy is conserved under Hamiltonian flow.
 * Coarse-grained entropy (computed on a coarse grid) increases
 * due to filamentation of the density in phase space — this is
 * the microscopic origin of the second law in classical mechanics.
 *
 * @param coarse_factor  Coarse-graining factor (e.g., 2 = 2x2 -> 1x1 blocks)
 */
double coarse_grained_entropy(const LiouvilleDensity *rho, int coarse_factor, double kB);

/* ================================================================
 * L7: Applications — ergodicity and mixing
 * ================================================================ */

/**
 * Estimate the mixing rate from an ensemble spread.
 *
 * Mixing: initial phase space region stretches and folds under the flow,
 * becoming uniformly distributed. Measure via the growth of inter-particle
 * distances in the ensemble.
 *
 * @param dist_mean  Output: mean pairwise phase space distance
 * @param dist_std   Output: standard deviation of pairwise distances
 */
void mixing_proxy(size_t n, size_t n_particles,
                   const double *ensemble_q, const double *ensemble_p,
                   double *dist_mean, double *dist_std);

/**
 * Test ergodicity: long-time average = ensemble average.
 *
 * Birkhoff's ergodic theorem (1931):
 *   For an ergodic flow, the time average of any integrable observable
 *   equals its phase space average for almost all initial conditions.
 *
 *   lim_{T→∞} (1/T)∫_0^T f(phi_t(x)) dt = ∫_P f dmu
 *
 * @param observable  Function to average
 * @param t_end       Integration time for time averaging
 * @return            Time-averaged value
 */
double ergodicity_time_average(const HamiltonianSystem *sys,
                                const double *q0, const double *p0,
                                Observable observable,
                                double t_end, double dt);

/**
 * Poincare recurrence check.
 *
 * Theorem (Poincare, 1890):
 *   For a volume-preserving flow on a bounded phase space domain,
 *   almost every trajectory returns arbitrarily close to its initial
 *   condition infinitely often.
 *
 * This function returns the first recurrence time (< t_max) or INFINITY.
 *
 * @param eps_recurrence  Distance threshold for "returned"
 * @return  First recurrence time, or INFINITY if none found within t_max
 */
double poincare_recurrence_time(const HamiltonianSystem *sys,
                                 const double *q0, const double *p0,
                                 double t_max, double dt,
                                 double eps_recurrence);

/**
 * Estimate Kolmogorov-Sinai (KS) entropy from a phase space ensemble.
 *
 * KS entropy measures the rate at which information about initial
 * conditions is produced by the dynamics.
 *
 *   KS = sum_{lambda_i > 0} lambda_i   (Pesin's theorem)
 *
 * where lambda_i are the Lyapunov exponents.
 * KS > 0 indicates chaos; KS = 0 for integrable/regular motion.
 *
 * This function provides a finite-time estimate using the divergence
 * of nearby trajectories.
 */
double ks_entropy_estimate(const HamiltonianSystem *sys,
                            size_t n, size_t n_particles,
                            const double *ensemble_q,
                            const double *ensemble_p,
                            double dt);

#ifdef __cplusplus
}
#endif

#endif /* HAMILTONIAN_LIOUVILLE_H */
