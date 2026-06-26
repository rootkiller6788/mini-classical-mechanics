/**
 * @file hamiltonian_hamilton_jacobi.h
 * @brief Hamilton-Jacobi theory and integrable systems.
 *
 * Reference: Goldstein Ch.10, Arnold Ch.9, Landau-Lifshitz Vol.1 Ch.7
 *            Evans, Partial Differential Equations — viscosity solutions
 *
 * L2 Core Concepts: Hamilton-Jacobi equation, Hamilton's principal function
 * L3 Mathematical Structures: First-order nonlinear PDE, characteristics method
 * L4 Fundamental Laws: Equivalence of H-J, Hamilton, and Lagrangian formulations
 * L6 Canonical Systems: Separation of variables, action-angle from H-J
 *
 * University mapping:
 *   MIT 8.012 — Hamilton-Jacobi theory
 *   Stanford PHYSICS 230 — H-J and semiclassical limit
 *   Cambridge Part II — Integrable PDEs and separation of variables
 *   Caltech Ph 106 — Hamilton-Jacobi and old quantum theory
 */

#ifndef HAMILTONIAN_HAMILTON_JACOBI_H
#define HAMILTONIAN_HAMILTON_JACOBI_H

#include "hamiltonian_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * Hamilton-Jacobi equation
 * ================================================================ */

/**
 * Solve the time-independent Hamilton-Jacobi equation for 1D systems.
 *
 * H(q, dW/dq) = E  →  p = dW/dq
 *
 * W(q, E) = integral_{q0}^{q} p(q', E) dq'
 *
 * where p(q, E) is obtained by solving H(q, p) = E for p (taking the
 * positive branch for the forward-moving part of the orbit).
 *
 * W(q, E) is Hamilton's characteristic function. It generates the
 * canonical transformation to action-angle variables.
 *
 * @param H_func     H(q_vec, p_vec) -> scalar
 * @param q_range    [q_min, q_max] for evaluation
 * @param E          Energy level
 * @param n_points   Number of grid points
 * @param q_vals     Output: q grid values [n_points]
 * @param W_vals     Output: W(q, E) values [n_points]
 */
void hamilton_jacobi_1d(double (*H_func)(const double *q, const double *p, size_t n),
                         double q_min, double q_max,
                         double E, size_t n_points,
                         double *q_vals, double *W_vals);

/**
 * Solve p from H(q, p) = E using bisection for 1D systems.
 *
 * Assumes H(q, p) is monotonic in |p| (typical for kinetic energy T(p)).
 * Searches in [0, p_max] for the magnitude; sign determined by orbit phase.
 *
 * @param q      Coordinate value [1]
 * @param E      Target energy
 * @param p_max  Maximum search range for |p|
 * @param tol    Convergence tolerance
 * @return       p solution (positive branch)
 */
double solve_p_from_H_1d(double (*H_func)(const double *q, const double *p, size_t n),
                          double q, double E, double p_max, double tol);

/**
 * Compute Hamilton's characteristic function for the harmonic oscillator.
 *
 * H = p^2/(2m) + (1/2) m omega^2 q^2 = E
 * p = sqrt(2mE - m^2 omega^2 q^2)
 *
 * Analytic integration:
 *   W(q, E) = integral p dq
 *           = (E/omega) * arcsin(q/a) + (q/2) * sqrt(a^2 - q^2)
 *   where a = sqrt(2E/(m omega^2)) is the turning point.
 *
 * @param q_vals   Coordinate values for evaluation [n_points]
 * @param W_vals   Output: W(q, E) values [n_points] (NaN outside classical region)
 */
void harmonic_oscillator_HJ(double m, double omega, double E,
                              const double *q_vals, size_t n_points,
                              double *W_vals);

/**
 * Solve the time-dependent Hamilton-Jacobi equation.
 *
 *   H(q, dS/dq, t) + dS/dt = 0
 *
 * For time-independent H: S(q, alpha, t) = W(q, alpha) - E(alpha)·t
 * where alpha are separation constants (conserved momenta).
 *
 * The complete integral S(q, alpha, t) generates a canonical transformation
 * to constant new momenta alpha. The new coordinates beta = dS/dalpha
 * are also constant: beta(t) = beta_0 + omega(alpha)·t.
 *
 * @param H_func    H(q, p, t) with explicit time dependence
 * @param q0        Initial coordinate
 * @param t_end     Integration end time
 * @param dt        Time step
 * @param S_vals    Output: S(q, t) along the trajectory
 * @param n_steps   Maximum number of steps
 */
void time_dependent_HJ(double (*H_func)(const double *q, const double *p, double t, size_t n),
                        const double *q0, const double *p0, size_t n,
                        double t_end, double dt,
                        double *S_vals, size_t *n_steps_out, size_t max_steps);

/* ================================================================
 * Separation of variables
 * ================================================================ */

/**
 * Test whether a Hamiltonian is separable in the given coordinates.
 *
 * A Hamiltonian is separable if it can be written as:
 *   H(q, p) = sum_i H_i(q_i, p_i)
 *
 * or more generally if the H-J equation separates:
 *   S(q, alpha) = sum_i S_i(q_i, alpha)
 *
 * Condition for separability in coordinates q_i:
 *   d^2H / dq_i dp_j = 0  for all i != j
 *
 * Separability is the key property that makes a system integrable.
 * Stackel's theorem gives necessary and sufficient conditions.
 *
 * @param eps  Finite difference step for cross-derivative test
 * @return     true if cross-derivatives vanish within tolerance
 */
bool is_separable(double (*H_func)(const double *q, const double *p, size_t n),
                   size_t n, double eps);

/**
 * Check Stackel's separability condition for orthogonal coordinates.
 *
 * A system is separable in orthogonal coordinates u_i if there exists
 * a Stackel matrix Phi_{ij}(u_i) and a vector psi_i(u_i) such that:
 *
 *   (1/2) sum_i p_i^2/g_{ii} + V = E
 *
 * can be written as: sum_j Phi^{-1}_{1j}(H_j - E) = 0
 *
 * where H_j = (1/2) p_j^2 + V_j(u_j) for each coordinate.
 *
 * This is the classical condition for separation in the H-J sense.
 * All known integrable systems in classical mechanics are separable
 * in some coordinate system (e.g., ellipsoidal coordinates for the
 * Neumann problem).
 */
bool is_separable_stacked(size_t n,
                           const double *Phi_matrix,  /* n×n, row-major */
                           const double *psi_vector,  /* length n */
                           double eps);

/* ================================================================
 * Connections to quantum mechanics
 * ================================================================ */

/**
 * Compute the WKB phase integral (semiclassical approximation).
 *
 *   integral_{x1}^{x2} sqrt(2m(E - V(x))) dx
 *
 * This is proportional to Hamilton's characteristic function W(x, E).
 *
 * Bohr-Sommerfeld quantization condition (old quantum theory):
 *   J = (1/2pi) ∮ p dq = n hbar    →    ∮ p dq = 2pi n hbar = n h
 *
 * The WKB connection formula at turning points gives the corrected
 * quantization: ∮ p dq = (n + 1/2) h  (for soft turning points).
 *
 * This bridges classical Hamiltonian mechanics and quantum mechanics,
 * forming the foundation of the semiclassical (WKB) approximation.
 */
double wkb_phase_integral(double (*V)(double x),
                           double m, double E,
                           double x1, double x2,
                           size_t n_points);

/**
 * Compute the Bohr-Sommerfeld energy levels for a 1D potential.
 *
 * For a given quantum number n (0, 1, 2, ...), find E such that:
 *   ∮ p dq = (n + gamma) h
 *
 * where gamma = 1/2 for soft turning points (standard WKB),
 *       gamma = 0   for the old quantum theory,
 *       gamma = 3/4 for a hard wall + soft turning point.
 *
 * This uses the action function J(E) and inverts it numerically.
 *
 * Reference: Landau-Lifshitz Vol.3 §46-50, Messiah Quantum Mechanics Ch.6
 */
double bohr_sommerfeld_energy(double (*V)(double x), double m,
                               int n, double hbar,
                               double x_min, double x_max,
                               double gamma, size_t n_points);

#ifdef __cplusplus
}
#endif

#endif /* HAMILTONIAN_HAMILTON_JACOBI_H */
