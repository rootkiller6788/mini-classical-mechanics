/**
 * @file hamiltonian_action_angle.h
 * @brief Action-angle variables for integrable Hamiltonian systems.
 *
 * Reference: Goldstein Ch.10, Arnold Ch.10, Landau-Lifshitz Vol.1 Ch.7
 *            Born, The Mechanics of the Atom (1927) — historical origin
 *
 * L1 Definitions: Action J and angle theta variables
 * L2 Core Concepts: Integrability and invariant tori
 * L3 Mathematical Structures: Arnold-Liouville theorem, homology cycles
 * L6 Canonical Systems: Harmonic oscillator, pendulum, Kepler, Toda lattice
 * L7 Applications: Adiabatic invariants, EBK quantization, celestial mechanics
 *
 * University mapping:
 *   MIT 8.012 — Action-angle variables for 1D systems
 *   Princeton PHY 505 — Arnold-Liouville integrability theorem
 *   Cambridge Part II — Old quantum theory and EBK quantization
 *   ETH 402-0800 — Adiabatic invariants in classical mechanics
 */

#ifndef HAMILTONIAN_ACTION_ANGLE_H
#define HAMILTONIAN_ACTION_ANGLE_H

#include "hamiltonian_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1/L2: Action-angle variable definitions
 * ================================================================ */

/**
 * Compute the action J for a 1D system with potential U(q).
 *
 * Definition (Born-Sommerfeld):
 *   J(E) = (1/2pi) * contour_integral p dq
 *        = (1/pi)  * integral_{q_min}^{q_max} sqrt(2m(E - U(q))) dq
 *
 * The contour is taken along the energy surface H(q, p) = E.
 * For bounded motion, this forms a closed loop in phase space.
 *
 * Physical interpretation:
 *   J is the phase space area enclosed by the orbit divided by 2pi.
 *   In the old quantum theory, J = n*hbar (Bohr-Sommerfeld quantization).
 *
 * @param U         Potential energy function U(q)
 * @param m         Mass parameter
 * @param E         Energy level
 * @param q_min     Lower coordinate bound for integration
 * @param q_max     Upper coordinate bound for integration
 * @param n_points  Number of quadrature points
 * @return          Action J(E)
 */
double numerical_action_1d(double (*U)(double q),
                            double m, double E,
                            double q_min, double q_max,
                            size_t n_points);

/**
 * Compute the frequency omega(E) = dE/dJ from the action.
 *
 * For 1D systems: omega(E) = 2pi / T(E) where T(E) is the period.
 * Using the relation dJ/dE = T/(2pi), we have omega = dE/dJ.
 *
 * Numerically: omega(E) ≈ Delta_E / (J(E+Delta_E/2) - J(E-Delta_E/2))
 *
 * @param dE  Finite difference step for numerical derivative
 */
double numerical_frequency_1d(double (*U)(double q),
                               double m, double E,
                               double q_min, double q_max,
                               size_t n_points, double dE);

/**
 * Compute the canonical transformation (q, p) -> (theta, J) for 1D.
 *
 * For a given energy E, the angle variable evolves as theta = omega*t + theta0.
 * At any point on the orbit:
 *   theta(q, E) = (pi/J(E)) * integral_{q_min}^{q} p(q', E) dq'
 *   (for the upper half of the orbit; sign-adjusted for the lower half)
 */
void qp_to_action_angle_1d(double (*U)(double q), double m,
                            double q, double p,
                            double *J, double *theta,
                            size_t n_points);

/* ================================================================
 * L6: Canonical systems — specific action-angle solutions
 * ================================================================ */

/**
 * Harmonic oscillator action-angle variables.
 *
 * H = p^2/(2m) + (1/2) m omega_0^2 q^2
 *
 * Analytic solution:
 *   J = E / omega_0
 *   theta = arctan(m omega_0 q / p)
 *
 * Inverse: q = sqrt(2J/(m omega_0)) * sin(theta)
 *          p = sqrt(2J m omega_0) * cos(theta)
 */
void harmonic_oscillator_action_angle(double m, double omega0,
                                       double q, double p,
                                       double *J, double *theta);

void harmonic_oscillator_from_action_angle(double m, double omega0,
                                            double J, double theta,
                                            double *q, double *p);

/**
 * Pendulum action-angle variables.
 *
 * H = p^2/(2m L^2) - m g L cos(theta)
 *
 * Three regimes:
 *   E < -m g L:  no classical motion (forbidden)
 *   -m g L < E < m g L:  libration (oscillation)
 *   E > m g L:  rotation (full circle)
 *
 * The separatrix energy is E_sep = m g L.
 * Near the separatrix, the period diverges logarithmically.
 */
typedef enum { PENDULUM_LIBRATION, PENDULUM_ROTATION, PENDULUM_SEPARATRIX } PendulumRegime;

PendulumRegime pendulum_classify(double m, double L, double g, double E);

double pendulum_action(double m, double L, double g, double E, size_t n_points);

double pendulum_frequency(double m, double L, double g, double E, size_t n_points);

/**
 * Kepler problem radial action.
 *
 * For the Kepler problem H = p_r^2/(2m) + L^2/(2m r^2) - k/r,
 * the radial action is:
 *
 *   J_r = -L + k * sqrt(m / (2|E|))    (for E < 0, elliptical orbits)
 *
 * Combined with the angular action J_theta = L, the total degeneracy
 * J_r + J_theta = k * sqrt(m/(2|E|)) gives the energy:
 *
 *   E = -m k^2 / (2 (J_r + J_theta)^2)
 *
 * This degeneracy (frequencies omega_r = omega_theta) is responsible
 * for closed elliptical orbits (Bertrand's theorem).
 */
double kepler_radial_action(double m, double k, double E, double L);

/**
 * Infinite square well potential: V = 0 for |q| < a, V = infinity otherwise.
 *
 * Analytic action: J = (2a/pi) * sqrt(2mE)
 * Frequency: omega = pi * sqrt(2E/m) / a
 */
double infinite_well_action(double m, double a, double E);

double infinite_well_frequency(double m, double a, double E);

/**
 * Poschl-Teller potential: V(q) = -V0 * sech^2(alpha q).
 *
 * This is a reflectionless potential important in soliton theory.
 * Analytic action: J = (sqrt(2mV0) - sqrt(-2mE)) / alpha   (for E < 0)
 */
double poschl_teller_action(double m, double V0, double alpha, double E);

/* ================================================================
 * L7: Adiabatic invariants
 * ================================================================ */

/**
 * Simulate an adiabatic invariant: slowly varying parameter lambda(t).
 *
 * Theorem: For a Hamiltonian H(q, p; lambda(t)) with lambda varying
 * slowly compared to the orbital period (T * dlambda/dt << 1),
 * the action J is an adiabatic invariant:
 *
 *   J(t) = (1/2pi) * contour_integral p dq ≈ constant
 *
 * The error is exponentially small: Delta_J ~ exp(-1/epsilon) where
 * epsilon = T * dlambda/dt.
 *
 * Historical significance:
 *   — Einstein (1911): adiabatic hypothesis for old quantum theory
 *   — Landau & Lifshitz: adiabatic invariants in plasma physics
 *   — Hannay angle: geometric phase associated with adiabatic cycles
 *
 * Reference: Landau-Lifshitz Vol.1 §49, Arnold Ch.10
 */
void adiabatic_invariant_simulation(
    double (*H_func)(double q, double p, double lambda),
    void (*grad_H_func)(double q, double p, double lambda,
                         double *dH_dq, double *dH_dp),
    double q0, double p0,
    double (*lambda_func)(double t),
    double t_end, double dt,
    double *J_history, size_t *n_history,
    size_t max_history);

/* ================================================================
 * L8: KAM theory and resonances
 * ================================================================ */

/**
 * Check the Diophantine condition for a frequency vector.
 *
 * KAM theorem: For a nearly-integrable Hamiltonian H = H0(J) + eps H1(J,theta),
 * invariant tori with frequencies omega satisfying the Diophantine condition
 * survive for sufficiently small eps.
 *
 * Diophantine condition:
 *   |k · omega| ≥ gamma / |k|^tau    for all integer vectors k != 0
 *
 * where gamma > 0 and tau > n-1 (n = number of DOF).
 *
 * Tori with rational frequency ratios (resonances) are destroyed first.
 *
 * @param k_max  Maximum |k| to check
 * @param gamma  Diophantine constant
 * @param tau    Diophantine exponent (typically n-1)
 */
bool is_diophantine(const double *omega, size_t n,
                     int k_max, double gamma, double tau);

/**
 * Check if a frequency vector is resonant: k·omega ≈ 0.
 *
 * Resonance condition: there exists an integer vector k != 0 such that
 * |k·omega| < tol.
 *
 * Resonances cause the breakdown of invariant tori and the formation
 * of island chains (Poincare-Birkhoff theorem).
 */
bool is_resonant(const double *omega, size_t n,
                  const int *k, double tol);

/**
 * Compute the width of a resonance island (pendulum approximation).
 *
 * For a single resonance k·omega = 0, the dynamics near the resonance
 * is approximately a pendulum with half-width:
 *
 *   Delta_omega ≈ 2 * sqrt(eps * |H1_k|)
 *
 * where H1_k is the Fourier amplitude of the perturbation at wavevector k.
 *
 * Chirikov criterion: overlap of resonance islands -> global chaos.
 */
double resonance_island_width(double epsilon, double H1_amplitude);

#ifdef __cplusplus
}
#endif

#endif /* HAMILTONIAN_ACTION_ANGLE_H */
