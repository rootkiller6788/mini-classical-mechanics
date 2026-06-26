/**
 * rigid_body_tops.h — Heavy Symmetric Top, Lagrange Top, Gyroscope Theory
 *
 * Reference:
 *   Goldstein §5.7: Heavy symmetric top, precession, nutation
 *   Landau §35: Top motion in gravity
 *   MIT 8.012 Lecture 26: Gyroscopes
 *
 * The heavy symmetric top (Lagrange top) has I₁ = I₂ ≠ I₃.
 * Constants of motion:
 *   - L_z  (space-fixed z-axis angular momentum, from ∂L/∂φ̇)
 *   - L_3  (body-fixed z-axis angular momentum, from ∂L/∂ψ̇)
 *   - E    (total energy, from time-translation invariance)
 *
 * Effective 1D motion in θ:
 *   ½ I₁ θ̇² = E - V_eff(θ)
 *   V_eff(θ) = (L_z - L₃ cosθ)²/(2 I₁ sin²θ) + M g l cosθ + L₃²/(2 I₃)
 */

#ifndef RIGID_BODY_TOPS_H
#define RIGID_BODY_TOPS_H

#include "rigid_body_types.h"

/* ============================================================================
 * Heavy Symmetric Top — Effective Potential
 * ============================================================================ */

/**
 * Evaluate the effective potential for θ-motion of a heavy symmetric top.
 *
 * Goldstein (5.58):
 *   V_eff(θ) = (Lz - L3 cosθ)² / (2 I₁ sin²θ) + M g l cosθ + L3²/(2 I₃)
 *
 * @param I1   Transverse moment of inertia (I₁ = I₂)
 * @param I3   Axial moment of inertia
 * @param M    Mass of the top
 * @param g    Gravitational acceleration
 * @param l    Distance from pivot to CM
 * @param Lz   Conserved space-frame z angular momentum
 * @param L3   Conserved body-frame z angular momentum (I₃ ω₃)
 * @param theta Nutation angle (radians)
 * @return     V_eff(θ)
 */
double symmetric_top_eff_potential(double I1, double I3, double M, double g, double l,
                                   double Lz, double L3, double theta);

/**
 * Evaluate effective potential on a grid of θ values.
 *
 * @param I1, I3, M, g, l, Lz, L3  As above
 * @param n_pts    Number of grid points
 * @param thetas   Output: θ grid [n_pts] (caller allocates)
 * @param Veff     Output: V_eff(θ) [n_pts] (caller allocates)
 */
void symmetric_top_eff_potential_grid(double I1, double I3, double M, double g, double l,
                                      double Lz, double L3, int n_pts,
                                      double *thetas, double *Veff);

/**
 * Find turning points of θ motion: roots of E - V_eff(θ) = 0.
 * These are the bounds of nutation (θ_min, θ_max).
 *
 * @param I1, I3, M, g, l, Lz, L3  Parameters
 * @param E      Total energy
 * @param n_pts  Number of search grid points
 * @param n_tp_out  Output: number of turning points found (0, 1, or 2)
 * @param tp_out    Output: turning point θ values [2] (caller allocates)
 */
void nutation_turning_points(double I1, double I3, double M, double g, double l,
                             double Lz, double L3, double E, int n_pts,
                             int *n_tp_out, double tp_out[2]);


/* ============================================================================
 * Heavy Symmetric Top — Angular Rates
 * ============================================================================ */

/**
 * Compute dθ/dt (nutation rate) at a given θ.
 *
 * From energy conservation: ½ I₁ θ̇² = E - V_eff(θ)
 * ⇒ θ̇ = ±√(2(E - V_eff(θ))/I₁)     (Goldstein 5.55)
 *
 * Returns the positive root (magnitude). Sign depends on phase.
 */
double nutation_rate(double I1, double I3, double M, double g, double l,
                     double E, double Lz, double L3, double theta);

/**
 * Compute precession rate φ̇ at a given θ.
 *
 * Goldstein (5.57): φ̇ = (Lz - L3 cosθ) / (I₁ sin²θ)
 */
double precession_rate_top(double I1, double Lz, double L3, double theta);

/**
 * Compute spin rate ψ̇ at a given θ.
 *
 * ψ̇ = L3/I₃ - φ̇ cosθ = L3/I₃ - (Lz - L3 cosθ) cosθ / (I₁ sin²θ)
 */
double spin_rate_top(double I1, double I3, double Lz, double L3, double theta);


/* ============================================================================
 * Steady Precession
 *
 * Steady precession: θ = const (θ̇ = 0, i.e., E = V_eff(θ₀) and dV_eff/dθ|_θ₀ = 0)
 * ============================================================================ */

/**
 * Find steady precession angle θ₀ by solving dV_eff/dθ = 0 via bisection.
 *
 * Goldstein (5.59): The condition for steady precession is that θ is at an
 * extremum of V_eff.
 *
 * @param I1, I3, M, g, l, Lz, L3  Parameters
 * @param bracket  Search bracket [θ_low, θ_high]
 * @param theta0   Output: steady precession angle
 * @return         0 if root found, -1 if no sign change in bracket
 */
int steady_precession_angle(double I1, double I3, double M, double g, double l,
                            double Lz, double L3, const double bracket[2],
                            double *theta0);

/**
 * Steady precession rate at angle θ₀.
 */
double steady_precession_rate(double I1, double Lz, double L3, double theta0);

/**
 * Fast top approximation: when ω₃ is large, the gyroscopic effect dominates gravity.
 *
 * φ̇ ≈ M g l / (I₃ ω₃)     (Goldstein 5.66)
 *
 * This is the "gyroscopic approximation" — precession is slow and regular.
 */
double fast_top_precession_rate(double M, double g, double l, double I3, double omega3);

/**
 * Slow precession approximation: when Lz ≈ L3 cosθ.
 *
 * φ̇_slow ≈ M g l / L3
 */
double slow_precession_rate(double M, double g, double l, double L3);


/* ============================================================================
 * Sleeping Top Stability
 *
 * A "sleeping top" is one spinning perfectly upright (θ = 0).
 * Goldstein (5.63): The sleeping top is stable against small perturbations iff
 *   ω₃² > 4 M g l I₁ / I₃²
 *
 * When ω₃ falls below this critical value, the top "wakes up" and begins to nutate.
 * ============================================================================ */

/**
 * Analyze sleeping top stability.
 *
 * @param M, g, l, I1, I3  Parameters
 * @param omega3           Axial spin angular velocity
 * @param omega_crit_out   Output: critical spin rate
 * @return                 1 if stable (sleeping), 0 if unstable (wakes up)
 */
int sleeping_top_is_stable(double M, double g, double l, double I1, double I3,
                           double omega3, double *omega_crit_out);

/**
 * Compute the critical spin rate for sleeping top stability.
 * ω₃_crit = √(4 M g l I₁) / I₃
 */
double sleeping_top_critical_spin(double M, double g, double l, double I1, double I3);


/* ============================================================================
 * Lagrange Top — Full Equations of Motion (1st-order ODE system)
 *
 * State vector: y = [θ, φ, ψ, p_θ, p_φ, p_ψ]
 * where p_θ = I₁ θ̇, p_φ = Lz (const), p_ψ = L3 (const)
 *
 * Hamilton's equations for the heavy symmetric top (Goldstein §5.7).
 * ============================================================================ */

/**
 * Right-hand side of the Lagrange top ODE system.
 * Computes dy/dt from state y.
 *
 * @param t      Time (autonomous system, but included for ODE solvers)
 * @param y      State vector [θ, φ, ψ, p_θ, p_φ, p_ψ]
 * @param dydt   Output: time derivatives
 * @param params Pointer to struct {I1, I3, M, g, l}
 */
void lagrange_top_ode(double t, const double y[], double dydt[], const void *params);

/** Parameters for Lagrange top ODE */
typedef struct {
    double I1, I3;   /**< Transverse and axial moments of inertia */
    double M;        /**< Mass */
    double g;        /**< Gravitational acceleration */
    double l;        /**< Distance from pivot to CM */
} LagrangeTopParams;

/**
 * Simulate Lagrange top dynamics using RK4 integration.
 *
 * @param params     Top parameters
 * @param theta0     Initial nutation angle
 * @param phi0       Initial precession angle
 * @param psi0       Initial spin angle
 * @param theta_dot0 Initial nutation angular velocity
 * @param Lz         Conserved space-frame z angular momentum
 * @param L3         Conserved body-frame z angular momentum
 * @param t_end      Total simulation time
 * @param dt         Time step
 * @param n_steps_out Output: number of time steps
 * @param times_out   Output: time array [n_steps] (caller allocates)
 * @param thetas_out, phis_out, psis_out  Output: angle histories (caller allocates each)
 *
 * @return 0 on success
 */
int simulate_lagrange_top(const LagrangeTopParams *params,
                          double theta0, double phi0, double psi0,
                          double theta_dot0, double Lz, double L3,
                          double t_end, double dt,
                          int *n_steps_out, double *times_out,
                          double *thetas_out, double *phis_out, double *psis_out);


/* ============================================================================
 * Nutation Analysis
 * ============================================================================ */

/**
 * Small-amplitude nutation period.
 * For a nearly-vertical top with small θ oscillations:
 *   T_nut ≈ 2π √(I₁ / (M g l))
 */
double small_nutation_period(double I1, double M, double g, double l);

/**
 * Classify precession type from nutation amplitude range.
 *
 * Goldstein Fig.5.7 identifies three types:
 *   - MONOTONIC: φ̇ never changes sign, precession is uni-directional
 *   - LOOPING:   φ̇ changes sign, producing loops in the θ-φ projection
 *   - CUSPED:    φ̇ = 0 at turning points, producing cusps
 *
 * @param I1, Lz, L3  Parameters
 * @param theta_min, theta_max  Nutation bounds
 * @return            0=MONOTONIC, 1=LOOPING, 2=CUSPED
 */
typedef enum {
    PRECESSION_MONOTONIC = 0,
    PRECESSION_LOOPING   = 1,
    PRECESSION_CUSPED    = 2
} PrecessionType;

PrecessionType classify_precession(double I1, double Lz, double L3,
                                   double theta_min, double theta_max);


/* ============================================================================
 * Gyroscope Principles
 * ============================================================================ */

/**
 * Compute gyroscopic torque: N = ω_forced × L_spin.
 *
 * When a spinning body is forced to rotate about an axis not aligned with
 * its spin axis, it experiences a torque perpendicular to both axes (Coriolis effect).
 *
 * @param L_spin       Spin angular momentum vector
 * @param omega_forced Forced precession angular velocity vector
 * @param torque       Output: gyroscopic reaction torque
 */
void gyroscopic_torque(const double L_spin[3], const double omega_forced[3], double torque[3]);

/**
 * Compute gyroscope precession rate under constant external torque.
 *
 * For a single-degree-of-freedom gyroscope:
 *   Ω = N / L_spin     (steady precession rate)
 *
 * @param N       Magnitude of external torque (perpendicular to spin axis)
 * @param L_spin  Magnitude of spin angular momentum
 * @return        Steady precession angular rate
 */
double gyroscope_precession_rate(double N, double L_spin);

/**
 * Analyze a 2-axis gyroscope response to base rotation.
 * Models the full coupled dynamics of a gyroscope mounted on a rotating platform.
 *
 * @param L_spin       Spin angular momentum (magnitude)
 * @param omega_base   Base rotation angular velocity vector
 * @param I_gimbal     Gimbal inertia about output axis
 * @param K_spring     Torsion spring constant (if any, 0 for free gyro)
 * @param D_damping    Viscous damping coefficient
 * @param dt           Time step
 * @param n_steps      Number of steps
 * @param theta_out    Output: gimbal deflection angle history [n_steps]
 * @param omega_out    Output: gimbal angular velocity history [n_steps]
 */
void gyroscope_2axis_simulate(double L_spin, const double omega_base[3],
                              double I_gimbal, double K_spring, double D_damping,
                              double dt, int n_steps,
                              double *theta_out, double *omega_out);

/**
 * Compute nutation frequency of a free gyroscope.
 *
 * When perturbed, a free gyroscope nutates at frequency:
 *   ω_nut = L / I_transverse
 *
 * @param L_spin     Spin angular momentum
 * @param I_trans    Transverse moment of inertia
 * @return           Nutation angular frequency
 */
double free_gyro_nutation_frequency(double L_spin, double I_trans);

#endif /* RIGID_BODY_TOPS_H */
