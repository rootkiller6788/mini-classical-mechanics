/**
 * rigid_body_energy.h — Rotational Energy, Angular Momentum, Stability Analysis
 *
 * Reference:
 *   Goldstein §5.1-5.2, §5.6: Rotational kinetic energy, angular momentum,
 *                             Tennis racket theorem (intermediate axis instability)
 *   Landau §36-37: Stability of rotation, Euler equations
 *   MIT 8.012 Lecture 24-26
 *
 * Key results:
 *   - Rotational KE: T = ½ ωᵀ · I · ω
 *   - Angular momentum: L = I · ω
 *   - In body frame with principal axes: T = ½(I1 ω₁² + I2 ω₂² + I3 ω₃²)
 *   - Tennis Racket Theorem: rotation about intermediate principal axis is unstable
 *   - Energy-momentum inequality: T_min = L²/(2I₁) ≤ T ≤ L²/(2I₃) = T_max
 */

#ifndef RIGID_BODY_ENERGY_H
#define RIGID_BODY_ENERGY_H

#include "rigid_body_types.h"

/* ============================================================================
 * Rotational Kinetic Energy
 * ============================================================================ */

/**
 * Rotational kinetic energy using general inertia tensor.
 * T = ½ ωᵀ · I · ω    (Goldstein 5.9)
 */
double rotational_ke(const InertiaTensor *I, const double omega[3]);

/**
 * Rotational kinetic energy in principal axis frame (faster, no matrix multiply).
 * T = ½ (I₁ ω₁² + I₂ ω₂² + I₃ ω₃²)
 */
double rotational_ke_principal(const double I_principal[3], const double omega[3]);

/**
 * Total kinetic energy of a rigid body: translation + rotation.
 * T_total = ½ M |v_cm|² + ½ ωᵀ · I_cm · ω   (König's theorem for rigid bodies)
 *
 * @param M      Total mass
 * @param v_cm   CM velocity vector
 * @param I_cm   Inertia tensor about CM
 * @param omega  Angular velocity (body frame)
 * @param T_trans Output: translational kinetic energy
 * @param T_rot   Output: rotational kinetic energy
 * @return Total kinetic energy T_trans + T_rot
 */
double total_kinetic_energy(double M, vec3 v_cm, const InertiaTensor *I_cm,
                            const double omega[3], double *T_trans, double *T_rot);


/* ============================================================================
 * Angular Momentum
 * ============================================================================ */

/**
 * Compute angular momentum from inertia tensor and angular velocity.
 * L = I · ω = [I_xx ω_x + I_xy ω_y + I_xz ω_z, ...]   (Goldstein 5.10)
 *
 * @param I      Inertia tensor
 * @param omega  Angular velocity (body frame)
 * @param L      Output angular momentum (body frame)
 */
void angular_momentum_from_inertia(const InertiaTensor *I, const double omega[3], double L[3]);

/**
 * Angular momentum in principal axis frame.
 * L = (I₁ ω₁, I₂ ω₂, I₃ ω₃) — diagonal case
 */
void angular_momentum_principal(const double I_principal[3], const double omega[3], double L[3]);

/**
 * Total angular momentum: orbital + spin.
 * L_total = r_cm × (M v_cm) + I_cm · ω
 *
 * @param r_cm     CM position vector
 * @param M        Total mass
 * @param v_cm     CM velocity
 * @param I_cm     Inertia tensor about CM
 * @param omega    Angular velocity
 * @param L_total  Output: total angular momentum
 * @param L_orbital Output: orbital angular momentum
 * @param L_spin    Output: spin angular momentum
 */
void total_angular_momentum(vec3 r_cm, double M, vec3 v_cm,
                           const InertiaTensor *I_cm, const double omega[3],
                           double L_total[3], double L_orbital[3], double L_spin[3]);

/**
 * Compute angular velocity from angular momentum and inertia tensor.
 * ω = I⁻¹ · L   (inverse inertia tensor)
 *
 * For principal frame: ω_i = L_i / I_i
 *
 * @param I_principal  Principal moments
 * @param L            Angular momentum (principal frame)
 * @param omega        Output: angular velocity
 */
void angular_velocity_from_momentum(const double I_principal[3], const double L[3], double omega[3]);


/* ============================================================================
 * Intermediate Axis Theorem (Tennis Racket Theorem / Dzhanibekov Effect)
 *
 * Goldstein §5.6: For a rigid body with distinct principal moments I₁ > I₂ > I₃,
 * rotation about the intermediate axis (I₂) is unstable under small perturbations.
 *
 * Linearized analysis (Landau §37):
 *   Let ω = ω₀ e_k + δω. The perturbation evolves as:
 *     δω_i ∝ exp(λ t) where λ² = ω₀² (I_k - I_j)(I_i - I_k) / (I_i I_j)
 *   - For k=1 (max): λ² < 0 → oscillatory → stable
 *   - For k=2 (mid): λ² > 0 → exponential → unstable
 *   - For k=3 (min): λ² < 0 → oscillatory → stable
 * ============================================================================ */

/** Stability result type */
typedef enum {
    STABILITY_STABLE = 0,     /**< Small perturbations oscillate, no growth */
    STABILITY_UNSTABLE = 1,   /**< Small perturbations grow exponentially */
    STABILITY_DEGENERATE = 2  /**< Symmetric case, marginal stability */
} StabilityResult;

/**
 * Analyze linear stability of rotation about a specified principal axis.
 *
 * @param I_principal  Principal moments [I1, I2, I3], I1 ≥ I2 ≥ I3
 * @param axis         Which principal axis: 0→I1, 1→I2, 2→I3
 * @param omega_mag    Magnitude of angular velocity about that axis
 * @param growth_rate  Output: exponential growth rate (0 if stable, λ if unstable)
 * @return             STABILITY_STABLE, STABILITY_UNSTABLE, or STABILITY_DEGENERATE
 */
StabilityResult axis_stability_analysis(const double I_principal[3], int axis,
                                        double omega_mag, double *growth_rate);

/**
 * Full stability characterization of all three principal axes.
 *
 * @param I_principal  Principal moments
 * @param results      Output: stability for axes 0,1,2
 */
void full_stability_analysis(const double I_principal[3], StabilityResult results[3]);

/**
 * Estimate the flipping period for intermediate axis rotation.
 *
 * When a body rotates about its intermediate axis, small perturbations cause
 * periodic 180° flips. The half-period (time between flips) is approximately:
 *   T_flip ≈ ln(1/δω₀) / λ
 * where δω₀ is the initial perturbation magnitude and λ is the growth rate.
 *
 * @param I_principal   Principal moments
 * @param omega_mag     Nominal angular velocity about intermediate axis
 * @param delta_omega0  Initial perturbation magnitude (fraction of omega_mag)
 * @return              Estimated half-flip period (seconds)
 */
double flipping_period_estimate(const double I_principal[3], double omega_mag,
                                double delta_omega0);

/**
 * Simulate the Dzhanibekov effect by integrating Euler equations.
 * Starting near the intermediate axis, the body undergoes periodic flips.
 *
 * @param I_principal  Principal moments
 * @param omega0       Initial angular velocity (near intermediate axis)
 * @param t_end        Total simulation time
 * @param dt           Time step
 * @param n_steps_out  Output: number of steps
 * @param traj_out     Output: ω trajectory [n_steps * 3]
 * @param flip_times   Output: times of flip events (caller allocates ~20)
 * @param n_flips_out  Output: number of flips detected
 */
void dzhanibekov_simulation(const double I_principal[3], const double omega0[3],
                            double t_end, double dt,
                            int *n_steps_out, double *traj_out,
                            double *flip_times, int *n_flips_out);


/* ============================================================================
 * Energy-Momentum Geometry on the Angular Momentum Sphere
 *
 * For a given |L|² = L², the kinetic energy for all possible ω on the
 * angular momentum sphere (in ω-space) takes values in [T_min, T_max]:
 *   T_min = L² / (2 I₁)    (rotation about axis of largest inertia)
 *   T_max = L² / (2 I₃)    (rotation about axis of smallest inertia)
 * ============================================================================ */

/**
 * Compute minimum and maximum possible kinetic energy for given |L|².
 *
 * @param I_principal  Principal moments
 * @param L2           Angular momentum squared |L|²
 * @param T_min        Output: minimum possible KE
 * @param T_max        Output: maximum possible KE
 */
void energy_extrema_for_L2(const double I_principal[3], double L2,
                           double *T_min, double *T_max);

/**
 * Sample kinetic energy values on the angular momentum sphere.
 * For visualization of the Binet ellipsoid intersection with momentum sphere.
 *
 * @param I_principal  Principal moments
 * @param L2           Angular momentum squared |L|²
 * @param n_theta      Polar grid points
 * @param n_phi        Azimuthal grid points
 * @param T_map        Output: KE values [n_theta * n_phi] (caller allocates)
 * @param omega1,omega2,omega3  Output: sampled ω components (same size as T_map)
 */
void energy_on_momentum_sphere(const double I_principal[3], double L2,
                               int n_theta, int n_phi,
                               double *T_map, double *omega1, double *omega2, double *omega3);

/**
 * Compute the separatrix values: energies that separate regions of different
 * topological type on the energy-momentum diagram.
 *
 * Goldstein Fig.5.5: The separatrix passes through the intermediate axis,
 * where the polhode is self-intersecting.
 *
 * @param I_principal  Principal moments
 * @param L2           Angular momentum squared
 * @param T_sep        Output: separatrix energy = L²/(2 I₂)
 * @return             0 if separatrix is accessible, -1 if not (L²/(2I₂) outside [T_min, T_max])
 */
int separatrix_energy(const double I_principal[3], double L2, double *T_sep);


/* ============================================================================
 * Torque, Power, and Impulse
 * ============================================================================ */

/**
 * Instantaneous power of a torque: P = N · ω.
 * This is the rate at which the torque does work on the rigid body.
 */
double torque_power(const double torque[3], const double omega[3]);

/**
 * Angular impulse (torque integrated over time interval).
 * ΔL = ∫ N dt ≈ N · Δt  (for approximately constant torque)
 */
void torque_impulse(const double torque[3], double dt, double dL[3]);

/**
 * Work done by torque over a trajectory.
 * W = ∫ N(t) · ω(t) dt   approximated via trapezoidal rule.
 *
 * @param n_steps   Number of time steps
 * @param times     Time array [n_steps]
 * @param torques   Torque array [n_steps * 3]
 * @param omegas    Angular velocity array [n_steps * 3]
 * @return          Total work done
 */
double torque_work(int n_steps, const double *times, const double *torques, const double *omegas);

/**
 * Compute the angular momentum change from torque history.
 * ΔL = Σ N_i Δt_i (Euler-forward approximation).
 *
 * @param n_steps    Number of time intervals
 * @param torques    Torque history [n_steps * 3]
 * @param dt         Time step (constant)
 * @param delta_L    Output: total angular momentum change
 */
void torque_angular_impulse_total(int n_steps, const double *torques, double dt, double delta_L[3]);

#endif /* RIGID_BODY_ENERGY_H */
