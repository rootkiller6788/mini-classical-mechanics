/**
 * rigid_body_euler.h — Euler Equations of Rigid Body Dynamics
 *
 * Reference:
 *   Goldstein §5.4-5.5: Torque-free Euler equations, Poinsot construction
 *   Landau §36-37: Euler's equations, stability of rotation
 *   MIT 8.012 Lecture 25-26: Rigid body dynamics
 *
 * The Euler equations describe rotation of a rigid body in the body frame:
 *   I1 ω̇₁ - (I2 - I3) ω₂ ω₃ = N1
 *   I2 ω̇₂ - (I3 - I1) ω₃ ω₁ = N2
 *   I3 ω̇₃ - (I1 - I2) ω₁ ω₂ = N3
 *
 * For torque-free motion (N=0), two constants of motion exist:
 *   T = ½ Σ I_i ω_i²   (kinetic energy)
 *   L² = Σ (I_i ω_i)²  (angular momentum squared)
 */

#ifndef RIGID_BODY_EULER_H
#define RIGID_BODY_EULER_H

#include "rigid_body_types.h"

/* ============================================================================
 * Euler Equations — Right-Hand Side Evaluation
 * ============================================================================ */

/**
 * Torque-free Euler equations: dω/dt = f(ω) in principal axis frame.
 *
 * Goldstein (5.44) with I1, I2, I3 along principal axes:
 *   dω₁/dt = (I2 - I3)/I1 · ω₂ ω₃
 *   dω₂/dt = (I3 - I1)/I2 · ω₃ ω₁
 *   dω₃/dt = (I1 - I2)/I3 · ω₁ ω₂
 *
 * @param I       Principal moments [I1, I2, I3] (must be diagonal)
 * @param omega   Current angular velocity [ω₁, ω₂, ω₃] (body frame)
 * @param domega  Output angular acceleration [dω₁/dt, dω₂/dt, dω₃/dt]
 */
void euler_free_derivative(const double I[3], const double omega[3], double domega[3]);

/**
 * Euler equations with external torque: dω/dt = f(ω) + I⁻¹·N.
 *
 * Goldstein (5.41-5.43) generalized:
 *   dω₁/dt = ((I2 - I3) ω₂ ω₃ + N1) / I1
 *   dω₂/dt = ((I3 - I1) ω₃ ω₁ + N2) / I2
 *   dω₃/dt = ((I1 - I2) ω₁ ω₂ + N3) / I3
 *
 * @param I       Principal moments [I1, I2, I3]
 * @param omega   Current angular velocity (body frame)
 * @param torque  External torque vector (body frame components)
 * @param domega  Output angular acceleration
 */
void euler_with_torque(const double I[3], const double omega[3],
                       const double torque[3], double domega[3]);

/**
 * General Euler equations in arbitrary (non-principal) frame.
 * Uses full inertia matrix: I · dω/dt + ω × (I·ω) = N
 *   → dω/dt = I⁻¹ · (N - ω × (I·ω))
 *
 * @param Imat   Full 3×3 inertia matrix (arbitrary orientation)
 * @param omega  Angular velocity
 * @param torque External torque
 * @param domega Output angular acceleration
 */
void euler_general_frame(const mat3 *Imat, const double omega[3],
                         const double torque[3], double domega[3]);


/* ============================================================================
 * Constants of Motion (Torque-Free Case)
 * ============================================================================ */

/**
 * Compute the two constants of motion for torque-free rigid body.
 *
 * Goldstein (5.45):
 *   T  = ½ (I1 ω₁² + I2 ω₂² + I3 ω₃²)    (kinetic energy)
 *   L² = I1² ω₁² + I2² ω₂² + I3² ω₃²      (angular momentum squared)
 *
 * @param I     Principal moments [I1, I2, I3]
 * @param omega Angular velocity
 * @param T     Output kinetic energy
 * @param L2    Output angular momentum squared |L|²
 */
void motion_constants_euler(const double I[3], const double omega[3],
                            double *T, double *L2);

/**
 * Verify that a given angular velocity satisfies the polhode constraint:
 * the angular velocity vector lies on the intersection of the kinetic energy
 * ellipsoid and the angular momentum sphere.
 *
 * @param I     Principal moments
 * @param omega Angular velocity to check
 * @param tol   Tolerance for constant drift
 * @return      1 if consistent, 0 otherwise
 */
int polhode_constraint_check(const double I[3], const double omega[3], double tol);

/**
 * Compute the Binet ellipsoid intersection constraint value:
 * f(ω) = |ω|² - 2T·(I⁻¹)_avg   (a derived invariant for polhode analysis)
 */
double binet_ellipsoid_eval(const double I[3], const double omega[3]);


/* ============================================================================
 * Numerical Integrators for Body-Frame Angular Velocity
 * ============================================================================ */

/**
 * Forward Euler step (O(dt), 1st order).
 * ω_{n+1} = ω_n + dt · f(ω_n)
 */
void euler_step_omega(const double I[3], const double omega[3],
                      double dt, double omega_next[3]);

/**
 * RK4 (Runge-Kutta 4th order) step — recommended for accuracy.
 * Classical 4-stage RK4 with local error O(dt⁵), global O(dt⁴).
 */
void rk4_step_omega(const double I[3], const double omega[3],
                     double dt, double omega_next[3]);

/**
 * Midpoint method (RK2, 2nd order) — balance of speed and stability.
 */
void midpoint_step_omega(const double I[3], const double omega[3],
                         double dt, double omega_next[3]);

/**
 * Velocity Verlet adapted for rotational dynamics (symplectic-like).
 *
 * Half-step acceleration → full-step velocity → new acceleration → half-step correction.
 * Good long-term energy behavior for Hamiltonian-like systems.
 */
void verlet_step_omega(const double I[3], const double omega[3],
                       double dt, double omega_next[3]);

/**
 * Dormand-Prince 5(4) adaptive step (embedded RK pair).
 * Takes one step with error estimation for adaptive step-size control.
 *
 * @param I        Principal moments
 * @param omega    Input/output: angular velocity (updated in place)
 * @param dt       Suggested step size
 * @param dt_next  Output: recommended next step size
 * @param error    Output: estimated local truncation error
 */
void dopri54_step_omega(const double I[3], double omega[3],
                        double dt, double *dt_next, double *error);


/* ============================================================================
 * Full Trajectory Simulation
 * ============================================================================ */

/** Integration method selector */
typedef enum {
    INTEGRATOR_EULER = 0,
    INTEGRATOR_RK4   = 1,
    INTEGRATOR_MIDPOINT = 2,
    INTEGRATOR_VERLET = 3,
    INTEGRATOR_DOPRI54 = 4
} IntegratorMethod;

/**
 * Simulate torque-free rigid body motion: body-frame angular velocity evolution.
 *
 * @param I           Principal moments [I1, I2, I3]
 * @param omega0      Initial angular velocity (body frame)
 * @param t_end       Total simulation time
 * @param dt          Time step (fixed for non-adaptive methods)
 * @param method      Integration method
 * @param n_steps_out Output: number of steps in trajectory
 * @param traj_out    Output: trajectory array [[ωx,ωy,ωz]_0, [ωx,ωy,ωz]_1, ...]
 *                    (caller allocates: max_steps = ceil(t_end/dt) + 1)
 *
 * @return 0 on success, -1 on error
 */
int simulate_free_rigid_body(const double I[3], const double omega0[3],
                             double t_end, double dt, IntegratorMethod method,
                             int *n_steps_out, double *traj_out);

/**
 * Simulate rigid body with external torque function.
 *
 * @param I           Principal moments
 * @param omega0      Initial angular velocity
 * @param torque_fn   External torque N(t, ω) → torque[3]
 * @param params      User data passed to torque_fn
 * @param t_end       Total simulation time
 * @param dt          Time step
 * @param n_steps_out Output: number of steps
 * @param times_out   Output: time array [n_steps]
 * @param omega_out   Output: angular velocity array [n_steps * 3]
 *
 * @return 0 on success, -1 on error
 */
typedef void (*torque_function)(double t, const double omega[3], double torque[3], void *params);

int simulate_rigid_body_torque(const double I[3], const double omega0[3],
                               torque_function torque_fn, void *params,
                               double t_end, double dt,
                               int *n_steps_out, double *times_out, double *omega_out);

/**
 * Simulate full rigid body including spatial orientation (Euler angles).
 * Couples angular velocity evolution (Euler eqn) with kinematic equations
 * for Euler angle evolution: d(euler)/dt = f(ω, euler).
 *
 * @param I             Principal moments
 * @param omega0        Initial angular velocity (body frame)
 * @param euler0        Initial Euler angles
 * @param t_end         Total simulation time
 * @param dt            Time step
 * @param n_steps_out   Output: number of steps
 * @param times_out     Output: time array
 * @param omega_out     Output: omega trajectory [n_steps * 3]
 * @param euler_out     Output: Euler angles trajectory [n_steps * 3]
 *
 * @return 0 on success
 */
int simulate_full_rigid_body(const double I[3], const double omega0[3],
                             const EulerAngles *euler0,
                             double t_end, double dt,
                             int *n_steps_out, double *times_out,
                             double *omega_out, double *euler_out);


/* ============================================================================
 * Motion Constant Monitoring
 * ============================================================================ */

/**
 * Monitor kinetic energy and angular momentum drift along a trajectory.
 *
 * @param I         Principal moments
 * @param n_steps   Number of time steps
 * @param traj      Angular velocity trajectory [n_steps * 3]
 * @param T_hist    Output kinetic energy history [n_steps]
 * @param L2_hist   Output |L|² history [n_steps]
 */
void monitor_motion_constants(const double I[3], int n_steps, const double *traj,
                              double *T_hist, double *L2_hist);

/**
 * Compute drift statistics from monitored constants.
 *
 * @param n_steps    Number of steps
 * @param T_hist     Kinetic energy history
 * @param L2_hist    |L|² history
 * @param T_drift    Output: max fractional drift in T
 * @param L2_drift   Output: max fractional drift in L²
 */
void constant_drift_report(int n_steps, const double *T_hist, const double *L2_hist,
                           double *T_drift, double *L2_drift);

/**
 * Compute body-to-space angular velocity transformation along trajectory.
 * ω_space(t) = R(t) · ω_body(t) where R(t) is from Euler angles.
 *
 * @param n_steps        Number of steps
 * @param omega_body     Body-frame omega trajectory [n_steps * 3]
 * @param euler_traj     Euler angle trajectory [n_steps * 3]
 * @param omega_space    Output: space-frame omega trajectory [n_steps * 3]
 */
void body_to_space_trajectory(int n_steps, const double *omega_body,
                              const double *euler_traj, double *omega_space);


/* ============================================================================
 * Poinsot Geometric Construction
 *
 * The Poinsot construction (1834) gives a geometric interpretation:
 * The inertia ellipsoid rolls without slipping on the "invariable plane"
 * perpendicular to the (constant) angular momentum vector.
 * 
 * - Polhode:  curve traced by ω on the inertia ellipsoid (body frame)
 * - Herpolhode: curve traced by ω on the invariable plane (space frame)
 * ============================================================================ */

/**
 * Compute the invariable plane normal (direction of angular momentum L in space).
 * For torque-free motion, L is constant in space frame.
 *
 * @param I       Principal moments
 * @param omega   Angular velocity in body frame (at any time)
 * @param R       Rotation matrix from body to space frame
 * @param normal  Output: unit normal to invariable plane (space frame)
 * @param dist    Output: distance of invariable plane from origin = 2T/|L|
 */
void invariable_plane(const double I[3], const double omega[3],
                      const mat3 *R, vec3 *normal, double *dist);

/**
 * Compute polhode tangent direction at a point on the inertia ellipsoid.
 * The polhode is the intersection of the energy ellipsoid and momentum sphere.
 *
 * @param I       Principal moments
 * @param omega   Angular velocity point on polhode
 * @param tangent Output: tangent direction to polhode curve
 */
void polhode_tangent(const double I[3], const double omega[3], double tangent[3]);

#endif /* RIGID_BODY_EULER_H */
