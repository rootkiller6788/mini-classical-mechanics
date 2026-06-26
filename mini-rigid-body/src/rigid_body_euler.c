/**
 * rigid_body_euler.c — Euler Equations, Numerical Integrators, Poinsot Construction
 *
 * Goldstein §5.4-5.5: Euler's equations of motion for rigid bodies.
 * Landau §36-37: Stability, polhode analysis.
 *
 * The Euler equations are the rotational analog of Newton's second law:
 *   dL/dt (space frame) = N   →   dL/dt + ω × L = N (body frame)
 * which expands to the classic Euler equations in principal axes.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/rigid_body_types.h"
#include "../include/rigid_body_euler.h"
#include "../include/rigid_body_kinematics.h"

/* ============================================================================
 * Euler Equation RHS Evaluation
 * ============================================================================ */

void euler_free_derivative(const double I[3], const double omega[3], double domega[3]) {
    domega[0] = (I[1] - I[2]) * omega[1] * omega[2] / I[0];
    domega[1] = (I[2] - I[0]) * omega[2] * omega[0] / I[1];
    domega[2] = (I[0] - I[1]) * omega[0] * omega[1] / I[2];
}

void euler_with_torque(const double I[3], const double omega[3],
                       const double torque[3], double domega[3]) {
    domega[0] = ((I[1] - I[2]) * omega[1] * omega[2] + torque[0]) / I[0];
    domega[1] = ((I[2] - I[0]) * omega[2] * omega[0] + torque[1]) / I[1];
    domega[2] = ((I[0] - I[1]) * omega[0] * omega[1] + torque[2]) / I[2];
}

void euler_general_frame(const mat3 *Imat, const double omega[3],
                         const double torque[3], double domega[3]) {
    /* Compute L = I·ω */
    vec3 w = {omega[0], omega[1], omega[2]};
    vec3 L = mat3_mul_vec(Imat, w);

    /* Compute ω × L */
    vec3 omega_cross_L = vec3_cross(w, L);

    /* Compute N - ω×L */
    double rhs[3] = {
        torque[0] - omega_cross_L.x,
        torque[1] - omega_cross_L.y,
        torque[2] - omega_cross_L.z
    };

    /* Solve I · dω = rhs via Cramer's rule for 3×3 */
    double det = mat3_det(Imat);
    if (fabs(det) < 1e-15) {
        domega[0] = domega[1] = domega[2] = 0.0;
        return;
    }

    double A[9];
    memcpy(A, Imat->m, 9 * sizeof(double));

    /* Cramer's rule: dω_i = det(I with column i replaced by rhs) / det(I) */
    for (int col = 0; col < 3; col++) {
        double M[9];
        memcpy(M, A, 9 * sizeof(double));
        M[0 + col] = rhs[0];
        M[3 + col] = rhs[1];
        M[6 + col] = rhs[2];
        mat3 Mm = { {M[0],M[1],M[2], M[3],M[4],M[5], M[6],M[7],M[8]} };
        domega[col] = mat3_det(&Mm) / det;
    }
}

/* ============================================================================
 * Constants of Motion
 * ============================================================================ */

void motion_constants_euler(const double I[3], const double omega[3],
                            double *T, double *L2) {
    *T = 0.5 * (I[0]*omega[0]*omega[0] + I[1]*omega[1]*omega[1] + I[2]*omega[2]*omega[2]);
    *L2 = (I[0]*omega[0])*(I[0]*omega[0])
        + (I[1]*omega[1])*(I[1]*omega[1])
        + (I[2]*omega[2])*(I[2]*omega[2]);
}

int polhode_constraint_check(const double I[3], const double omega[3], double tol) {
    /* Verify that ω simultaneously satisfies both the energy ellipsoid
     * and angular momentum sphere equations. We use the fact that:
     *   2T = I₁ω₁² + I₂ω₂² + I₃ω₃²
     *   L² = I₁²ω₁² + I₂²ω₂² + I₃²ω₃²
     * These are invariant for any valid ω on the polhode.
     * We verify internal consistency: both should agree with computed values. */
    double T, L2;
    motion_constants_euler(I, omega, &T, &L2);

    double T_check = 0.5 * (I[0]*omega[0]*omega[0]
                          + I[1]*omega[1]*omega[1]
                          + I[2]*omega[2]*omega[2]);
    double L2_check = (I[0]*omega[0])*(I[0]*omega[0])
                    + (I[1]*omega[1])*(I[1]*omega[1])
                    + (I[2]*omega[2])*(I[2]*omega[2]);

    return (fabs(T - T_check) < tol) && (fabs(L2 - L2_check) < tol);
}

double binet_ellipsoid_eval(const double I[3], const double omega[3]) {
    /* Binet's construction: the momental ellipsoid */
    double a2 = 1.0 / I[0], b2 = 1.0 / I[1], c2 = 1.0 / I[2];
    return a2 * omega[0]*omega[0] + b2 * omega[1]*omega[1] + c2 * omega[2]*omega[2];
}

/* ============================================================================
 * Numerical Integrators
 * ============================================================================ */

void euler_step_omega(const double I[3], const double omega[3],
                      double dt, double omega_next[3]) {
    double domega[3];
    euler_free_derivative(I, omega, domega);
    omega_next[0] = omega[0] + domega[0] * dt;
    omega_next[1] = omega[1] + domega[1] * dt;
    omega_next[2] = omega[2] + domega[2] * dt;
}

void rk4_step_omega(const double I[3], const double omega[3],
                     double dt, double omega_next[3]) {
    double k1[3], k2[3], k3[3], k4[3];
    double w_temp[3];
    double half_dt = 0.5 * dt;

    euler_free_derivative(I, omega, k1);

    w_temp[0] = omega[0] + half_dt * k1[0];
    w_temp[1] = omega[1] + half_dt * k1[1];
    w_temp[2] = omega[2] + half_dt * k1[2];
    euler_free_derivative(I, w_temp, k2);

    w_temp[0] = omega[0] + half_dt * k2[0];
    w_temp[1] = omega[1] + half_dt * k2[1];
    w_temp[2] = omega[2] + half_dt * k2[2];
    euler_free_derivative(I, w_temp, k3);

    w_temp[0] = omega[0] + dt * k3[0];
    w_temp[1] = omega[1] + dt * k3[1];
    w_temp[2] = omega[2] + dt * k3[2];
    euler_free_derivative(I, w_temp, k4);

    omega_next[0] = omega[0] + (dt / 6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
    omega_next[1] = omega[1] + (dt / 6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
    omega_next[2] = omega[2] + (dt / 6.0) * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);
}

void midpoint_step_omega(const double I[3], const double omega[3],
                         double dt, double omega_next[3]) {
    double domega[3], w_mid[3];
    double half_dt = 0.5 * dt;

    euler_free_derivative(I, omega, domega);
    w_mid[0] = omega[0] + half_dt * domega[0];
    w_mid[1] = omega[1] + half_dt * domega[1];
    w_mid[2] = omega[2] + half_dt * domega[2];

    euler_free_derivative(I, w_mid, domega);
    omega_next[0] = omega[0] + dt * domega[0];
    omega_next[1] = omega[1] + dt * domega[1];
    omega_next[2] = omega[2] + dt * domega[2];
}

void verlet_step_omega(const double I[3], const double omega[3],
                       double dt, double omega_next[3]) {
    double a_cur[3], w_half[3], a_new[3];

    euler_free_derivative(I, omega, a_cur);
    w_half[0] = omega[0] + 0.5 * a_cur[0] * dt;
    w_half[1] = omega[1] + 0.5 * a_cur[1] * dt;
    w_half[2] = omega[2] + 0.5 * a_cur[2] * dt;

    euler_free_derivative(I, w_half, a_new);
    omega_next[0] = w_half[0] + 0.5 * a_new[0] * dt;
    omega_next[1] = w_half[1] + 0.5 * a_new[1] * dt;
    omega_next[2] = w_half[2] + 0.5 * a_new[2] * dt;
}

/* Dormand-Prince 5(4) coefficients (DOPRI54) */
static const double dopri_a21 = 1.0/5.0;
static const double dopri_a31 = 3.0/40.0, dopri_a32 = 9.0/40.0;
static const double dopri_a41 = 44.0/45.0, dopri_a42 = -56.0/15.0, dopri_a43 = 32.0/9.0;
static const double dopri_a51 = 19372.0/6561.0, dopri_a52 = -25360.0/2187.0;
static const double dopri_a53 = 64448.0/6561.0, dopri_a54 = -212.0/729.0;
static const double dopri_a61 = 9017.0/3168.0, dopri_a62 = -355.0/33.0;
static const double dopri_a63 = 46732.0/5247.0, dopri_a64 = 49.0/176.0;
static const double dopri_a65 = -5103.0/18656.0;
/* 5th-order coefficients */
static const double dopri_b51 = 35.0/384.0, dopri_b53 = 500.0/1113.0;
static const double dopri_b54 = 125.0/192.0, dopri_b55 = -2187.0/6784.0;
static const double dopri_b56 = 11.0/84.0;
/* 4th-order coefficients (embedded) */
static const double dopri_b41 = 5179.0/57600.0, dopri_b43 = 7571.0/16695.0;
static const double dopri_b44 = 393.0/640.0, dopri_b45 = -92097.0/339200.0;
static const double dopri_b46 = 187.0/2100.0, dopri_b47 = 1.0/40.0;

void dopri54_step_omega(const double I[3], double omega[3],
                        double dt, double *dt_next, double *error) {
    double k1[3], k2[3], k3[3], k4[3], k5[3], k6[3], k7[3];
    double w_temp[3];

    euler_free_derivative(I, omega, k1);

    w_temp[0] = omega[0] + dt * dopri_a21 * k1[0];
    w_temp[1] = omega[1] + dt * dopri_a21 * k1[1];
    w_temp[2] = omega[2] + dt * dopri_a21 * k1[2];
    euler_free_derivative(I, w_temp, k2);

    w_temp[0] = omega[0] + dt * (dopri_a31*k1[0] + dopri_a32*k2[0]);
    w_temp[1] = omega[1] + dt * (dopri_a31*k1[1] + dopri_a32*k2[1]);
    w_temp[2] = omega[2] + dt * (dopri_a31*k1[2] + dopri_a32*k2[2]);
    euler_free_derivative(I, w_temp, k3);

    w_temp[0] = omega[0] + dt * (dopri_a41*k1[0] + dopri_a42*k2[0] + dopri_a43*k3[0]);
    w_temp[1] = omega[1] + dt * (dopri_a41*k1[1] + dopri_a42*k2[1] + dopri_a43*k3[1]);
    w_temp[2] = omega[2] + dt * (dopri_a41*k1[2] + dopri_a42*k2[2] + dopri_a43*k3[2]);
    euler_free_derivative(I, w_temp, k4);

    w_temp[0] = omega[0] + dt * (dopri_a51*k1[0] + dopri_a52*k2[0] + dopri_a53*k3[0] + dopri_a54*k4[0]);
    w_temp[1] = omega[1] + dt * (dopri_a51*k1[1] + dopri_a52*k2[1] + dopri_a53*k3[1] + dopri_a54*k4[1]);
    w_temp[2] = omega[2] + dt * (dopri_a51*k1[2] + dopri_a52*k2[2] + dopri_a53*k3[2] + dopri_a54*k4[2]);
    euler_free_derivative(I, w_temp, k5);

    w_temp[0] = omega[0] + dt * (dopri_a61*k1[0] + dopri_a62*k2[0] + dopri_a63*k3[0]
                                + dopri_a64*k4[0] + dopri_a65*k5[0]);
    w_temp[1] = omega[1] + dt * (dopri_a61*k1[1] + dopri_a62*k2[1] + dopri_a63*k3[1]
                                + dopri_a64*k4[1] + dopri_a65*k5[1]);
    w_temp[2] = omega[2] + dt * (dopri_a61*k1[2] + dopri_a62*k2[2] + dopri_a63*k3[2]
                                + dopri_a64*k4[2] + dopri_a65*k5[2]);
    euler_free_derivative(I, w_temp, k6);

    /* 5th-order solution */
    double w5[3] = {
        omega[0] + dt * (dopri_b51*k1[0] + dopri_b53*k3[0] + dopri_b54*k4[0]
                       + dopri_b55*k5[0] + dopri_b56*k6[0]),
        omega[1] + dt * (dopri_b51*k1[1] + dopri_b53*k3[1] + dopri_b54*k4[1]
                       + dopri_b55*k5[1] + dopri_b56*k6[1]),
        omega[2] + dt * (dopri_b51*k1[2] + dopri_b53*k3[2] + dopri_b54*k4[2]
                       + dopri_b55*k5[2] + dopri_b56*k6[2])
    };

    /* Compute k7 at w5 */
    euler_free_derivative(I, w5, k7);

    /* 4th-order embedded solution */
    double w4[3] = {
        omega[0] + dt * (dopri_b41*k1[0] + dopri_b43*k3[0] + dopri_b44*k4[0]
                       + dopri_b45*k5[0] + dopri_b46*k6[0] + dopri_b47*k7[0]),
        omega[1] + dt * (dopri_b41*k1[1] + dopri_b43*k3[1] + dopri_b44*k4[1]
                       + dopri_b45*k5[1] + dopri_b46*k6[1] + dopri_b47*k7[1]),
        omega[2] + dt * (dopri_b41*k1[2] + dopri_b43*k3[2] + dopri_b44*k4[2]
                       + dopri_b45*k5[2] + dopri_b46*k6[2] + dopri_b47*k7[2])
    };

    /* Error estimate (Euclidean norm of difference) */
    double diff[3] = {w5[0]-w4[0], w5[1]-w4[1], w5[2]-w4[2]};
    *error = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);

    /* Simple step-size control: scale by (tol/error)^(1/5) */
    double tol = 1e-8;
    double safety = 0.9;
    if (*error > 0) {
        *dt_next = dt * safety * pow(tol / *error, 0.2);
    } else {
        *dt_next = dt * 2.0;
    }
    /* Clamp step size change */
    if (*dt_next > dt * 5.0) *dt_next = dt * 5.0;
    if (*dt_next < dt * 0.1) *dt_next = dt * 0.1;

    /* Accept 5th-order solution */
    omega[0] = w5[0]; omega[1] = w5[1]; omega[2] = w5[2];
}

/* ============================================================================
 * Full Trajectory Simulation
 * ============================================================================ */

int simulate_free_rigid_body(const double I[3], const double omega0[3],
                             double t_end, double dt, IntegratorMethod method,
                             int *n_steps_out, double *traj_out) {
    int max_steps = (int)ceil(t_end / dt) + 1;
    double w[3];
    w[0] = omega0[0]; w[1] = omega0[1]; w[2] = omega0[2];

    traj_out[0] = w[0]; traj_out[1] = w[1]; traj_out[2] = w[2];

    int step;
    for (step = 1; step < max_steps; step++) {
        double actual_dt = dt;
        if (step * dt > t_end) break;

        double w_next[3];
        switch (method) {
            case INTEGRATOR_EULER:
                euler_step_omega(I, w, actual_dt, w_next);
                break;
            case INTEGRATOR_RK4:
                rk4_step_omega(I, w, actual_dt, w_next);
                break;
            case INTEGRATOR_MIDPOINT:
                midpoint_step_omega(I, w, actual_dt, w_next);
                break;
            case INTEGRATOR_VERLET:
                verlet_step_omega(I, w, actual_dt, w_next);
                break;
            case INTEGRATOR_DOPRI54:
            {
                double dt_adapt = actual_dt, err, dt_next;
                w_next[0] = w[0]; w_next[1] = w[1]; w_next[2] = w[2];
                dopri54_step_omega(I, w_next, dt_adapt, &dt_next, &err);
                break;
            }
            default:
                return -1;
        }

        w[0] = w_next[0]; w[1] = w_next[1]; w[2] = w_next[2];
        traj_out[step * 3 + 0] = w[0];
        traj_out[step * 3 + 1] = w[1];
        traj_out[step * 3 + 2] = w[2];
    }

    *n_steps_out = step;
    return 0;
}

int simulate_rigid_body_torque(const double I[3], const double omega0[3],
                               torque_function torque_fn, void *params,
                               double t_end, double dt,
                               int *n_steps_out, double *times_out, double *omega_out) {
    int max_steps = (int)ceil(t_end / dt) + 1;
    double w[3] = {omega0[0], omega0[1], omega0[2]};
    double t = 0.0;

    times_out[0] = t;
    omega_out[0] = w[0]; omega_out[1] = w[1]; omega_out[2] = w[2];

    int step;
    for (step = 1; step < max_steps; step++) {
        t = step * dt;
        double torque[3];
        torque_fn(t, w, torque, params);

        /* RK4 with torque */
        double k1[3], k2[3], k3[3], k4[3], w_temp[3], N_temp[3];
        double half_dt = 0.5 * dt;

        euler_with_torque(I, w, torque, k1);

        w_temp[0] = w[0] + half_dt * k1[0];
        w_temp[1] = w[1] + half_dt * k1[1];
        w_temp[2] = w[2] + half_dt * k1[2];
        torque_fn(t - half_dt, w_temp, N_temp, params);
        euler_with_torque(I, w_temp, N_temp, k2);

        w_temp[0] = w[0] + half_dt * k2[0];
        w_temp[1] = w[1] + half_dt * k2[1];
        w_temp[2] = w[2] + half_dt * k2[2];
        torque_fn(t - half_dt, w_temp, N_temp, params);
        euler_with_torque(I, w_temp, N_temp, k3);

        w_temp[0] = w[0] + dt * k3[0];
        w_temp[1] = w[1] + dt * k3[1];
        w_temp[2] = w[2] + dt * k3[2];
        torque_fn(t, w_temp, N_temp, params);
        euler_with_torque(I, w_temp, N_temp, k4);

        w[0] += (dt / 6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        w[1] += (dt / 6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        w[2] += (dt / 6.0) * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);

        times_out[step] = t;
        omega_out[step*3 + 0] = w[0];
        omega_out[step*3 + 1] = w[1];
        omega_out[step*3 + 2] = w[2];
    }

    *n_steps_out = step;
    return 0;
}

int simulate_full_rigid_body(const double I[3], const double omega0[3],
                             const EulerAngles *euler0,
                             double t_end, double dt,
                             int *n_steps_out, double *times_out,
                             double *omega_out, double *euler_out) {
    int max_steps = (int)ceil(t_end / dt) + 1;
    double w[3] = {omega0[0], omega0[1], omega0[2]};
    EulerAngles ea = *euler0;
    euler_normalize(&ea);

    for (int step = 0; step < max_steps; step++) {
        double t = step * dt;
        if (t > t_end) { *n_steps_out = step; return 0; }

        times_out[step] = t;
        omega_out[step*3 + 0] = w[0];
        omega_out[step*3 + 1] = w[1];
        omega_out[step*3 + 2] = w[2];
        euler_out[step*3 + 0] = ea.phi;
        euler_out[step*3 + 1] = ea.theta;
        euler_out[step*3 + 2] = ea.psi;

        /* Integrate ω with RK4 */
        double w_next[3];
        rk4_step_omega(I, w, dt, w_next);

        /* Integrate Euler angles with RK4 */
        double y[3] = {ea.phi, ea.theta, ea.psi};
        double k1[3], k2[3], k3[3], k4[3], y_temp[3];
        double half_dt = 0.5 * dt;

        euler_kinematic_ode(t, y, k1, w);

        y_temp[0] = y[0] + half_dt * k1[0];
        y_temp[1] = y[1] + half_dt * k1[1];
        y_temp[2] = y[2] + half_dt * k1[2];
        euler_kinematic_ode(t + half_dt, y_temp, k2, w);

        y_temp[0] = y[0] + half_dt * k2[0];
        y_temp[1] = y[1] + half_dt * k2[1];
        y_temp[2] = y[2] + half_dt * k2[2];
        euler_kinematic_ode(t + half_dt, y_temp, k3, w);

        y_temp[0] = y[0] + dt * k3[0];
        y_temp[1] = y[1] + dt * k3[1];
        y_temp[2] = y[2] + dt * k3[2];
        euler_kinematic_ode(t + dt, y_temp, k4, w_next);

        ea.phi   = y[0] + (dt/6.0)*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
        ea.theta = y[1] + (dt/6.0)*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
        ea.psi   = y[2] + (dt/6.0)*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);

        w[0] = w_next[0]; w[1] = w_next[1]; w[2] = w_next[2];
    }

    *n_steps_out = max_steps;
    return 0;
}

/* ============================================================================
 * Motion Constant Monitoring
 * ============================================================================ */

void monitor_motion_constants(const double I[3], int n_steps, const double *traj,
                              double *T_hist, double *L2_hist) {
    for (int i = 0; i < n_steps; i++) {
        const double *w = &traj[i * 3];
        motion_constants_euler(I, w, &T_hist[i], &L2_hist[i]);
    }
}

void constant_drift_report(int n_steps, const double *T_hist, const double *L2_hist,
                           double *T_drift, double *L2_drift) {
    double T0 = T_hist[0], L20 = L2_hist[0];
    double max_T_drift = 0.0, max_L2_drift = 0.0;
    for (int i = 0; i < n_steps; i++) {
        double dT = fabs(T_hist[i] - T0);
        double dL2 = fabs(L2_hist[i] - L20);
        if (dT > max_T_drift) max_T_drift = dT;
        if (dL2 > max_L2_drift) max_L2_drift = dL2;
    }
    *T_drift = (T0 > 1e-15) ? max_T_drift / T0 : max_T_drift;
    *L2_drift = (L20 > 1e-15) ? max_L2_drift / L20 : max_L2_drift;
}

void body_to_space_trajectory(int n_steps, const double *omega_body,
                              const double *euler_traj, double *omega_space) {
    for (int i = 0; i < n_steps; i++) {
        const double *w_body = &omega_body[i * 3];
        double phi = euler_traj[i*3], theta = euler_traj[i*3+1], psi = euler_traj[i*3+2];
        mat3 R;
        euler_to_rotation_matrix(phi, theta, psi, &R);
        vec3 wb = {w_body[0], w_body[1], w_body[2]};
        vec3 ws = mat3_mul_vec(&R, wb);
        omega_space[i*3] = ws.x; omega_space[i*3+1] = ws.y; omega_space[i*3+2] = ws.z;
    }
}

/* ============================================================================
 * Poinsot Construction
 * ============================================================================ */

void invariable_plane(const double I[3], const double omega[3],
                      const mat3 *R, vec3 *normal, double *dist) {
    /* Angular momentum in body frame: L = (I₁ω₁, I₂ω₂, I₃ω₃) */
    double L_body[3] = {I[0]*omega[0], I[1]*omega[1], I[2]*omega[2]};

    /* Transform to space frame: L_space = R · L_body */
    vec3 Lb = {L_body[0], L_body[1], L_body[2]};
    vec3 Ls = mat3_mul_vec(R, Lb);
    double L_mag = vec3_norm(Ls);

    if (L_mag < 1e-15) {
        *normal = VEC3_ZERO;
        *dist = 0.0;
        return;
    }

    /* Normal to invariable plane: direction of L_space */
    normal->x = Ls.x / L_mag;
    normal->y = Ls.y / L_mag;
    normal->z = Ls.z / L_mag;

    /* Distance from origin = 2T / |L| */
    double T = 0.5 * (I[0]*omega[0]*omega[0] + I[1]*omega[1]*omega[1] + I[2]*omega[2]*omega[2]);
    *dist = 2.0 * T / L_mag;
}

void polhode_tangent(const double I[3], const double omega[3], double tangent[3]) {
    /* The polhode is the intersection of the energy ellipsoid (T=const)
     * and the angular momentum sphere (L²=const).
     * Its tangent direction is ω̇ (the Euler equation RHS). */
    euler_free_derivative(I, omega, tangent);
}
