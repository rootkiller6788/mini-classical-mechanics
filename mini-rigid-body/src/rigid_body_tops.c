/**
 * rigid_body_tops.c — Heavy Symmetric Top, Lagrange Top, Gyroscope Theory
 *
 * Goldstein §5.7: Heavy symmetric top with one point fixed in gravity.
 * Landau §35-36: Motion of a top.
 *
 * The Lagrange top is a symmetric top (I₁ = I₂ ≠ I₃) with center of mass
 * at distance l from a fixed pivot, under uniform gravity.
 *
 * Constants of motion for Lagrange top:
 *   p_φ = ∂L/∂φ̇ = I₁ φ̇ sin²θ + I₃ (ψ̇ + φ̇ cosθ) cosθ = Lz  (space z-axis ang. mom.)
 *   p_ψ = ∂L/∂ψ̇ = I₃ (ψ̇ + φ̇ cosθ) = L3                  (body z-axis ang. mom.)
 *   E   = T + V (total energy)
 *
 * Reduced to 1D problem in θ:
 *   V_eff(θ) = (Lz - L₃ cosθ)²/(2 I₁ sin²θ) + M g l cosθ + L₃²/(2 I₃)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/rigid_body_types.h"
#include "../include/rigid_body_tops.h"

/* ============================================================================
 * Heavy Symmetric Top — Effective Potential
 * ============================================================================ */

double symmetric_top_eff_potential(double I1, double I3, double M, double g, double l,
                                   double Lz, double L3, double theta) {
    double st = sin(theta), ct = cos(theta);
    if (fabs(st) < 1e-15) {
        /* θ ≈ 0 or π, use limiting behavior */
        double V1;
        if (fabs(theta) < 1e-10) {
            /* θ→0: (Lz-L3)²/(2I1 θ²) + Mgl. Centrifugal barrier. */
            V1 = (fabs(Lz - L3) > 1e-10) ? INFINITY : 0.0;
        } else {
            /* θ→π: similar */
            V1 = (fabs(Lz + L3) > 1e-10) ? INFINITY : 0.0;
        }
        return V1 + M * g * l * ct + L3*L3/(2.0*I3);
    }

    double V1 = (Lz - L3*ct) * (Lz - L3*ct) / (2.0 * I1 * st * st);
    double V2 = M * g * l * ct;
    double V3 = L3 * L3 / (2.0 * I3);
    return V1 + V2 + V3;
}

void symmetric_top_eff_potential_grid(double I1, double I3, double M, double g, double l,
                                      double Lz, double L3, int n_pts,
                                      double *thetas, double *Veff) {
    for (int i = 0; i < n_pts; i++) {
        thetas[i] = 1e-4 + (M_PI - 2e-4) * i / (n_pts - 1);
        Veff[i] = symmetric_top_eff_potential(I1, I3, M, g, l, Lz, L3, thetas[i]);
    }
}

void nutation_turning_points(double I1, double I3, double M, double g, double l,
                             double Lz, double L3, double E, int n_pts,
                             int *n_tp_out, double tp_out[2]) {
    int n_tp = 0;
    double prev_val = 0.0;
    int prev_valid = 0;

    for (int i = 0; i < n_pts; i++) {
        double theta = 1e-4 + (M_PI - 2e-4) * i / (n_pts - 1);
        double Veff = symmetric_top_eff_potential(I1, I3, M, g, l, Lz, L3, theta);
        double diff = E - Veff;

        if (prev_valid && prev_val * diff < 0.0) {
            /* Sign change detected — root (turning point) between theta_{i-1} and theta_i */
            if (n_tp < 2) {
                /* Linear interpolation for root */
                double theta_prev = 1e-4 + (M_PI - 2e-4) * (i - 1) / (n_pts - 1);
                double t = -prev_val / (diff - prev_val);
                tp_out[n_tp] = theta_prev + t * (theta - theta_prev);
                n_tp++;
            }
        }
        /* Also detect near-zero crossings */
        if (fabs(diff) < 1e-12 && n_tp < 2) {
            int is_dup = 0;
            for (int j = 0; j < n_tp; j++) {
                if (fabs(theta - tp_out[j]) < 1e-6) { is_dup = 1; break; }
            }
            if (!is_dup) tp_out[n_tp++] = theta;
        }

        prev_val = diff;
        prev_valid = 1;
    }

    *n_tp_out = n_tp;
}

/* ============================================================================
 * Heavy Symmetric Top — Angular Rates
 * ============================================================================ */

double nutation_rate(double I1, double I3, double M, double g, double l,
                     double E, double Lz, double L3, double theta) {
    double Veff = symmetric_top_eff_potential(I1, I3, M, g, l, Lz, L3, theta);
    double arg = (E - Veff) * 2.0 / I1;
    if (arg < 0.0) return 0.0; /* Classically forbidden region */
    return sqrt(arg);
}

double precession_rate_top(double I1, double Lz, double L3, double theta) {
    double st = sin(theta);
    if (fabs(st) < 1e-15) return 0.0;
    return (Lz - L3 * cos(theta)) / (I1 * st * st);
}

double spin_rate_top(double I1, double I3, double Lz, double L3, double theta) {
    double phidot = precession_rate_top(I1, Lz, L3, theta);
    return L3 / I3 - phidot * cos(theta);
}

/* ============================================================================
 * Steady Precession
 * ============================================================================ */

static double dVeff_dtheta(double I1, double M, double g, double l,
                           double Lz, double L3, double theta) {
    double st = sin(theta), ct = cos(theta);
    if (fabs(st) < 1e-15) return 0.0;
    /* d/dθ [ (Lz - L3 ct)²/(2 I1 st²) + Mgl ct ]
     * = -(Lz - L3 ct)(Lz ct - L3)/(I1 st³) - Mgl st */
    double term1 = -(Lz - L3*ct) * (Lz*ct - L3) / (I1 * st*st*st);
    double term2 = -M * g * l * st;
    return term1 + term2;
}

int steady_precession_angle(double I1, double I3, double M, double g, double l,
                            double Lz, double L3, const double bracket[2],
                            double *theta0) {
    (void)I3; /* Not needed for dVeff/dθ = 0, but kept for API consistency */
    double a = bracket[0], b = bracket[1];
    double fa = dVeff_dtheta(I1, M, g, l, Lz, L3, a);
    double fb = dVeff_dtheta(I1, M, g, l, Lz, L3, b);

    if (fa * fb > 0.0) return -1; /* No sign change */

    for (int iter = 0; iter < 60; iter++) {
        double mid = 0.5 * (a + b);
        double fmid = dVeff_dtheta(I1, M, g, l, Lz, L3, mid);
        if (fabs(fmid) < 1e-12 || (b - a) < 1e-14) {
            *theta0 = mid;
            return 0;
        }
        if (fa * fmid < 0.0) {
            b = mid; fb = fmid;
        } else {
            a = mid; fa = fmid;
        }
    }
    *theta0 = 0.5 * (a + b);
    return 0;
}

double steady_precession_rate(double I1, double Lz, double L3, double theta0) {
    return precession_rate_top(I1, Lz, L3, theta0);
}

double fast_top_precession_rate(double M, double g, double l, double I3, double omega3) {
    /* φ̇ ≈ Mgl / (I₃ ω₃) — gyroscopic approximation */
    if (fabs(omega3) < 1e-15) return 0.0;
    return M * g * l / (I3 * omega3);
}

double slow_precession_rate(double M, double g, double l, double L3) {
    if (fabs(L3) < 1e-15) return 0.0;
    return M * g * l / L3;
}

/* ============================================================================
 * Sleeping Top
 * ============================================================================ */

int sleeping_top_is_stable(double M, double g, double l, double I1, double I3,
                           double omega3, double *omega_crit_out) {
    double omega_crit = sqrt(4.0 * M * g * l * I1) / I3;
    if (omega_crit_out) *omega_crit_out = omega_crit;
    return (omega3 > omega_crit) ? 1 : 0;
}

double sleeping_top_critical_spin(double M, double g, double l, double I1, double I3) {
    return sqrt(4.0 * M * g * l * I1) / I3;
}

/* ============================================================================
 * Lagrange Top — Full ODE System (6D state space)
 *
 * State: [θ, φ, ψ, p_θ, p_φ, p_ψ]
 *
 * Equations of motion:
 *   θ̇ = p_θ / I1
 *   φ̇ = (p_φ - p_ψ cosθ) / (I1 sin²θ)
 *   ψ̇ = p_ψ/I3 - φ̇ cosθ
 *   ṗ_θ = (p_φ - p_ψ cosθ)(p_φ cosθ - p_ψ) / (I1 sin³θ) - Mgl sinθ
 *   ṗ_φ = 0   (Lz conserved)
 *   ṗ_ψ = 0   (L3 conserved)
 * ============================================================================ */

void lagrange_top_ode(double t, const double y[], double dydt[], const void *params) {
    (void)t;
    const LagrangeTopParams *p = (const LagrangeTopParams *)params;
    double I1 = p->I1, I3 = p->I3, M = p->M, g = p->g, l = p->l;

    double theta = y[0]; /* phi = y[1], psi = y[2] */
    double ptheta = y[3];
    double Lz = y[4];    /* p_φ — conserved */
    double L3 = y[5];    /* p_ψ — conserved */

    double st = sin(theta), ct = cos(theta);

    /* θ̇ = p_θ / I1 */
    dydt[0] = ptheta / I1;

    /* φ̇ = (Lz - L3 cosθ) / (I1 sin²θ) */
    if (fabs(st) > 1e-15) {
        dydt[1] = (Lz - L3 * ct) / (I1 * st * st);
    } else {
        dydt[1] = 0.0;
    }
    double phidot = dydt[1];

    /* ψ̇ = L3/I3 - φ̇ cosθ */
    dydt[2] = L3 / I3 - phidot * ct;

    /* ṗ_θ = (Lz - L3 cosθ)(Lz cosθ - L3)/(I1 sin³θ) - Mgl sinθ */
    if (fabs(st) > 1e-15) {
        dydt[3] = (Lz - L3*ct) * (Lz*ct - L3) / (I1 * st*st*st) - M*g*l*st;
    } else {
        dydt[3] = 0.0;
    }

    /* Constants of motion */
    dydt[4] = 0.0; /* Lz conserved */
    dydt[5] = 0.0; /* L3 conserved */
}

int simulate_lagrange_top(const LagrangeTopParams *params,
                          double theta0, double phi0, double psi0,
                          double theta_dot0, double Lz, double L3,
                          double t_end, double dt,
                          int *n_steps_out, double *times_out,
                          double *thetas_out, double *phis_out, double *psis_out) {
    int max_steps = (int)ceil(t_end / dt) + 1;
    double ptheta0 = params->I1 * theta_dot0;
    double y[6] = {theta0, phi0, psi0, ptheta0, Lz, L3};

    for (int step = 0; step < max_steps; step++) {
        double t = step * dt;
        if (t > t_end) { *n_steps_out = step; return 0; }

        times_out[step] = t;
        thetas_out[step] = y[0];
        phis_out[step] = y[1];
        psis_out[step] = y[2];

        /* RK4 integration */
        double k1[6], k2[6], k3[6], k4[6], y_temp[6];
        double half_dt = 0.5 * dt;

        lagrange_top_ode(t, y, k1, params);

        for (int j = 0; j < 6; j++) y_temp[j] = y[j] + half_dt * k1[j];
        lagrange_top_ode(t + half_dt, y_temp, k2, params);

        for (int j = 0; j < 6; j++) y_temp[j] = y[j] + half_dt * k2[j];
        lagrange_top_ode(t + half_dt, y_temp, k3, params);

        for (int j = 0; j < 6; j++) y_temp[j] = y[j] + dt * k3[j];
        lagrange_top_ode(t + dt, y_temp, k4, params);

        for (int j = 0; j < 6; j++) {
            y[j] += (dt / 6.0) * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
        }
    }

    *n_steps_out = max_steps;
    return 0;
}

/* ============================================================================
 * Nutation Analysis
 * ============================================================================ */

double small_nutation_period(double I1, double M, double g, double l) {
    /* T_nut ≈ 2π √(I₁/(Mgl)) — small-amplitude nutation approximation */
    if (M * g * l < 1e-15) return INFINITY;
    return M_2PI * sqrt(I1 / (M * g * l));
}

PrecessionType classify_precession(double I1, double Lz, double L3,
                                   double theta_min, double theta_max) {
    /* Check φ̇ at both turning points and interior */
    if (theta_max - theta_min < 1e-10) return PRECESSION_MONOTONIC;

    double phi_min = precession_rate_top(I1, Lz, L3, theta_min);
    double phi_max = precession_rate_top(I1, Lz, L3, theta_max);
    double phi_mid = precession_rate_top(I1, Lz, L3, 0.5*(theta_min + theta_max));

    /* Check if φ̇ changes sign */
    if (phi_min * phi_mid < 0.0 || phi_mid * phi_max < 0.0 || phi_min * phi_max < 0.0) {
        return PRECESSION_LOOPING;
    }

    /* Check if φ̇ ≈ 0 at turning points (cusp) */
    if (fabs(phi_min) < 1e-12 || fabs(phi_max) < 1e-12) {
        return PRECESSION_CUSPED;
    }

    return PRECESSION_MONOTONIC;
}

/* ============================================================================
 * Gyroscope Principles
 * ============================================================================ */

void gyroscopic_torque(const double L_spin[3], const double omega_forced[3], double torque[3]) {
    /* N = -ω_forced × L_spin (the torque the gyro exerts on its support)
     * or equivalently, N_needed = ω_forced × L_spin (torque needed to precess) */
    vec3 L = {L_spin[0], L_spin[1], L_spin[2]};
    vec3 w = {omega_forced[0], omega_forced[1], omega_forced[2]};
    vec3 N = vec3_cross(L, w); /* Gyroscopic torque N = ω × L convention */
    torque[0] = N.x; torque[1] = N.y; torque[2] = N.z;
}

double gyroscope_precession_rate(double N, double L_spin) {
    if (fabs(L_spin) < 1e-15) return 0.0;
    return N / L_spin;
}

void gyroscope_2axis_simulate(double L_spin, const double omega_base[3],
                              double I_gimbal, double K_spring, double D_damping,
                              double dt, int n_steps,
                              double *theta_out, double *omega_out) {
    /* 2nd-order ODE for gimbal angle θ:
     *   I_gimbal θ̈ + D θ̇ + K θ = L_spin × ω_base (projected onto gimbal axis)
     * Simplified: using the component of ω_base perpendicular to spin axis. */
    double theta = 0.0, dtheta = 0.0;

    /* Driving torque from base rotation perpendicular to spin axis */
    double N_drive = L_spin * omega_base[1]; /* Assume spin along z, base ω along y */

    for (int i = 0; i < n_steps; i++) {
        /* Euler integration (for demonstration; RK4 would be better) */
        double ddtheta = (N_drive - K_spring * theta - D_damping * dtheta) / I_gimbal;
        dtheta += ddtheta * dt;
        theta += dtheta * dt;

        theta_out[i] = theta;
        omega_out[i] = dtheta;
    }
}

double free_gyro_nutation_frequency(double L_spin, double I_trans) {
    if (I_trans < 1e-15) return 0.0;
    return L_spin / I_trans;
}
