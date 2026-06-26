/**
 * @file hamiltonian_symplectic.c
 * @brief Symplectic integrator implementations for Hamiltonian systems.
 *
 * Implements: Stormer-Verlet, Symplectic Euler, Ruth-4, Yoshida-6, RK4,
 * Hamiltonian flow solver, energy diagnostics.
 *
 * Reference: Hairer, Lubich, Wanner — Geometric Numerical Integration (2006)
 *            Forest & Ruth, Physica D 43 (1990)
 *            Yoshida, Phys. Lett. A 150 (1990)
 *
 * Knowledge coverage:
 *   L5: Symplectic integration methods (4 explicit + 1 reference method)
 *
 * University mapping:
 *   ETH 402-0800 — Geometric integration theory and practice
 *   MIT 8.012 — Numerical integration of Hamilton's equations
 *   Princeton PHY 505 — Long-time integration of dynamical systems
 */

#include "hamiltonian_types.h"
#include "hamiltonian_equations.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ──────────────────────────────────────────────────────────────
 * Stormer-Verlet (2nd order, separable Hamiltonian)
 *
 * For H = T(p) + V(q):
 *   p_{n+1/2} = p_n - (h/2)·grad_V(q_n)
 *   q_{n+1}   = q_n +  h ·grad_T(p_{n+1/2})
 *   p_{n+1}   = p_{n+1/2} - (h/2)·grad_V(q_{n+1})
 *
 * Properties:
 *   — Symplectic: preserves canonical 2-form exactly for quadratic H
 *   — Time-reversible: psi_{-h} ∘ psi_h = id
 *   — 2nd order accurate
 *
 * The method is the "gold standard" for molecular dynamics and
 * celestial mechanics due to its excellent long-term stability.
 * ────────────────────────────────────────────────────────────── */

void stormer_verlet_step(const HamiltonianSystem *sys,
                          const double *q, const double *p,
                          double dt,
                          double *q_new, double *p_new) {
    if (!sys || !q || !p || !q_new || !p_new) return;
    size_t n = sys->n;
    double *dH_dq = (double *)malloc(n * sizeof(double));
    double *dH_dp = (double *)malloc(n * sizeof(double));
    double *p_half = (double *)malloc(n * sizeof(double));
    if (!dH_dq || !dH_dp || !p_half) {
        free(dH_dq); free(dH_dp); free(p_half); return;
    }

    /* Half-step in p using q_n */
    sys->grad_H(q, p, n, dH_dq, dH_dp);
    for (size_t i = 0; i < n; i++) {
        p_half[i] = p[i] - 0.5 * dt * dH_dq[i];
    }

    /* Full step in q using p_{n+1/2} */
    sys->grad_H(q, p_half, n, dH_dq, dH_dp);
    for (size_t i = 0; i < n; i++) {
        q_new[i] = q[i] + dt * dH_dp[i];
    }

    /* Complete p-step using q_{n+1} */
    sys->grad_H(q_new, p_half, n, dH_dq, dH_dp);
    for (size_t i = 0; i < n; i++) {
        p_new[i] = p_half[i] - 0.5 * dt * dH_dq[i];
    }

    free(dH_dq); free(dH_dp); free(p_half);
}

/* ──────────────────────────────────────────────────────────────
 * Symplectic Euler (1st order, q-first variant)
 *
 *   q_{n+1} = q_n + h·∂H/∂p(q_n, p_n)
 *   p_{n+1} = p_n - h·∂H/∂q(q_{n+1}, p_n)
 *
 * This is the simplest symplectic method. While only first-order
 * accurate, its symplectic property ensures bounded energy error
 * for integrable systems (unlike non-symplectic Euler).
 *
 * When composed with its adjoint (p-first variant), one recovers
 * the 2nd-order Stormer-Verlet via Trotter-Suzuki:
 *   psi_h^{Verlet} = phi_{h/2}^q ∘ phi_{h/2}^p
 * ────────────────────────────────────────────────────────────── */

void symplectic_euler_step(const HamiltonianSystem *sys,
                            const double *q, const double *p,
                            double dt,
                            double *q_new, double *p_new) {
    if (!sys || !q || !p || !q_new || !p_new) return;
    size_t n = sys->n;
    double *dH_dq = (double *)malloc(n * sizeof(double));
    double *dH_dp = (double *)malloc(n * sizeof(double));
    if (!dH_dq || !dH_dp) { free(dH_dq); free(dH_dp); return; }

    /* q-step using (q_n, p_n) */
    sys->grad_H(q, p, n, dH_dq, dH_dp);
    for (size_t i = 0; i < n; i++) {
        q_new[i] = q[i] + dt * dH_dp[i];
    }

    /* p-step using (q_{n+1}, p_n) */
    sys->grad_H(q_new, p, n, dH_dq, dH_dp);
    for (size_t i = 0; i < n; i++) {
        p_new[i] = p[i] - dt * dH_dq[i];
    }

    free(dH_dq); free(dH_dp);
}

/* ──────────────────────────────────────────────────────────────
 * Ruth 4th-order symmetric composition integrator
 *
 * Forest-Ruth (1990): minimum-stage 4th-order symplectic method.
 * Uses 4 stages with coefficients involving the cube root of 2.
 *
 * Coefficients for the symmetric composition:
 *   c1 = c4 = 1/(2(2-2^{1/3})),  c2 = c3 = (1-2^{1/3})/(2(2-2^{1/3}))
 *   d1 = d3 = 1/(2-2^{1/3}),      d2 = -2^{1/3}/(2-2^{1/3}),  d4 = 0
 *
 * Note: d2 is negative, which is characteristic of higher-order
 * symplectic composition methods (the backward step enables
 * cancellation of lower-order errors).
 *
 * The method requires 3 force evaluations per step despite being
 * 4th order (vs. 4 for standard RK4).
 * ────────────────────────────────────────────────────────────── */

void ruth4_step(const HamiltonianSystem *sys,
                 const double *q, const double *p,
                 double dt,
                 double *q_new, double *p_new) {
    if (!sys || !q || !p || !q_new || !p_new) return;
    size_t n = sys->n;

    const double cbrt2 = 1.2599210498948732;  /* 2^{1/3} */
    const double w = 1.0 / (2.0 - cbrt2);
    const double c[4] = { w/2.0, (1.0-cbrt2)*w/2.0, (1.0-cbrt2)*w/2.0, w/2.0 };
    const double d[4] = { w, -cbrt2*w, w, 0.0 };

    double *q_tmp = (double *)malloc(n * sizeof(double));
    double *p_tmp = (double *)malloc(n * sizeof(double));
    double *dH_dq = (double *)malloc(n * sizeof(double));
    double *dH_dp = (double *)malloc(n * sizeof(double));
    if (!q_tmp || !p_tmp || !dH_dq || !dH_dp) {
        free(q_tmp); free(p_tmp); free(dH_dq); free(dH_dp); return;
    }

    memcpy(q_tmp, q, n * sizeof(double));
    memcpy(p_tmp, p, n * sizeof(double));

    for (int stage = 0; stage < 4; stage++) {
        /* p-drift stage */
        if (fabs(d[stage]) > 1e-15) {
            sys->grad_H(q_tmp, p_tmp, n, dH_dq, dH_dp);
            for (size_t i = 0; i < n; i++) {
                p_tmp[i] += dt * d[stage] * (-dH_dq[i]);
            }
        }
        /* q-kick stage */
        if (fabs(c[stage]) > 1e-15) {
            sys->grad_H(q_tmp, p_tmp, n, dH_dq, dH_dp);
            for (size_t i = 0; i < n; i++) {
                q_tmp[i] += dt * c[stage] * dH_dp[i];
            }
        }
    }

    memcpy(q_new, q_tmp, n * sizeof(double));
    memcpy(p_new, p_tmp, n * sizeof(double));

    free(q_tmp); free(p_tmp); free(dH_dq); free(dH_dp);
}

/* ──────────────────────────────────────────────────────────────
 * Yoshida 6th-order symmetric composition integrator
 *
 * Constructed by composing 3 Stormer-Verlet (leapfrog) steps with
 * irrationally-related step sizes. The coefficients w_i satisfy
 * sum w_i = 1 and produce a 6th-order method.
 *
 * w = [w3, w2, w1, w0, w0, w1, w2, w3] (symmetric)
 *
 * The negative coefficients (w1 < 0) are unavoidable for explicit
 * symplectic methods beyond 2nd order (Sheng-Suzuki theorem).
 *
 * Reference: Yoshida, Phys. Lett. A 150, 262-268 (1990)
 * ────────────────────────────────────────────────────────────── */

void yoshida6_step(const HamiltonianSystem *sys,
                    const double *q, const double *p,
                    double dt,
                    double *q_new, double *p_new) {
    if (!sys || !q || !p || !q_new || !p_new) return;
    size_t n = sys->n;

    /* Yoshida (1990) 6th-order coefficients */
    const double w[8] = {
         0.392256805238780,   /* w3 */
         0.510043411918458,   /* w2 */
        -0.471053385409757,   /* w1 (negative — required beyond 2nd order) */
         0.068753168252520,   /* w0 */
         0.068753168252520,   /* w0 */
        -0.471053385409757,   /* w1 */
         0.510043411918458,   /* w2 */
         0.392256805238780    /* w3 */
    };

    double *q_tmp = (double *)malloc(n * sizeof(double));
    double *p_tmp = (double *)malloc(n * sizeof(double));
    if (!q_tmp || !p_tmp) { free(q_tmp); free(p_tmp); return; }

    memcpy(q_tmp, q, n * sizeof(double));
    memcpy(p_tmp, p, n * sizeof(double));

    /* Compose 8 Verlet substeps */
    for (int s = 0; s < 8; s++) {
        double h = w[s] * dt;
        stormer_verlet_step(sys, q_tmp, p_tmp, h, q_tmp, p_tmp);
    }

    memcpy(q_new, q_tmp, n * sizeof(double));
    memcpy(p_new, p_tmp, n * sizeof(double));

    free(q_tmp); free(p_tmp);
}

/* ──────────────────────────────────────────────────────────────
 * Classical 4th-order Runge-Kutta (non-symplectic reference)
 *
 * Standard RK4: 4 stages, 4 function evaluations, 4th-order accurate.
 * Unlike symplectic methods, RK4 does NOT conserve symplectic structure,
 * leading to systematic energy drift in long-time integrations.
 *
 * The method is included for educational comparison:
 *   — RK4: small short-time error, but energy drift over many periods
 *   — Verlet: larger short-time error, but no long-term energy drift
 *
 * This illustrates the key insight of geometric integration:
 * preserving geometric structure can be more important than
 * maximizing formal order of accuracy.
 * ────────────────────────────────────────────────────────────── */

void rk4_hamiltonian_step(const HamiltonianSystem *sys,
                           const double *q, const double *p,
                           double dt,
                           double *q_new, double *p_new) {
    if (!sys || !q || !p || !q_new || !p_new) return;
    size_t n = sys->n;
    size_t n2 = 2 * n;

    double *y = (double *)malloc(n2 * sizeof(double));
    double *k1 = (double *)malloc(n2 * sizeof(double));
    double *k2 = (double *)malloc(n2 * sizeof(double));
    double *k3 = (double *)malloc(n2 * sizeof(double));
    double *k4 = (double *)malloc(n2 * sizeof(double));
    double *ytmp = (double *)malloc(n2 * sizeof(double));
    if (!y || !k1 || !k2 || !k3 || !k4 || !ytmp) {
        free(y); free(k1); free(k2); free(k3); free(k4); free(ytmp); return;
    }

    /* y = [q; p] */
    memcpy(y, q, n * sizeof(double));
    memcpy(y + n, p, n * sizeof(double));

    /* Stage 1 */
    hamiltons_rhs(sys, y, y + n, k1, k1 + n);
    /* Stage 2 */
    for (size_t i = 0; i < n2; i++) ytmp[i] = y[i] + 0.5 * dt * k1[i];
    hamiltons_rhs(sys, ytmp, ytmp + n, k2, k2 + n);
    /* Stage 3 */
    for (size_t i = 0; i < n2; i++) ytmp[i] = y[i] + 0.5 * dt * k2[i];
    hamiltons_rhs(sys, ytmp, ytmp + n, k3, k3 + n);
    /* Stage 4 */
    for (size_t i = 0; i < n2; i++) ytmp[i] = y[i] + dt * k3[i];
    hamiltons_rhs(sys, ytmp, ytmp + n, k4, k4 + n);

    /* Combine */
    for (size_t i = 0; i < n2; i++) {
        y[i] += (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

    memcpy(q_new, y, n * sizeof(double));
    memcpy(p_new, y + n, n * sizeof(double));

    free(y); free(k1); free(k2); free(k3); free(k4); free(ytmp);
}

/* ──────────────────────────────────────────────────────────────
 * Hamiltonian flow solver with fixed-step integration
 *
 * Integrates from t=0 to t=t_end with constant step size dt.
 * Records trajectory points at regular intervals.
 *
 * The choice of integrator profoundly affects qualitative behavior:
 *   — Symplectic methods: bounded energy error, correct phase portrait
 *   — Non-symplectic methods: systematic energy drift, distorted orbits
 *
 * For long-time integrations (many orbital periods), symplectic
 * methods are strongly preferred despite their sometimes lower
 * formal order of accuracy.
 * ────────────────────────────────────────────────────────────── */

size_t solve_hamiltonian(const HamiltonianSystem *sys,
                          const double *q0, const double *p0,
                          double t_end, double dt,
                          IntegratorType integrator,
                          PhaseTrajectory *traj,
                          size_t record_every) {
    if (!sys || !q0 || !p0 || !traj || dt <= 0.0) return 0;
    size_t n = sys->n;

    double *q = (double *)malloc(n * sizeof(double));
    double *p = (double *)malloc(n * sizeof(double));
    double *q_new = (double *)malloc(n * sizeof(double));
    double *p_new = (double *)malloc(n * sizeof(double));
    if (!q || !p || !q_new || !p_new) {
        free(q); free(p); free(q_new); free(p_new); return 0;
    }

    memcpy(q, q0, n * sizeof(double));
    memcpy(p, p0, n * sizeof(double));

    double t = 0.0;
    size_t step = 0;
    double H_val = sys->H(q, p, n);
    phase_trajectory_append(traj, t, q, p, H_val);

    while (t < t_end) {
        double h = (dt < t_end - t) ? dt : (t_end - t);
        switch (integrator) {
            case INTEGRATOR_STORMER_VERLET:
                stormer_verlet_step(sys, q, p, h, q_new, p_new); break;
            case INTEGRATOR_SYMPLECTIC_EULER:
                symplectic_euler_step(sys, q, p, h, q_new, p_new); break;
            case INTEGRATOR_RUTH_4TH:
                ruth4_step(sys, q, p, h, q_new, p_new); break;
            case INTEGRATOR_YOSHIDA_6TH:
                yoshida6_step(sys, q, p, h, q_new, p_new); break;
            case INTEGRATOR_RK4:
                rk4_hamiltonian_step(sys, q, p, h, q_new, p_new); break;
            default:
                free(q); free(p); free(q_new); free(p_new); return step;
        }
        memcpy(q, q_new, n * sizeof(double));
        memcpy(p, p_new, n * sizeof(double));
        t += h;
        step++;

        if (step % record_every == 0) {
            H_val = sys->H(q, p, n);
            phase_trajectory_append(traj, t, q, p, H_val);
        }
    }

    /* Record final point */
    if (traj->length == 0 || traj->times[traj->length - 1] < t_end - 1e-12) {
        H_val = sys->H(q, p, n);
        phase_trajectory_append(traj, t, q, p, H_val);
    }

    free(q); free(p); free(q_new); free(p_new);
    return step;
}

/* ──────────────────────────────────────────────────────────────
 * Integrator descriptor factory
 *
 * Creates metadata for each integrator type:
 *   order, stage count, symplectic property, time-reversibility.
 * ────────────────────────────────────────────────────────────── */

IntegratorDescriptor integrator_descriptor_create(IntegratorType type) {
    IntegratorDescriptor desc;
    memset(&desc, 0, sizeof(desc));
    desc.type = type;
    desc.coeffs = NULL;
    switch (type) {
        case INTEGRATOR_SYMPLECTIC_EULER:
            desc.order = 1; desc.n_stages = 1;
            desc.is_symplectic = true; desc.is_time_reversible = false;
            break;
        case INTEGRATOR_STORMER_VERLET:
            desc.order = 2; desc.n_stages = 2;
            desc.is_symplectic = true; desc.is_time_reversible = true;
            break;
        case INTEGRATOR_RUTH_4TH:
            desc.order = 4; desc.n_stages = 4;
            desc.is_symplectic = true; desc.is_time_reversible = true;
            break;
        case INTEGRATOR_YOSHIDA_6TH:
            desc.order = 6; desc.n_stages = 8;
            desc.is_symplectic = true; desc.is_time_reversible = true;
            break;
        case INTEGRATOR_RK4:
            desc.order = 4; desc.n_stages = 4;
            desc.is_symplectic = false; desc.is_time_reversible = false;
            break;
    }
    return desc;
}

void integrator_descriptor_free(IntegratorDescriptor *desc) {
    if (desc) { free(desc->coeffs); desc->coeffs = NULL; }
}

/* ──────────────────────────────────────────────────────────────
 * Energy diagnostics
 *
 * These functions monitor the quality of numerical integration
 * by tracking conservation laws.
 * ────────────────────────────────────────────────────────────── */

double energy_drift(const PhaseTrajectory *traj) {
    if (!traj || traj->length < 2) return 0.0;
    double H0 = traj->H_values[0];
    if (fabs(H0) < 1e-15) return 0.0;
    double max_drift = 0.0;
    for (size_t i = 1; i < traj->length; i++) {
        double drift = fabs(traj->H_values[i] - H0) / fabs(H0);
        if (drift > max_drift) max_drift = drift;
    }
    return max_drift;
}

double action_drift(const PhaseTrajectory *traj,
                     const ActionAngleSystem *aas) {
    if (!traj || !aas || traj->length < 2) return 0.0;
    size_t n = aas->n;
    double *J0 = (double *)malloc(n * sizeof(double));
    double *Jf = (double *)malloc(n * sizeof(double));
    double *theta = (double *)malloc(n * sizeof(double));
    if (!J0 || !Jf || !theta) {
        free(J0); free(Jf); free(theta); return DBL_MAX;
    }
    aas->to_action_angle(traj->q_array[0], traj->p_array[0], n, J0, theta);
    aas->to_action_angle(traj->q_array[traj->length - 1],
                          traj->p_array[traj->length - 1], n, Jf, theta);
    double max_drift = 0.0;
    for (size_t i = 0; i < n; i++) {
        double drift = fabs(Jf[i] - J0[i]);
        if (drift > max_drift) max_drift = drift;
    }
    free(J0); free(Jf); free(theta);
    return max_drift;
}

double energy_error(const HamiltonianSystem *sys,
                     const double *q, const double *p,
                     double E_expected) {
    if (!sys || !q || !p) return DBL_MAX;
    double H_val = sys->H(q, p, sys->n);
    if (fabs(E_expected) < 1e-15) return fabs(H_val - E_expected);
    return fabs(H_val - E_expected) / fabs(E_expected);
}
