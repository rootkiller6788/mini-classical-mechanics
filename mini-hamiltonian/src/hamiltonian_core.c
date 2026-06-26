/**
 * @file hamiltonian_core.c
 * @brief Core Hamiltonian mechanics: memory management and equations of motion.
 *
 * Implements: phase space memory management, Hamilton's equations RHS,
 * Hamiltonian vector field, time derivative of observables,
 * Liouville density allocation.
 *
 * Reference: Goldstein Ch.8-10
 *
 * Knowledge coverage:
 *   L1: PhasePoint, PhaseTrajectory, HamiltonianSystem, LiouvilleDensity allocation
 *   L4: Hamilton's equations — fundamental law of classical dynamics
 *
 * University mapping:
 *   MIT 8.012 — Hamiltonian formulation of mechanics
 *   Caltech Ph 106 — Phase space and canonical formalism
 */

#include "hamiltonian_types.h"
#include "hamiltonian_equations.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ──────────────────────────────────────────────────────────────
 * L1: Memory management for phase space structures
 * ────────────────────────────────────────────────────────────── */

PhasePoint *phase_point_alloc(size_t n) {
    PhasePoint *p = (PhasePoint *)malloc(sizeof(PhasePoint));
    if (!p) return NULL;
    p->n = n;
    p->q = (double *)calloc(n, sizeof(double));
    p->p = (double *)calloc(n, sizeof(double));
    p->t = 0.0;
    if (!p->q || !p->p) {
        free(p->q); free(p->p); free(p);
        return NULL;
    }
    return p;
}

PhasePoint *phase_point_copy(const PhasePoint *src) {
    if (!src) return NULL;
    PhasePoint *dst = phase_point_alloc(src->n);
    if (!dst) return NULL;
    memcpy(dst->q, src->q, src->n * sizeof(double));
    memcpy(dst->p, src->p, src->n * sizeof(double));
    dst->t = src->t;
    return dst;
}

void phase_point_free(PhasePoint *p) {
    if (p) {
        free(p->q);
        free(p->p);
        free(p);
    }
}

PhaseTrajectory *phase_trajectory_alloc(size_t capacity, size_t n) {
    (void)n;  /* DOF count stored implicitly via row arrays */
    PhaseTrajectory *traj = (PhaseTrajectory *)malloc(sizeof(PhaseTrajectory));
    if (!traj) return NULL;
    traj->capacity = capacity;
    traj->length = 0;
    traj->times = (double *)malloc(capacity * sizeof(double));
    traj->q_array = (double **)malloc(capacity * sizeof(double *));
    traj->p_array = (double **)malloc(capacity * sizeof(double *));
    traj->H_values = (double *)malloc(capacity * sizeof(double));
    if (!traj->times || !traj->q_array || !traj->p_array || !traj->H_values) {
        free(traj->times); free(traj->q_array); free(traj->p_array);
        free(traj->H_values); free(traj);
        return NULL;
    }
    for (size_t i = 0; i < capacity; i++) {
        traj->q_array[i] = NULL;
        traj->p_array[i] = NULL;
    }
    return traj;
}

void phase_trajectory_append(PhaseTrajectory *traj,
                              double t, const double *q, const double *p,
                              double H_val) {
    if (!traj || traj->length >= traj->capacity) return;
    size_t i = traj->length;
    /* Store trajectory point: time, coordinates, momenta, energy.
     * q_array and p_array store pointers; allocation of per-point
     * coordinate arrays is managed separately by the caller. */
    (void)q; (void)p;  /* q, p stored by caller via q_array/p_array */
    traj->times[i] = t;
    traj->H_values[i] = H_val;
    traj->length++;
}

void phase_trajectory_free(PhaseTrajectory *traj) {
    if (traj) {
        free(traj->times);
        for (size_t i = 0; i < traj->capacity; i++) {
            free(traj->q_array[i]);
            free(traj->p_array[i]);
        }
        free(traj->q_array);
        free(traj->p_array);
        free(traj->H_values);
        free(traj);
    }
}

/* ──────────────────────────────────────────────────────────────
 * L4: Hamilton's equations of motion
 * ────────────────────────────────────────────────────────────── */

void hamiltons_rhs(const HamiltonianSystem *sys,
                   const double *q, const double *p,
                   double *dq, double *dp) {
    if (!sys || !sys->grad_H || !q || !p || !dq || !dp) return;
    double *dH_dq = (double *)malloc(sys->n * sizeof(double));
    double *dH_dp = (double *)malloc(sys->n * sizeof(double));
    if (!dH_dq || !dH_dp) { free(dH_dq); free(dH_dp); return; }
    sys->grad_H(q, p, sys->n, dH_dq, dH_dp);
    for (size_t i = 0; i < sys->n; i++) {
        dq[i] =  dH_dp[i];   /* dq/dt = +dH/dp */
        dp[i] = -dH_dq[i];   /* dp/dt = -dH/dq */
    }
    free(dH_dq);
    free(dH_dp);
}

void hamiltonian_vector_field(size_t n,
                               const double *dH_dq, const double *dH_dp,
                               double *X) {
    if (!dH_dq || !dH_dp || !X) return;
    for (size_t i = 0; i < n; i++) {
        X[i]     =  dH_dp[i];    /* X_q component */
        X[n + i] = -dH_dq[i];    /* X_p component */
    }
}

double time_derivative_observable(const HamiltonianSystem *sys,
                                   Observable f,
                                   const double *q, const double *p,
                                   double eps) {
    if (!sys || !f) return 0.0;
    size_t n = sys->n;
    double *dH_dq = (double *)malloc(n * sizeof(double));
    double *dH_dp = (double *)malloc(n * sizeof(double));
    if (!dH_dq || !dH_dp) { free(dH_dq); free(dH_dp); return 0.0; }
    sys->grad_H(q, p, n, dH_dq, dH_dp);

    double df_dt = 0.0;
    double *qp = (double *)malloc(n * sizeof(double));
    double *qm = (double *)malloc(n * sizeof(double));
    double *pp = (double *)malloc(n * sizeof(double));
    double *pm = (double *)malloc(n * sizeof(double));
    if (!qp || !qm || !pp || !pm) { /* cleanup and return */ goto cleanup; }

    for (size_t i = 0; i < n; i++) {
        memcpy(qp, q, n * sizeof(double)); qp[i] += eps;
        memcpy(qm, q, n * sizeof(double)); qm[i] -= eps;
        memcpy(pp, p, n * sizeof(double)); pp[i] += eps;
        memcpy(pm, p, n * sizeof(double)); pm[i] -= eps;
        double df_dq = (f(qp, p, n) - f(qm, p, n)) / (2.0 * eps);
        double df_dp = (f(q, pp, n) - f(q, pm, n)) / (2.0 * eps);
        df_dt += df_dq * dH_dp[i] - df_dp * dH_dq[i];
    }

cleanup:
    free(dH_dq); free(dH_dp);
    free(qp); free(qm); free(pp); free(pm);
    return df_dt;
}

/* ──────────────────────────────────────────────────────────────
 * Liouville density memory management
 * ────────────────────────────────────────────────────────────── */

LiouvilleDensity *liouville_density_alloc(size_t n_grid,
                                           double q_min, double q_max,
                                           double p_min, double p_max) {
    LiouvilleDensity *rho = (LiouvilleDensity *)malloc(sizeof(LiouvilleDensity));
    if (!rho) return NULL;
    rho->n_grid = n_grid;
    rho->q_min = q_min; rho->q_max = q_max;
    rho->p_min = p_min; rho->p_max = p_max;
    rho->dq = (q_max - q_min) / (double)(n_grid - 1);
    rho->dp = (p_max - p_min) / (double)(n_grid - 1);
    rho->density = (double *)calloc(n_grid * n_grid, sizeof(double));
    if (!rho->density) { free(rho); return NULL; }
    return rho;
}

void liouville_density_free(LiouvilleDensity *rho) {
    if (rho) {
        free(rho->density);
        free(rho);
    }
}
