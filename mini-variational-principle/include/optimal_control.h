/**
 * optimal_control.h -- Optimal Control Theory (L5-L8)
 * Pontryagin Maximum Principle, LQR, Bellman DP, MPC, HJB.
 */
#ifndef OPTIMAL_CONTROL_H
#define OPTIMAL_CONTROL_H
#include <stddef.h>
#include <stdbool.h>
#ifdef __cplusplus
extern "C" {
#endif
double pontryagin_hamiltonian(
    double (*L)(const double *x, const double *u, double t),
    void (*f)(const double *x, const double *u, double t, double *dxdt),
    const double *x, const double *u, const double *p,
    int nx, int nu, double t);
void adjoint_rhs(
    double (*H)(const double *x, const double *u, const double *p, double t),
    const double *x, const double *u, const double *p,
    int n, double t, double dx, double *dp);
double switching_function(
    double (*H)(const double *x, double u, const double *p, double t),
    const double *x, const double *p, double t,
    const double *u_vals, int nu, int nx, double *best_u);
int lqr_controller(const double *A, const double *B, const double *Q,
    const double *R, int n, int m, double *K, double *P,
    int max_iter, double tol);
double lqr_cost(const double *x0, const double *P, int n);
int lqr_tracking_gain(const double *A, const double *B, const double *Q,
    const double *R, int n, int m, double *K, double *P);
void bellman_value_iteration(
    void (*f)(const double *x, const double *u, double *xn),
    double (*g)(const double *x, const double *u),
    const double *xg, int nx, int dx, const double *ug, int nu, int du,
    int ns, double *V, int *policy);
double solve_optimal_control_bvp(
    void (*f)(const double *x, const double *u, double t, double *dxdt),
    double (*L)(const double *x, const double *u, double t),
    const double *x0, double t0, double tF, int n,
    int ng, double *p0_out);
double bang_bang_switch_curve(double x, double v);
double minimum_time_control(double x, double v, double xt, double vt);
double minimum_time_double_integrator(double x0, double v0, double umax);
void mpc_step(const double *A, const double *B, const double *Q,
    const double *R, const double *xc, int n, int m, int hz, double *u);
bool hjb_lqr_verify(const double *P, const double *A, const double *B,
    const double *Q, const double *R, int n, int m, double tol);
#ifdef __cplusplus
}
#endif
#endif
