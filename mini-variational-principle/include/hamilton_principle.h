/**
 * hamilton_principle.h -- Hamilton's Variational Principle
 */
#ifndef HAMILTON_PRINCIPLE_H
#define HAMILTON_PRINCIPLE_H
#include <stddef.h>
#include <stdbool.h>
#ifdef __cplusplus
extern "C" {
#endif
double phase_space_action(
    double (*H)(const double *q, const double *p, int n),
    const double *t, const double *q, const double *qdot,
    const double *p, int n_dof, int n_steps);
void phase_space_el_residual(
    double (*H)(const double *q, const double *p, int n),
    void (*grad_H_q)(const double *q, const double *p, int n, double *out),
    void (*grad_H_p)(const double *q, const double *p, int n, double *out),
    const double *q, const double *p, const double *qdot, const double *pdot,
    int n, double *residual);
void discrete_hamilton_principle(
    double (*L_d)(const double *qk, const double *qk1, double h),
    const double *qs, int n_steps, int n_dof, double h, double *residuals);
double symplectic_2form(const double *dxi, const double *deta, int n);
double momentum_map(const double *G, const double *q, const double *p, int n);
bool symplectic_reduction(double J_val, double mu, double tol);
void hamiltonian_vector_field(
    void (*grad_H)(const double *q, const double *p, int n, double *out),
    const double *q, const double *p, int n, double *out);
double poisson_bracket(
    double (*F)(const double *q, const double *p, int n),
    double (*G)(const double *q, const double *p, int n),
    const double *q, const double *p, int n, double eps);
double casimir_check(
    double (*C)(const double *q, const double *p, int n),
    double (*H)(const double *q, const double *p, int n),
    const double *q, const double *p, int n, double eps);
void multisymplectic_residual(
    const double *K, const double *M,
    const double *z, const double *dz_dt, const double *dz_dx,
    int n, double *residual);
void least_action_path(
    double (*L)(double t, const double *q, const double *qdot, int n),
    void (*gqL)(double t, const double *q, const double *qdot, int n, double *o),
    void (*gqdL)(double t, const double *q, const double *qdot, int n, double *o),
    const double *q0, const double *qT,
    double t0, double tF, int n_dof, int n_steps,
    double *q_path, int max_iter, double lr);
double noether_charge_energy(
    double (*L)(double t, const double *q, const double *qdot, int n),
    double t, const double *q, const double *qdot, int n, double eps);
#ifdef __cplusplus
}
#endif
#endif
