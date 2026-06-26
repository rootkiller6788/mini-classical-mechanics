/**
 * numerical_variational.h -- Numerical Methods for Variational Problems (L5)
 */
#ifndef NUMERICAL_VARIATIONAL_H
#define NUMERICAL_VARIATIONAL_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
void fem_stiffness_1d(const double *xn, int n, double (*D)(double), double *K);
void fem_mass_1d(const double *xn, int n, double (*rho)(double), double *M);
double integrate_over_element(double (*f)(double), double h, int ng);
void rayleigh_ritz_solve(const double *K, const double *f, int n,
    const int *bc_dofs, const double *bc_vals, int nbc, double *u);
int generalized_eigenvalue(const double *K, const double *M, int n,
    int nm, double *evals, double *evecs);
void spectral_galerkin_fourier(double (*op)(int m, double x, double L),
    int nm, double L, double *Km);
void variational_fd_laplacian(int n, double dx, double *K);
void obstacle_problem_solve(const double *K, const double *f,
    const double *psi, int n, int ni, double *u);
double zienkiewicz_zhu_error_estimate(const double *fe, const double *sr, int n);
double h_refinement_criterion(double err, double tol, double h, int p);
double vms_stabilization_parameter(double h, double kappa);
double vms_advection_diffusion(double Pe, double h);
int morse_index_estimate(const double *H, int n, double tol);
void newton_kantorovich_step(
    void (*F)(const double *u, double *r),
    void (*J)(const double *u, double *jac),
    double *u, int n, double tol, int mi);
double adaptive_variational_stepsize(double eest, double tol, double h);
#ifdef __cplusplus
}
#endif
#endif
