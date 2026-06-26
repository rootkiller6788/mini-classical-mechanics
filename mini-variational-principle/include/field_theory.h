/**
 * field_theory.h -- Field Theory Variational Principles (L6-L8)
 */
#ifndef FIELD_THEORY_H
#define FIELD_THEORY_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
double ginzburg_landau_free_energy(const double *psi, const double *dpsi,
    double a, double b, double c, double x0, double x1, int n);
void solve_ginzburg_landau(double a, double b, double c, double L,
    double bcL, double bcR, int nx, int niter, double omega,
    double *x, double *psi);
double gl_coherence_length(double a, double c);
double gl_critical_field(double a, double b);
void allen_cahn_rhs(const double *phi, double eps, double M,
    double dx, int n, double *dp);
void cahn_hilliard_rhs(const double *phi, double eps, double M,
    double dx, int n, double *dp);
double interface_energy(const double *phi, double eps, double dx, int n);
double helmholtz_free_energy(const double *strain, const double *C,
    double alpha, double dT, double lam, double mu);
double gibbs_free_energy(const double *stress, const double *S);
double em_action_density(const double *E, const double *B,
    double rho, const double *J, const double *A, double phi);
void stress_energy_tensor(
    double (*L)(const double *f, const double *g, int n),
    const double *f, const double *g, int n, double *T);
double topological_charge_density(const double *phi, double dx, int n);
#ifdef __cplusplus
}
#endif
#endif
