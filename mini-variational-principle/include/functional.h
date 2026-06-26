/**
 * functional.h -- Functional Analysis Foundations (L1-L4)
 *
 * References:
 *   Gelfand & Fomin "Calculus of Variations" (Dover 2000)
 *   Courant & Hilbert "Methods of Mathematical Physics" Vol. I
 *   Goldstein, Poole & Safko "Classical Mechanics" Ch.2 (2002)
 *
 * Course alignment:
 *   MIT 8.012, Stanford PHYSICS 230, Cambridge Part II,
 *   Caltech Ph 106, Berkeley PHYS 242, ETH 402-0800
 */

#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Lagrangian density: F: R^3 -> R, (x, y, y') |-> F(x, y, y') */
typedef double (*LagrangianDensity)(double x, double y, double yp);

/** Discrete representation of a functional J[y] = integral F(x,y,y') dx */
typedef struct {
    int     n_nodes;
    double *x;
    double *y;
    double *yp;
    double  a, b;
} FunctionalDiscretization;

/** Functional evaluation context with partial derivatives */
typedef struct {
    LagrangianDensity F;
    LagrangianDensity dF_dy;
    LagrangianDensity dF_dyp;
    LagrangianDensity d2F_dyp2;
    double a, b;
} FunctionalContext;

/* === L2: Functional evaluation === */
double evaluate_functional(LagrangianDensity F, const double *x,
    const double *y, const double *yp, int n_nodes);
double evaluate_functional_simpson(LagrangianDensity F, const double *x,
    const double *y, const double *yp, int n_nodes);
double evaluate_functional_gauss4(LagrangianDensity F, const double *x,
    const double *y, const double *yp, int n_nodes);
double functional_norm(const double *y, const double *yp,
    double a, double b, int n_nodes);
double functional_L2_norm(const double *y, double a, double b, int n_nodes);
double functional_increment(LagrangianDensity F, const double *x,
    const double *y, const double *yp,
    const double *eta, const double *etap, int n_nodes);
void functional_energy_decomposition(LagrangianDensity F, const double *x,
    const double *y, const double *yp, int n_nodes,
    double *kinetic, double *potential);

/* === L3: Functional derivatives === */
double gateaux_derivative(LagrangianDensity F, const double *x,
    const double *y, const double *yp,
    const double *eta, const double *etap, int n_nodes, double epsilon);
void functional_derivative(LagrangianDensity dF_dy, LagrangianDensity dF_dyp,
    const double *x, const double *y, const double *yp,
    int n_nodes, double *result);
double functional_derivative_L2_norm(LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp, const double *x, const double *y,
    const double *yp, int n_nodes);
double frechet_derivative_test(LagrangianDensity F, const double *x,
    const double *y, const double *yp,
    const double *eta, const double *etap, int n_nodes);
void frechet_derivative_convergence(LagrangianDensity F, const double *x,
    const double *y, const double *yp,
    const double *eta, const double *etap,
    int n_nodes, int n_steps, double *rates);

/* === L4: Beltrami identity (first integrals) === */
double beltrami_identity(LagrangianDensity F, LagrangianDensity dF_dyp,
    double x, double y, double yp);
double verify_beltrami(LagrangianDensity F, LagrangianDensity dF_dyp,
    const double *x, const double *y, const double *yp,
    int n_nodes, double *values_out);
double cyclic_coordinate_integral(LagrangianDensity dF_dyp,
    double x, double y, double yp);

/* === L4: Weak and strong extrema === */
int weak_extremum_type(LagrangianDensity F, const double *x,
    const double *y, const double *yp,
    int n_nodes, double delta, int n_perturbations);
void compare_functionals(LagrangianDensity F, const double *x,
    const double *y1, const double *yp1,
    const double *y2, const double *yp2,
    int n_nodes, double *J1, double *J2, double *diff, double *rel_diff);

/* === L4: First and second variations === */
double first_variation_formula(LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp, const double *x,
    const double *y, const double *yp,
    const double *eta, const double *etap, int n_nodes);
double second_variation_formula(LagrangianDensity d2F_dyp2,
    LagrangianDensity d2F_dy2, const double *x,
    const double *eta, const double *etap, int n_nodes);
double second_variation_jacobi_form(LagrangianDensity d2F_dyp2,
    LagrangianDensity d2F_dy2, LagrangianDensity d2F_dydyp,
    const double *x, const double *y, const double *yp,
    const double *eta, const double *etap, int n_nodes);

/* === L3: Convexity === */
bool is_convex_functional(LagrangianDensity F, const double *x,
    const double *y1, const double *yp1,
    const double *y2, const double *yp2,
    int n_nodes, int n_theta);
bool is_quadratic_functional(LagrangianDensity F, const double *x,
    const double *y, const double *yp, int n_nodes);

/* === L5: Utility functions === */
void generate_uniform_nodes(double a, double b, int n, double *x);
void generate_chebyshev_nodes(double a, double b, int n, double *x);
void evaluate_function_samples(double (*y_func)(double),
    double (*yp_func)(double), const double *x, int n,
    double *y, double *yp);
void initialize_functional_context(FunctionalContext *ctx,
    LagrangianDensity F, LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp, LagrangianDensity d2F_dyp2,
    double a, double b);
void compute_derivative_fd(double (*y_func)(double),
    const double *x, int n, double *yp);
double linear_interpolate(const double *x, const double *y,
    int n, double x_query);

#ifdef __cplusplus
}
#endif

#endif /* FUNCTIONAL_H */
