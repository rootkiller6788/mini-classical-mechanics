/**
 * second_variation.h -- Second Variation Theory (L4-L6)
 *
 * Legendre condition, Jacobi condition, Weierstrass excess function,
 * conjugate points, sufficient conditions, Ritz and Galerkin methods.
 *
 * References:
 *   Gelfand & Fomin "Calculus of Variations" Ch.5, Ch.6
 *   Bolza "Lectures on the Calculus of Variations" (1904)
 *   Bliss "Lectures on the Calculus of Variations" (1946)
 *
 * Courses: MIT 8.012, Stanford PHYSICS 230, Cambridge Part II,
 *          Berkeley PHYS 242, ETH 402-0800
 */

#ifndef SECOND_VARIATION_H
#define SECOND_VARIATION_H

#include "functional.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Legendre condition */
bool legendre_condition(double d2F_dyp2_val);
bool strong_legendre(double d2F_dyp2_val);
void verify_legendre_along_extremal(
    LagrangianDensity d2F_dyp2,
    const double *x, const double *y, const double *yp,
    int n_nodes, bool *all_pass, double *min_val, double *max_val,
    double *p_vals_out);

/* Jacobi condition and conjugate points */
int find_conjugate_points(
    double a, double b,
    LagrangianDensity P_func, LagrangianDensity Q_func,
    const double *x_ext, const double *y_ext, const double *yp_ext,
    int n_grid, double *conjugate_pts_out, int max_pts);
bool jacobi_condition(const double *conjugate_pts, int n_pts, double a, double b);
double first_conjugate_point(const double *conjugate_pts, int n_pts, double a);
void solve_jacobi_equation(
    double a, double b,
    LagrangianDensity P_func, LagrangianDensity Q_func,
    const double *x_ext, const double *y_ext, const double *yp_ext,
    int n_grid, double u0prime, double *u_out, double *up_out);

/* Weierstrass E-function */
double weierstrass_excess(
    LagrangianDensity F, LagrangianDensity dF_dyp,
    double x, double y, double p_opt, double p);
bool weierstrass_condition(
    LagrangianDensity F, LagrangianDensity dF_dyp,
    double x, double y, double p_opt,
    const double *p_range, int n_p);
bool verify_weierstrass_along_extremal(
    LagrangianDensity F, LagrangianDensity dF_dyp,
    const double *x, const double *y, const double *yp,
    int n_nodes, const double *p_range, int n_p);

/* Hilbert invariant integral */
double hilbert_invariant_integral(
    LagrangianDensity F, LagrangianDensity dF_dyp,
    const double *x, const double *y, const double *yp,
    double (*p_field)(double, double), int n_nodes);
double hilbert_integral_test(
    LagrangianDensity F, LagrangianDensity dF_dyp,
    double (*p_field)(double, double),
    const double *x, const double *y1, const double *yp1,
    const double *y2, const double *yp2, int n_nodes);

/* Sufficient conditions */
typedef struct {
    bool legendre_ok; bool jacobi_ok;
    bool weierstrass_ok; bool sufficient;
} SufficientConditions;

SufficientConditions check_sufficient_conditions(
    LagrangianDensity F, LagrangianDensity dF_dyp,
    LagrangianDensity d2F_dyp2,
    LagrangianDensity P_func, LagrangianDensity Q_func,
    const double *x, const double *y, const double *yp,
    int n_nodes, const double *p_range, int n_p);

/* Direct methods: Ritz, Galerkin */
int ritz_method(
    LagrangianDensity d2F_dyp2, LagrangianDensity d2F_dy2,
    double a, double b,
    int n_basis,
    double (**basis)(double),
    double (**basis_deriv)(double),
    int n_quad, double *coeffs_out);

int galerkin_method(
    double (*el_operator)(double (*u)(double), double x),
    double a, double b,
    int n_basis,
    double (**basis)(double),
    int n_quad, double *coeffs_out);

/* Morse theory: conjugate point counting */
int morse_index(const double *conjugate_pts, int n_pts, double a, double b);
int count_jacobi_zeros(const double *u, int n);

#ifdef __cplusplus
}
#endif

#endif /* SECOND_VARIATION_H */
