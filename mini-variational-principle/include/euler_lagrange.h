/**
 * euler_lagrange.h -- Euler-Lagrange Equations (L1-L4)
 *
 * References:
 *   Gelfand & Fomin "Calculus of Variations" Ch.1-3
 *   Goldstein, Poole & Safko "Classical Mechanics" Ch.2 (2002)
 *   Courant & Hilbert "Methods of Mathematical Physics" Vol. I, Ch.4
 *
 * Course alignment:
 *   MIT 8.012, Stanford PHYSICS 230, Cambridge Part II,
 *   Caltech Ph 106, Berkeley PHYS 242, ETH 402-0800
 */

#ifndef EULER_LAGRANGE_H
#define EULER_LAGRANGE_H

#include "functional.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * L2: 1D Euler-Lagrange BVP Solver
 * ============================================================ */

/**
 * solve_el_bvp -- Shooting method for 1D Euler-Lagrange BVP
 *
 * Solves: d/dx(dF/dy') - dF/dy = 0 with y(a)=ya, y(b)=yb
 *
 * Algorithm: Bisection shooting
 *   1. Guess y'(a) = yp0
 *   2. Integrate EL equation forward to x=b using Forward Euler
 *   3. Compare y(b) with yb
 *   4. Bisect on yp0 until |y(b) - yb| < tol
 */
double solve_el_bvp(
    LagrangianDensity F,
    LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp,
    LagrangianDensity d2F_dyp2,
    double a, double b, double ya, double yb,
    int n_grid, double tol);

/**
 * solve_el_bvp_rk4 -- RK4 shooting method for EL BVP
 *
 * Uses 4th-order Runge-Kutta integration for higher accuracy.
 */
double solve_el_bvp_rk4(
    LagrangianDensity F,
    LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp,
    LagrangianDensity d2F_dyp2,
    double a, double b, double ya, double yb,
    int n_grid, double tol);

/**
 * solve_el_fem -- FEM for linear EL problems
 *
 * For quadratic Lagrangians, EL is linear.
 * Assembles stiffness matrix and solves Ku = f.
 * Output arrays x_out, y_out must be pre-allocated (n_elements+1 elements).
 */
void solve_el_fem(
    LagrangianDensity d2F_dyp2,
    LagrangianDensity d2F_dy2,
    double a, double b, double ya, double yb,
    int n_elements, double *x_out, double *y_out);

/**
 * get_el_solution -- Reconstruct full trajectory from shooting parameter
 *
 * x_out, y_out, yp_out must be pre-allocated with n_grid elements.
 */
void get_el_solution(
    LagrangianDensity F,
    LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp,
    LagrangianDensity d2F_dyp2,
    double a, double ya, double yp0,
    int n_grid, double *x_out, double *y_out, double *yp_out);

/* ============================================================
 * L4: Boundary Conditions
 * ============================================================ */

/** Natural boundary condition: dF/dy' = 0 at free endpoint */
double natural_boundary_condition(
    LagrangianDensity dF_dyp, double x, double y, double yp);

/**
 * Transversality condition when endpoint moves along curve phi(x):
 *   F + (phi'(x) - y') * dF/dy' = 0
 */
double transversality_condition(
    LagrangianDensity F,
    LagrangianDensity dF_dyp,
    double (*phi)(double),
    double (*phip)(double),
    double x, double y, double yp);

/** Boundary term: [eta * dF/dy']_a^b from integration by parts */
double boundary_term(
    LagrangianDensity dF_dyp,
    const double *eta,
    double a, double b,
    const double *y_sol, const double *yp_sol, int n_nodes);

/* ============================================================
 * L4: First Integrals
 * ============================================================ */

/** Cyclic coordinate: if dF/dy=0 then dF/dy'=const (conserved momentum) */
double first_integral_cyclic(
    LagrangianDensity dF_dyp, double x, double y, double yp);

/* ============================================================
 * L3: Multi-Variable Euler-Lagrange
 * ============================================================ */

/**
 * Multi-variable EL residual for N dependent variables.
 * residuals[i] = dF/dy_i - d/dx(dF/dy_i')
 */
void multi_variable_el_residual(
    int N,
    LagrangianDensity *dF_dy_i,
    LagrangianDensity *dF_dyp_i,
    double x, const double *y_vals, const double *yp_vals,
    double *residuals, double dx);

/** Jacobian of multi-variable EL system for Newton iteration */
void multi_variable_el_jacobian(
    int N,
    LagrangianDensity *d2F_dyidyj,
    LagrangianDensity *d2F_dyidypj,
    LagrangianDensity *d2F_dypidypj,
    double x, const double *y_vals, const double *yp_vals,
    double *jacobian, double dx);

/* ============================================================
 * L3: Higher-Order Euler-Lagrange
 * ============================================================ */

/**
 * Higher-order EL operator for F(x, y, y', ..., y^(m)):
 *   sum_{k=0}^m (-1)^k d^k/dx^k (dF/dy^(k)) = 0
 */
double higher_order_el_operator(
    int m,
    LagrangianDensity *dF_ders,
    double x, const double *y_ders, double dx);

/** k-th order numerical derivative via central differences */
double numerical_derivative_k(
    double (*f)(double), double x, int k, double dx);

/* ============================================================
 * L3: Parametric and Geometric Problems
 * ============================================================ */

/** EL residual for parametric curves x(t) */
double parametric_el_residual(
    LagrangianDensity dF_dx,
    LagrangianDensity dF_dxdot,
    double t, double x, double xdot, double dt);

/** Geodesic (shortest path) Lagrangian: sqrt(sum xdot_i^2) */
double geodesic_integrand(const double *xdot, int dim);

/** Minimal surface Lagrangian: sqrt(1 + ux^2 + uy^2) */
double minimal_surface_integrand(double ux, double uy);

/* ============================================================
 * L4: Field Theory Euler-Lagrange
 * ============================================================ */

/** 2D scalar field EL: dF/du - d/dx(dF/du_x) - d/dy(dF/du_y) = 0 */
double field_el_residual(
    LagrangianDensity dF_du,
    LagrangianDensity dF_dux,
    LagrangianDensity dF_duy,
    double x, double y,
    double u, double ux, double uy,
    double dxy);

/** 1D wave equation residual: u_tt - c^2*u_xx */
double wave_equation_residual(double c, double u_tt, double u_xx);

/* ============================================================
 * L4: Weierstrass-Erdmann Corner Conditions
 * ============================================================ */

/**
 * Corner conditions: dF/dy' and F-y'*dF/dy' must be continuous across
 * any corner (discontinuity in y').
 */
bool corner_condition(
    LagrangianDensity F, LagrangianDensity dF_dyp,
    double x_c,
    double y_left, double yp_left,
    double y_right, double yp_right,
    double tol);

/* ============================================================
 * L6: Canonical Lagrangian Densities
 * ============================================================ */

/** Brachistochrone: F = sqrt(1+y'^2)/sqrt(y) */
double brachistochrone_lagrangian(double x, double y, double yp, double g);

/** Catenary: F = y*sqrt(1+y'^2) */
double catenary_lagrangian(double x, double y, double yp);

/** Harmonic oscillator: L = (1/2)*m*y'^2 - (1/2)*k*y^2 */
double harmonic_oscillator_lagrangian(double t, double y, double yp,
    double m, double k);

/** Classical action: L = (1/2)*m*y'^2 - V(y) */
double classical_action_lagrangian(double x, double y, double yp,
    double m, double (*V_func)(double));

/** L^2 norm of EL residual along a curve */
double verify_el_solution(
    LagrangianDensity dF_dy, LagrangianDensity dF_dyp,
    const double *x, const double *y, const double *yp, int n_nodes);

/** RMS norm of EL residual */
double el_residual_norm(
    LagrangianDensity dF_dy, LagrangianDensity dF_dyp,
    const double *x, const double *y, const double *yp, int n_nodes);

#ifdef __cplusplus
}
#endif

#endif /* EULER_LAGRANGE_H */