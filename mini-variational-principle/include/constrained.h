/**
 * constrained.h -- Constrained Variational Problems (L4-L6)
 *
 * Isoperimetric problems, Lagrange multipliers, variable endpoints,
 * transversality conditions, obstacle problems, Bolza form.
 *
 * References:
 *   Gelfand & Fomin "Calculus of Variations" Ch.3, Ch.7
 *   Bliss "Lectures on the Calculus of Variations" (1946)
 *   Goldstein, Poole & Safko "Classical Mechanics" Ch.2 (2002)
 *
 * Courses: MIT 8.012, Stanford PHYSICS 230, Cambridge Part II
 */

#ifndef CONSTRAINED_H
#define CONSTRAINED_H

#include "functional.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * L4: Isoperimetric Problems
 * ============================================================ */

/**
 * isoperimetric_lagrange_multiplier -- Find Lagrange multiplier
 * for isoperimetric constraint: integral G(x,y,y') dx = c
 *
 * Theorem: For the problem min J[y] = integral F dx subject to
 * K[y] = integral G dx = c, the extremal satisfies EL for
 * augmented Lagrangian F* = F + lambda*G.
 *
 * Returns the lambda that satisfies the constraint.
 */
double isoperimetric_lagrange_multiplier(
    LagrangianDensity F, LagrangianDensity G,
    double c, double a, double b,
    double ya, double yb, double lambda_guess);

/**
 * dido_solution -- Dido's problem: maximum area with fixed perimeter
 *
 * Classical isoperimetric problem: maximize area enclosed by a
 * curve of given length L. Solution: circle of radius R = L/(2*pi).
 *
 * Output: x_vals[n], y_vals[n] on circle, *R_out = L/(2*pi)
 */
void dido_solution(double L, int n, double *x_vals, double *y_vals,
                   double *R_out);

/**
 * lagrange_multiplier_functional -- General Lagrange multiplier setup
 *
 * For m constraints integral G_i dx = c_i:
 *   L = F + sum_i lambda_i * G_i
 *
 * Output lambdas_out[m] satisfies all constraints (simplified).
 */
void lagrange_multiplier_functional(
    LagrangianDensity F,
    LagrangianDensity *constraints_G, double *constraints_c,
    int m,
    double a, double b, double ya, double yb,
    double *lambdas_out);

/**
 * solve_constrained_el_bvp -- Solve constrained EL BVP
 *
 * Uses augmented Lagrangian F + lambda*G with shooting method.
 * Returns the Lagrange multiplier lambda.
 */
double solve_constrained_el_bvp(
    LagrangianDensity F, LagrangianDensity G,
    double c, double a, double b,
    double ya, double yb, int n_grid, double tol);

/* ============================================================
 * L4: Variable Endpoint Conditions
 * ============================================================ */

/**
 * variable_left_endpoint_condition -- Natural BC at free endpoint
 *
 * If y(a) is free: dF/dy'|_{x=a} = 0
 */
double variable_left_endpoint_condition(
    LagrangianDensity F, LagrangianDensity dF_dyp,
    double a, double y_a, double yp_a);

/**
 * transversality_diff -- Transversality when endpoint moves
 *
 * (F - y'*dF/dy')*dx + dF/dy'*dy = 0
 * at a variable endpoint (x_end, y_end).
 */
double transversality_diff(
    LagrangianDensity F, LagrangianDensity dF_dyp,
    double x_end, double y_end, double yp_end,
    double dx, double dy);

/* ============================================================
 * L3: Multi-point and Obstacle Problems
 * ============================================================ */

/**
 * multi_point_el_system -- BVP with interior points
 *
 * Solve EL on [a,b] with conditions at interior points.
 * Corner conditions apply at each interior point.
 */
double multi_point_el_system(
    LagrangianDensity F,
    double a, double b,
    const double *interior_pts, int n_interior,
    double ya, double yb, int n_grid);

/**
 * obstacle_problem_residual -- Signorini-type obstacle
 *
 * min J[y] s.t. y(x) >= psi(x)
 *
 * Complementarity: y(x) > psi(x) => EL=0
 *                  y(x) = psi(x) => EL >= 0
 */
double obstacle_problem_residual(
    LagrangianDensity F, double (*psi)(double),
    const double *y, const double *yp,
    double x, int idx, double dx);

/* ============================================================
 * L4: Bolza Problem
 * ============================================================ */

/**
 * bolza_boundary_condition -- Bolza form BC
 *
 * min phi(y(a),y(b)) + integral_a^b F(x,y,y') dx
 *
 * BC: dphi/dy(a) + dF/dy'|_a = 0
 *     dphi/dy(b) - dF/dy'|_b = 0
 */
double bolza_boundary_condition(
    double (*phi)(double, double),
    LagrangianDensity dF_dyp,
    double a, double b,
    const double *y, const double *yp, int n_nodes);

/* ============================================================
 * L5: Numerical Methods for Constraints
 * ============================================================ */

/** Holonomic constraint: g(x,y)=0, Lagrange multiplier term */
double holonomic_constraint_el(
    double (*constraint)(double, double),
    double (*dconstraint_dy)(double, double),
    double x, double y, double yp, double lambda);

/** Penalty method: add (rho/2)*g^2 to Lagrangian */
double penalty_method_el(
    LagrangianDensity F,
    double (*constraint)(double, double),
    double x, double y, double yp, double penalty_param);

/** Augmented Lagrangian: L = F + lambda*g + (rho/2)*g^2 */
double augmented_lagrangian_el(
    LagrangianDensity F,
    double (*constraint)(double, double),
    double (*dconstraint_dy)(double, double),
    double x, double y, double yp,
    double lambda, double rho);

/** RMS constraint violation */
double constraint_violation_norm(
    double (*constraint)(double, double),
    const double *x, const double *y, int n_nodes);

/** Check if constraint is satisfied within tolerance */
bool verify_constraint_satisfaction(
    double (*constraint)(double, double),
    const double *x, const double *y, int n_nodes, double tol);

/** Identify active set for inequality constraints */
int active_set_el(
    double (*constraint)(double, double),
    const double *x, const double *y, int n_nodes,
    double tol, int *active_set);

#ifdef __cplusplus
}
#endif

#endif /* CONSTRAINED_H */
