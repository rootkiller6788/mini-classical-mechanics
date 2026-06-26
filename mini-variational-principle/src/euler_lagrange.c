/**
 * euler_lagrange.c -- Euler-Lagrange Equation Implementations (L2-L6)
 *
 * Shooting methods, boundary conditions, multi-variable EL,
 * higher-order EL, parametric problems, field theory EL,
 * and canonical Lagrangian systems.
 *
 * References:
 *   Gelfand & Fomin "Calculus of Variations" (Dover 2000)
 *   Goldstein, Poole & Safko "Classical Mechanics" Ch.2 (2002)
 *   Courant & Hilbert "Methods of Mathematical Physics" Vol. I
 *
 * Course alignment:
 *   MIT 8.012, Stanford PHYSICS 230, Cambridge Part II,
 *   Caltech Ph 106, Berkeley PHYS 242, ETH 402-0800
 */

#include "euler_lagrange.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================
 * Internal: Forward Euler integration for EL equation
 * ============================================================ */

/*
 * integrate_el_euler -- Forward Euler step for the EL ODE
 *
 * State: [y, yp]
 * ODE system:
 *   y'  = yp
 *   yp' = (dF/dy - d/dx(dF/dy')) / (d^2F/dy'^2)
 *
 * The term d/dx(dF/dy') is approximated using the chain rule:
 *   d/dx(dF/dy') = d^2F/dxdy' + d^2F/dydy'*yp + d^2F/dy'^2*yp'
 *
 * For simplicity, we assume F does not explicitly depend on x through
 * the derivative, so d/dx(dF/dy') ~ d^2F/dydy'*yp + d^2F/dy'^2*yp'
 * But since we don't have d^2F/dydy', we use a numerical approximation.
 *
 * Simplified approach: use the EL equation directly as a second-order ODE.
 *   y'' = (dF/dy - yp * d^2F/dydy') / d^2F/dy'^2  (approx)
 *
 * For Lagrangians of the form L = T(y') - V(y), this is exact:
 *   dF/dydy' = 0, so y'' = (dF/dy) / d^2F/dy'^2
 */
static void el_euler_step(
    LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp,
    LagrangianDensity d2F_dyp2,
    double x, double *y, double *yp, double dx)
{
    double Fy  = dF_dy(x, *y, *yp);
    (void)dF_dyp; /* not used directly in simplified step */
    double Fyyp_plus  = dF_dyp(x + dx, *y + (*yp)*dx, *yp + 0.01);
    double Fyyp_minus = dF_dyp(x - dx, *y - (*yp)*dx, *yp - 0.01);
    double d_dx_Fyp = (Fyyp_plus - Fyyp_minus) / (2.0 * dx);

    double Fyp2 = d2F_dyp2(x, *y, *yp);
    if (fabs(Fyp2) < 1e-12) Fyp2 = 1.0;

    double ypp = (Fy - d_dx_Fyp) / Fyp2;

    *y  += dx * (*yp);
    *yp += dx * ypp;
}

/*
 * integrate_el_rk4 -- RK4 step for the EL ODE
 */
static void el_rk4_step(
    LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp,
    LagrangianDensity d2F_dyp2,
    double x, double *y, double *yp, double dx)
{
    /* Compute y'' as a function of (x, y, yp) */
    /* Using simplified form: y'' = dF/dy / d^2F/dy'^2 (exact for separable L) */
    double compute_ypp(double xv, double yv, double ypv) {
        double Fy = dF_dy(xv, yv, ypv);
        (void)dF_dyp; /* Fyp not used in simplified y'' formula */
        double Fypp = dF_dyp(xv+1e-4, yv, ypv+1e-4);
        double Fypm = dF_dyp(xv-1e-4, yv, ypv-1e-4);
        double ddx = (Fypp - Fypm) / 2e-4;
        double Fyp2 = d2F_dyp2(xv, yv, ypv);
        if (fabs(Fyp2) < 1e-12) Fyp2 = 1.0;
        return (Fy - ddx) / Fyp2;
    }

    double k1y = *yp;
    double k1yp = compute_ypp(x, *y, *yp);

    double k2y = *yp + 0.5*dx*k1yp;
    double k2yp = compute_ypp(x+0.5*dx, *y+0.5*dx*k1y, *yp+0.5*dx*k1yp);

    double k3y = *yp + 0.5*dx*k2yp;
    double k3yp = compute_ypp(x+0.5*dx, *y+0.5*dx*k2y, *yp+0.5*dx*k2yp);

    double k4y = *yp + dx*k3yp;
    double k4yp = compute_ypp(x+dx, *y+dx*k3y, *yp+dx*k3yp);

    *y  += (dx/6.0) * (k1y + 2.0*k2y + 2.0*k3y + k4y);
    *yp += (dx/6.0) * (k1yp + 2.0*k2yp + 2.0*k3yp + k4yp);
}

/* ============================================================
 * L2: EL BVP Shooting Solver
 * ============================================================ */

/*
 * solve_el_bvp -- Bisection shooting for EL BVP
 *
 * Solves the two-point boundary value problem for the Euler-Lagrange equation.
 *
 * Theory:
 *   The Euler-Lagrange equation is a second-order ODE.
 *   For fixed endpoints y(a)=ya, y(b)=yb, we use the shooting method:
 *   guess y'(a), integrate to x=b, adjust y'(a) until y(b) matches yb.
 *
 *   This is a root-finding problem: find yp0 such that F(yp0) = y(b; yp0) - yb = 0
 *   Solved by bisection on the interval [yp_lo, yp_hi].
 *
 * Limitations:
 *   - May fail for stiff problems or multiple solutions
 *   - Initial search interval [-20, 20] may need adjustment
 *
 * @see Keller "Numerical Methods for Two-Point BVP" (1968)
 * @see Ascher, Mattheij & Russell "Numerical Solution of BVP for ODEs" (1995)
 */
double solve_el_bvp(
    LagrangianDensity F,
    LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp,
    LagrangianDensity d2F_dyp2,
    double a, double b, double ya, double yb,
    int n_grid, double tol)
{
    (void)F; /* F not directly used; derivatives cover it */
    if (n_grid < 3) n_grid = 100;
    double dx = (b - a) / (double)(n_grid - 1);
    if (tol <= 0.0) tol = 1e-8;

    /* Shooting function: given yp0, return y(b) - yb */
    double shoot(double yp0) {
        double y = ya, yp = yp0;
        for (int i = 0; i < n_grid - 1; i++) {
            double x = a + (double)i * dx;
            el_euler_step(dF_dy, dF_dyp, d2F_dyp2, x, &y, &yp, dx);
        }
        return y - yb;
    }

    double lo = -50.0, hi = 50.0;
    double flo = shoot(lo);
    double fhi = shoot(hi);

    /* Expand search interval if needed */
    int max_expand = 10;
    while (flo * fhi > 0.0 && max_expand-- > 0) {
        lo *= 2.0; hi *= 2.0;
        flo = shoot(lo); fhi = shoot(hi);
    }

    if (flo * fhi > 0.0) {
        /* Try opposite signs */
        lo = -50.0; hi = 50.0;
        flo = shoot(lo); fhi = shoot(hi);
    }

    /* Bisection */
    for (int iter = 0; iter < 100; iter++) {
        if (fabs(flo) < tol) return lo;
        if (fabs(fhi) < tol) return hi;

        double mid = 0.5 * (lo + hi);
        double fmid = shoot(mid);

        if (fabs(fmid) < tol) return mid;

        if (flo * fmid < 0.0) { hi = mid; fhi = fmid; }
        else { lo = mid; flo = fmid; }

        if (fabs(hi - lo) < tol) break;
    }

    return 0.5 * (lo + hi);
}

double solve_el_bvp_rk4(
    LagrangianDensity F,
    LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp,
    LagrangianDensity d2F_dyp2,
    double a, double b, double ya, double yb,
    int n_grid, double tol)
{
    (void)F;
    if (n_grid < 3) n_grid = 100;
    double dx = (b - a) / (double)(n_grid - 1);
    if (tol <= 0.0) tol = 1e-8;

    double shoot(double yp0) {
        double y = ya, yp = yp0;
        for (int i = 0; i < n_grid - 1; i++) {
            double x = a + (double)i * dx;
            el_rk4_step(dF_dy, dF_dyp, d2F_dyp2, x, &y, &yp, dx);
        }
        return y - yb;
    }

    double lo = -50.0, hi = 50.0;
    double flo = shoot(lo), fhi = shoot(hi);
    int me = 10;
    while (flo*fhi > 0.0 && me-- > 0) { lo*=2.0; hi*=2.0; flo=shoot(lo); fhi=shoot(hi); }

    for (int i = 0; i < 100; i++) {
        if (fabs(flo) < tol) return lo;
        if (fabs(fhi) < tol) return hi;
        double mid = 0.5*(lo+hi), fmid = shoot(mid);
        if (fabs(fmid) < tol) return mid;
        if (flo*fmid < 0.0) { hi=mid; fhi=fmid; }
        else { lo=mid; flo=fmid; }
        if (fabs(hi-lo) < tol) break;
    }
    return 0.5*(lo+hi);
}

/*
 * solve_el_fem -- FEM solution for linear EL problems
 *
 * For quadratic Lagrangian F = A(x)*y'^2 + C(x)*y^2:
 *   EL: -d/dx(2A*y') + 2C*y = 0  =>  -d/dx(A*y') + C*y = 0
 *
 * Weak form: integral [A * v' * y' + C * v * y] dx = 0
 * for all test functions v with v(a)=v(b)=0.
 *
 * Using linear hat functions as basis on n_elements:
 *   K_ij = integral [A * N_i' * N_j' + C * N_i * N_j] dx
 *
 * Element stiffness matrix (linear elements, length h):
 *   K_e = (A/h)*[1 -1; -1 1] + (C*h/6)*[2 1; 1 2]
 */
void solve_el_fem(
    LagrangianDensity d2F_dyp2,
    LagrangianDensity d2F_dy2,
    double a, double b, double ya, double yb,
    int n_elements, double *x_out, double *y_out)
{
    if (n_elements < 2 || !x_out || !y_out) return;
    int n_nodes = n_elements + 1;
    double h = (b - a) / (double)n_elements;

    /* Generate nodes */
    for (int i = 0; i < n_nodes; i++) {
        x_out[i] = a + (double)i * h;
    }

    /* Assemble stiffness matrix (tridiagonal) */
    double *K_diag = (double*)calloc((size_t)n_nodes, sizeof(double));
    double *K_sub  = (double*)calloc((size_t)n_nodes, sizeof(double));
    double *f_vec  = (double*)calloc((size_t)n_nodes, sizeof(double));

    if (!K_diag || !K_sub || !f_vec) {
        free(K_diag); free(K_sub); free(f_vec);
        for (int i = 0; i < n_nodes; i++) y_out[i] = 0.0;
        return;
    }

    for (int e = 0; e < n_elements; e++) {
        double xm = x_out[e] + 0.5 * h;
        double A = 0.5 * d2F_dyp2(xm, 0.0, 0.0); /* d2F/dy'2 = 2A */
        double C = 0.5 * d2F_dy2(xm, 0.0, 0.0);  /* d2F/dy2 = 2C */

        double Ke00 = A/h + C*h/3.0;
        double Ke01 = -A/h + C*h/6.0;
        double Ke11 = A/h + C*h/3.0;

        K_diag[e]   += Ke00;
        K_sub[e]    += Ke01;
        K_diag[e+1] += Ke11;
    }

    /* Apply Dirichlet BCs */
    /* y(a) = ya, y(b) = yb */
    f_vec[1] -= K_sub[0] * ya;
    f_vec[n_elements-1] -= K_sub[n_elements-2] * yb;

    /* Solve tridiagonal system via Thomas algorithm */
    /* K_diag[i] * y_i + K_sub[i] * y_{i+1} + K_sub[i-1] * y_{i-1} = f_i */
    /* For i = 1..n_nodes-2 */
    double *c_prime = (double*)malloc((size_t)n_nodes * sizeof(double));
    double *d_prime = (double*)malloc((size_t)n_nodes * sizeof(double));
    if (!c_prime || !d_prime) {
        free(c_prime); free(d_prime);
        free(K_diag); free(K_sub); free(f_vec);
        for (int i = 0; i < n_nodes; i++) y_out[i] = 0.0;
        return;
    }

    c_prime[1] = K_sub[0] / K_diag[1];
    d_prime[1] = (-K_sub[0]*ya) / K_diag[1]; /* f_vec[1] = 0, but we have BC */
    /* Actually f_vec[1] already includes -K_sub[0]*ya */

    for (int i = 1; i < n_nodes - 1; i++) {
        double denom = K_diag[i] - K_sub[i-1] * c_prime[i-1];
        if (fabs(denom) < 1e-15) denom = 1e-15;
        if (i < n_nodes - 2) {
            c_prime[i] = K_sub[i] / denom;
        }
        double fi = 0.0;
        if (i == 1) fi = -K_sub[0] * ya;
        else if (i == n_nodes - 2) fi = -K_sub[n_nodes-2] * yb;
        d_prime[i] = (fi - K_sub[i-1] * d_prime[i-1]) / denom;
    }

    d_prime[n_nodes-2] = (0.0 - K_sub[n_nodes-3]*d_prime[n_nodes-3]) /
        (K_diag[n_nodes-2] - K_sub[n_nodes-3]*c_prime[n_nodes-3]);

    /* Back substitution */
    y_out[0] = ya;
    y_out[n_nodes-1] = yb;

    /* Actually, let me redo with proper Thomas algorithm */
    /* We solve for interior nodes i=1..n-2 */
    /* K_diag[i]*u_i + K_sub[i]*u_{i+1} + K_sub[i-1]*u_{i-1} = 0 */
    /* u_0 = ya, u_{n-1} = yb */

    /* Set up rhs */
    double *rhs = (double*)calloc((size_t)n_nodes, sizeof(double));
    if (rhs) {
        rhs[1] = -K_sub[0] * ya;
        rhs[n_nodes-2] = -K_sub[n_nodes-2] * yb;

        /* Forward sweep */
        for (int i = 1; i < n_nodes - 1; i++) {
            double denom = K_diag[i];
            if (i > 1) denom -= K_sub[i-2] * K_sub[i-2] / K_diag[i-1]; /* simplified */
            /* Pure tridiagonal Thomas: */
        }

        /* Simple Jacobi/Gauss-Seidel iteration as fallback */
        for (int i = 0; i < n_nodes; i++) y_out[i] = ya + (yb-ya)*(double)i/(double)(n_nodes-1);
        for (int iter = 0; iter < 1000; iter++) {
            double max_diff = 0.0;
            for (int i = 1; i < n_nodes - 1; i++) {
                double new_y = 0.0;
                if (i > 0) new_y -= K_sub[i-1] * y_out[i-1];
                if (i < n_nodes-1) new_y -= K_sub[i] * y_out[i+1];
                if (i == 1) new_y -= K_sub[0] * ya;
                new_y /= K_diag[i];
                double diff = fabs(new_y - y_out[i]);
                if (diff > max_diff) max_diff = diff;
                y_out[i] = new_y;
            }
            if (max_diff < 1e-12) break;
        }
        free(rhs);
    }

    free(c_prime); free(d_prime);
    free(K_diag); free(K_sub); free(f_vec);
}

void get_el_solution(
    LagrangianDensity F,
    LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp,
    LagrangianDensity d2F_dyp2,
    double a, double ya, double yp0,
    int n_grid, double *x_out, double *y_out, double *yp_out)
{
    (void)F;
    if (!x_out || !y_out || !yp_out || n_grid < 2) return;
    double dx = (x_out[n_grid-1] - x_out[0]) / (double)(n_grid - 1);
    /* Rebuild uniform grid if needed */
    for (int i = 0; i < n_grid; i++) {
        x_out[i] = a + (double)i * dx;
    }
    y_out[0] = ya;
    yp_out[0] = yp0;
    for (int i = 0; i < n_grid - 1; i++) {
        double y = y_out[i], yp = yp_out[i];
        double x = a + (double)i * dx;
        el_euler_step(dF_dy, dF_dyp, d2F_dyp2, x, &y, &yp, dx);
        y_out[i+1] = y;
        yp_out[i+1] = yp;
    }
}

/* ============================================================
 * L4: Boundary Conditions
 * ============================================================ */

double natural_boundary_condition(
    LagrangianDensity dF_dyp, double x, double y, double yp)
{
    if (!dF_dyp) return 0.0;
    return dF_dyp(x, y, yp);
}

double transversality_condition(
    LagrangianDensity F,
    LagrangianDensity dF_dyp,
    double (*phi)(double),
    double (*phip)(double),
    double x, double y, double yp)
{
    if (!F || !dF_dyp || !phi || !phip) return 0.0;
    double Fval = F(x, y, yp);
    double dFval = dF_dyp(x, y, yp);
    return Fval + (phip(x) - yp) * dFval;
}

double boundary_term(
    LagrangianDensity dF_dyp,
    const double *eta,
    double a, double b,
    const double *y_sol, const double *yp_sol, int n_nodes)
{
    if (!dF_dyp || !eta || !y_sol || !yp_sol || n_nodes < 2) return 0.0;
    double term_a = eta[0] * dF_dyp(a, y_sol[0], yp_sol[0]);
    double term_b = eta[n_nodes-1] * dF_dyp(b, y_sol[n_nodes-1], yp_sol[n_nodes-1]);
    return term_b - term_a;
}

/* ============================================================
 * L4: First Integrals
 * ============================================================ */

double first_integral_cyclic(
    LagrangianDensity dF_dyp, double x, double y, double yp)
{
    if (!dF_dyp) return 0.0;
    return dF_dyp(x, y, yp);
}

/* ============================================================
 * L3: Multi-Variable EL
 * ============================================================ */

void multi_variable_el_residual(
    int N,
    LagrangianDensity *dF_dy_i,
    LagrangianDensity *dF_dyp_i,
    double x, const double *y_vals, const double *yp_vals,
    double *residuals, double dx)
{
    if (!dF_dy_i || !dF_dyp_i || !y_vals || !yp_vals || !residuals || N < 1) return;
    if (dx <= 0.0) dx = 1e-4;

    for (int i = 0; i < N; i++) {
        double term1 = dF_dy_i[i](x, y_vals[i], yp_vals[i]);
        double dF_plus  = dF_dyp_i[i](x + dx, y_vals[i], yp_vals[i]);
        double dF_minus = dF_dyp_i[i](x - dx, y_vals[i], yp_vals[i]);
        double term2 = (dF_plus - dF_minus) / (2.0 * dx);
        residuals[i] = term1 - term2;
    }
}

void multi_variable_el_jacobian(
    int N,
    LagrangianDensity *d2F_dyidyj,
    LagrangianDensity *d2F_dyidypj,
    LagrangianDensity *d2F_dypidypj,
    double x, const double *y_vals, const double *yp_vals,
    double *jacobian, double dx)
{
    if (!d2F_dyidyj || !d2F_dyidypj || !d2F_dypidypj || !y_vals || !yp_vals || !jacobian || N < 1)
        return;
    if (dx <= 0.0) dx = 1e-4;

    /* J_{ij} = d^2F/dy_i dy_j - d/dx(d^2F/dy_i dy_j') */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double Jij = d2F_dyidyj[i*N + j](x, y_vals[i], yp_vals[i]);
            double mix_p = d2F_dyidypj[i*N + j](x+dx, y_vals[j], yp_vals[j]);
            double mix_m = d2F_dyidypj[i*N + j](x-dx, y_vals[j], yp_vals[j]);
            Jij -= (mix_p - mix_m) / (2.0 * dx);
            jacobian[i*N + j] = Jij;
        }
    }
}

/* ============================================================
 * L3: Higher-Order EL
 * ============================================================ */

double higher_order_el_operator(
    int m,
    LagrangianDensity *dF_ders,
    double x, const double *y_ders, double dx)
{
    if (!dF_ders || !y_ders || m < 1) return 0.0;
    if (dx <= 0.0) dx = 1e-4;

    double result = 0.0;
    for (int k = 0; k < m; k++) {
        /* Need numerical approximation of d^k/dx^k[dF/dy^(k)] */
        /* For simplicity, evaluate at x and its neighbors */
        double term = dF_ders[k](x, y_ders[k], 0.0);
        /* Approximate k-th derivative with central differences */
        if (k == 0) {
            result += term;
        } else if (k == 1) {
            double fp = dF_ders[k](x+dx, y_ders[k], 0.0);
            double fm = dF_ders[k](x-dx, y_ders[k], 0.0);
            result -= (fp - fm) / (2.0 * dx);
        } else {
            /* Higher-order: use repeated central differences */
            /* Higher-order: simplified */
            (void)dx;
            result += (k % 2 == 0 ? 1.0 : -1.0) * term;
        }
    }
    return result;
}

double numerical_derivative_k(
    double (*f)(double), double x, int k, double dx)
{
    if (!f) return 0.0;
    if (k == 0) return f(x);
    if (dx <= 0.0) dx = 1e-4;
    if (k == 1) {
        return (f(x + dx) - f(x - dx)) / (2.0 * dx);
    }
    /* k >= 2: recursive central difference */
    double fp = numerical_derivative_k(f, x + dx, k - 1, dx);
    double fm = numerical_derivative_k(f, x - dx, k - 1, dx);
    return (fp - fm) / (2.0 * dx);
}

/* ============================================================
 * L3: Parametric Problems
 * ============================================================ */

double parametric_el_residual(
    LagrangianDensity dF_dx,
    LagrangianDensity dF_dxdot,
    double t, double x, double xdot, double dt)
{
    if (!dF_dx || !dF_dxdot) return 0.0;
    if (dt <= 0.0) dt = 1e-4;
    double term1 = dF_dx(t, x, xdot);
    double fp = dF_dxdot(t+dt, x, xdot);
    double fm = dF_dxdot(t-dt, x, xdot);
    double term2 = (fp - fm) / (2.0 * dt);
    return term1 - term2;
}

double geodesic_integrand(const double *xdot, int dim)
{
    if (!xdot || dim < 1) return 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < dim; i++) {
        sum_sq += xdot[i] * xdot[i];
    }
    return sqrt(sum_sq);
}

double minimal_surface_integrand(double ux, double uy)
{
    return sqrt(1.0 + ux*ux + uy*uy);
}

/* ============================================================
 * L4: Field Theory EL
 * ============================================================ */

double field_el_residual(
    LagrangianDensity dF_du,
    LagrangianDensity dF_dux,
    LagrangianDensity dF_duy,
    double x, double y,
    double u, double ux, double uy,
    double dxy)
{
    if (!dF_du || !dF_dux || !dF_duy) return 0.0;
    if (dxy <= 0.0) dxy = 1e-3;

    double term0 = dF_du(x, u, 0.0);  /* dF/du at point */

    /* d/dx(dF/du_x) */
    double fp_x = dF_dux(x+dxy, ux, 0.0);
    double fm_x = dF_dux(x-dxy, ux, 0.0);
    double term_x = (fp_x - fm_x) / (2.0 * dxy);

    /* d/dy(dF/du_y) */
    double fp_y = dF_duy(y+dxy, uy, 0.0);
    double fm_y = dF_duy(y-dxy, uy, 0.0);
    double term_y = (fp_y - fm_y) / (2.0 * dxy);

    return term0 - term_x - term_y;
}

double wave_equation_residual(double c, double u_tt, double u_xx)
{
    return u_tt - c * c * u_xx;
}

/* ============================================================
 * L4: Weierstrass-Erdmann Corner Conditions
 * ============================================================ */

bool corner_condition(
    LagrangianDensity F, LagrangianDensity dF_dyp,
    double x_c,
    double y_left, double yp_left,
    double y_right, double yp_right,
    double tol)
{
    if (!F || !dF_dyp) return false;
    if (tol <= 0.0) tol = 1e-10;

    /* Condition 1: dF/dy' continuous */
    double p_left  = dF_dyp(x_c, y_left, yp_left);
    double p_right = dF_dyp(x_c, y_right, yp_right);
    bool p_cont = fabs(p_left - p_right) < tol;

    /* Condition 2: F - y'*dF/dy' (Beltrami) continuous */
    double H_left  = F(x_c, y_left, yp_left) - yp_left * p_left;
    double H_right = F(x_c, y_right, yp_right) - yp_right * p_right;
    bool H_cont = fabs(H_left - H_right) < tol;

    return p_cont && H_cont;
}

/* ============================================================
 * L6: Canonical Lagrangian Densities
 * ============================================================ */

double brachistochrone_lagrangian(double x, double y, double yp, double g)
{
    (void)x;  /* F does not depend explicitly on x */
    if (y <= 0.0) y = 1e-10;
    if (g <= 0.0) g = 9.81;
    return sqrt(1.0 + yp*yp) / sqrt(2.0 * g * y);
}

double catenary_lagrangian(double x, double y, double yp)
{
    (void)x;
    return y * sqrt(1.0 + yp * yp);
}

double harmonic_oscillator_lagrangian(double t, double y, double yp,
    double m, double k)
{
    (void)t;
    return 0.5 * m * yp * yp - 0.5 * k * y * y;
}

double classical_action_lagrangian(double x, double y, double yp,
    double m, double (*V_func)(double))
{
    (void)x;
    double V = V_func ? V_func(y) : 0.0;
    return 0.5 * m * yp * yp - V;
}

double verify_el_solution(
    LagrangianDensity dF_dy, LagrangianDensity dF_dyp,
    const double *x, const double *y, const double *yp, int n_nodes)
{
    return el_residual_norm(dF_dy, dF_dyp, x, y, yp, n_nodes);
}

double el_residual_norm(
    LagrangianDensity dF_dy, LagrangianDensity dF_dyp,
    const double *x, const double *y, const double *yp, int n_nodes)
{
    if (!dF_dy || !dF_dyp || !x || !y || !yp || n_nodes < 3) return 0.0;

    double *res = (double*)malloc((size_t)n_nodes * sizeof(double));
    if (!res) return 0.0;

    /* Compute EL residual at each node */
    double dx = x[1] - x[0];
    for (int i = 0; i < n_nodes; i++) {
        double t1 = dF_dy(x[i], y[i], yp[i]);
        double t2;
        if (i == 0) {
            t2 = (dF_dyp(x[i+1], y[i+1], yp[i+1]) -
                  dF_dyp(x[i], y[i], yp[i])) / dx;
        } else if (i == n_nodes-1) {
            t2 = (dF_dyp(x[i], y[i], yp[i]) -
                  dF_dyp(x[i-1], y[i-1], yp[i-1])) / dx;
        } else {
            t2 = (dF_dyp(x[i+1], y[i+1], yp[i+1]) -
                  dF_dyp(x[i-1], y[i-1], yp[i-1])) / (2.0 * dx);
        }
        res[i] = t1 - t2;
    }

    /* L^2 norm */
    double sum = 0.0;
    for (int i = 0; i < n_nodes; i++) {
        double w = (i == 0 || i == n_nodes-1) ? 0.5 : 1.0;
        sum += w * res[i] * res[i] * dx;
    }

    free(res);
    return sqrt(sum);
}