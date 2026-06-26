/**
 * functional.c -- Functional Analysis Implementation (L1-L5)
 *
 * Implements: functional evaluation (trapezoidal, Simpson, Gauss4),
 * Gateaux/Frechet derivatives, Beltrami identity, first/second
 * variation formulas, convexity tests, and utility functions.
 *
 * References:
 *   Gelfand & Fomin "Calculus of Variations" (Dover, 2000)
 *   Goldstein, Poole & Safko "Classical Mechanics" Ch.2 (2002)
 *
 * Courses: MIT 8.012, Stanford PHYSICS 230, Cambridge Part II,
 *          Caltech Ph 106, Berkeley PHYS 242, ETH 402-0800
 */

#include "functional.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- Internal quadrature helpers ---- */
static double trap_int(const double *f, const double *x, int n) {
    if (n < 2) return 0.0;
    double s = 0.5 * (f[0] + f[n-1]);
    for (int i = 1; i < n-1; i++) s += f[i];
    return s * (x[1] - x[0]);
}
static double simp_int(const double *f, int n, double dx) {
    if (n < 3 || n % 2 == 0) {
        double s = 0.5 * (f[0] + f[n-1]);
        for (int i = 1; i < n-1; i++) s += f[i];
        return s * dx;
    }
    double s = f[0] + f[n-1];
    for (int i = 1; i < n-1; i += 2) s += 4.0 * f[i];
    for (int i = 2; i < n-2; i += 2) s += 2.0 * f[i];
    return s * dx / 3.0;
}

/* ---- L2: Functional Evaluation ---- */
double evaluate_functional(LagrangianDensity F,
    const double *x, const double *y, const double *yp, int n)
{
    if (!F || !x || !y || !yp || n < 2) return 0.0;
    double *fv = (double*)malloc((size_t)n * sizeof(double));
    if (!fv) return 0.0;
    for (int i = 0; i < n; i++) fv[i] = F(x[i], y[i], yp[i]);
    double r = trap_int(fv, x, n);
    free(fv);
    return r;
}
double evaluate_functional_simpson(LagrangianDensity F,
    const double *x, const double *y, const double *yp, int n)
{
    if (!F || !x || !y || !yp || n < 3) return 0.0;
    double dx = x[1] - x[0];
    double *fv = (double*)malloc((size_t)n * sizeof(double));
    if (!fv) return 0.0;
    for (int i = 0; i < n; i++) fv[i] = F(x[i], y[i], yp[i]);
    double r = simp_int(fv, n, dx);
    free(fv);
    return r;
}
double evaluate_functional_gauss4(LagrangianDensity F,
    const double *x, const double *y, const double *yp, int n)
{
    if (!F || !x || !y || !yp || n < 2) return 0.0;
    static const double gp[4] = {
        -0.8611363115940526, -0.3399810435848563,
         0.3399810435848563,  0.8611363115940526};
    static const double gw[4] = {
        0.3478548451374538, 0.6521451548625461,
        0.6521451548625461, 0.3478548451374538};
    double r = 0.0;
    for (int i = 0; i < n - 1; i++) {
        double xm = 0.5 * (x[i] + x[i+1]);
        double hh = 0.5 * (x[i+1] - x[i]);
        double sr = 0.0;
        for (int k = 0; k < 4; k++) {
            double xk = xm + hh * gp[k];
            double yl = y[i] + (y[i+1]-y[i]) * (xk-x[i]) / (x[i+1]-x[i]);
            double ypl = yp[i] + (yp[i+1]-yp[i]) * (xk-x[i]) / (x[i+1]-x[i]);
            sr += gw[k] * F(xk, yl, ypl);
        }
        r += hh * sr;
    }
    return r;
}
double functional_norm(const double *y, const double *yp,
    double a, double b, int n)
{
    if (!y || !yp || n < 2) return 0.0;
    double my = 0.0, myp = 0.0, l2 = 0.0;
    double dx = (b - a) / (double)(n - 1);
    for (int i = 0; i < n; i++) {
        double ay = fabs(y[i]), ayp = fabs(yp[i]);
        if (ay > my) my = ay;
        if (ayp > myp) myp = ayp;
        double w = (i == 0 || i == n - 1) ? 0.5 : 1.0;
        l2 += w * y[i] * y[i] * dx;
    }
    return my + myp + sqrt(l2);
}
double functional_L2_norm(const double *y, double a, double b, int n)
{
    if (!y || n < 2) return 0.0;
    double s = 0.0, dx = (b - a) / (double)(n - 1);
    for (int i = 0; i < n; i++) {
        double w = (i == 0 || i == n - 1) ? 0.5 : 1.0;
        s += w * y[i] * y[i] * dx;
    }
    return sqrt(s);
}
double functional_increment(LagrangianDensity F,
    const double *x, const double *y, const double *yp,
    const double *eta, const double *etap, int n)
{
    if (!F || !x || !y || !yp || !eta || !etap || n < 2) return 0.0;
    double *yn  = (double*)malloc((size_t)n * sizeof(double));
    double *ypn = (double*)malloc((size_t)n * sizeof(double));
    if (!yn || !ypn) { free(yn); free(ypn); return 0.0; }
    for (int i = 0; i < n; i++) {
        yn[i]  = y[i]  + eta[i];
        ypn[i] = yp[i] + etap[i];
    }
    double J0 = evaluate_functional(F, x, y, yp, n);
    double J1 = evaluate_functional(F, x, yn, ypn, n);
    free(yn); free(ypn);
    return J1 - J0;
}
void functional_energy_decomposition(LagrangianDensity F,
    const double *x, const double *y, const double *yp, int n,
    double *kinetic, double *potential)
{
    double T = 0.0, V = 0.0;
    if (!F || !x || !y || !yp || n < 2) goto done;
    double dx = x[1] - x[0], eps = 1e-6;
    for (int i = 0; i < n; i++) {
        double w = (i == 0 || i == n - 1) ? 0.5 : 1.0;
        double dF = (F(x[i],y[i],yp[i]+eps) - F(x[i],y[i],yp[i]-eps)) / (2.0*eps);
        T += w * dF * yp[i] * dx * 0.5;
        V += w * (-F(x[i], y[i], 0.0)) * dx;
    }
done:
    if (kinetic)   *kinetic   = T;
    if (potential) *potential = V;
}

/* ---- L3: Gateaux and Frechet Derivatives ---- */
double gateaux_derivative(LagrangianDensity F,
    const double *x, const double *y, const double *yp,
    const double *eta, const double *etap, int n, double epsilon)
{
    if (!F || !x || !y || !yp || !eta || !etap || n < 2) return 0.0;
    if (epsilon <= 0.0) epsilon = 1e-6;
    double *buf = (double*)malloc((size_t)n * 4 * sizeof(double));
    if (!buf) return 0.0;
    double *ypa  = buf,      *yma  = buf + n;
    double *yppa = buf + 2*n, *ypma = buf + 3*n;
    for (int i = 0; i < n; i++) {
        ypa[i]  = y[i]  + epsilon * eta[i];  yppa[i] = yp[i] + epsilon * etap[i];
        yma[i]  = y[i]  - epsilon * eta[i];  ypma[i] = yp[i] - epsilon * etap[i];
    }
    double Jp = evaluate_functional(F, x, ypa, yppa, n);
    double Jm = evaluate_functional(F, x, yma, ypma, n);
    free(buf);
    return (Jp - Jm) / (2.0 * epsilon);
}
void functional_derivative(LagrangianDensity dF_dy, LagrangianDensity dF_dyp,
    const double *x, const double *y, const double *yp, int n, double *result)
{
    if (!dF_dy || !dF_dyp || !x || !y || !yp || !result || n < 3) {
        if (result && n > 0)
            for (int i = 0; i < n; i++) result[i] = 0.0;
        return;
    }
    double dx = x[1] - x[0];
    for (int i = 0; i < n; i++) {
        double t1 = dF_dy(x[i], y[i], yp[i]), t2;
        if (i == 0)
            t2 = (dF_dyp(x[i+1],y[i+1],yp[i+1]) - dF_dyp(x[i],y[i],yp[i])) / dx;
        else if (i == n - 1)
            t2 = (dF_dyp(x[i],y[i],yp[i]) - dF_dyp(x[i-1],y[i-1],yp[i-1])) / dx;
        else
            t2 = (dF_dyp(x[i+1],y[i+1],yp[i+1]) - dF_dyp(x[i-1],y[i-1],yp[i-1])) / (2.0*dx);
        result[i] = t1 - t2;
    }
}
double functional_derivative_L2_norm(LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp, const double *x, const double *y,
    const double *yp, int n)
{
    if (!x || !y || !yp || n < 2) return 0.0;
    double *res = (double*)malloc((size_t)n * sizeof(double));
    if (!res) return 0.0;
    functional_derivative(dF_dy, dF_dyp, x, y, yp, n, res);
    double s = 0.0, dx = x[1] - x[0];
    for (int i = 0; i < n; i++) {
        double w = (i == 0 || i == n - 1) ? 0.5 : 1.0;
        s += w * res[i] * res[i] * dx;
    }
    free(res);
    return sqrt(s);
}
double frechet_derivative_test(LagrangianDensity F,
    const double *x, const double *y, const double *yp,
    const double *eta, const double *etap, int n)
{
    if (!F || !x || !y || !yp || !eta || !etap || n < 2) return 0.0;
    double Jy = evaluate_functional(F, x, y, yp, n);
    double *buf = (double*)malloc((size_t)n * 2 * sizeof(double));
    if (!buf) return 0.0;
    double *ya = buf, *yb = buf + n;
    for (int i = 0; i < n; i++) { ya[i] = y[i] + eta[i]; yb[i] = yp[i] + etap[i]; }
    double Jp = evaluate_functional(F, x, ya, yb, n);
    double dJ = gateaux_derivative(F, x, y, yp, eta, etap, n, 1e-6);
    free(buf);
    double R = fabs(Jp - Jy - dJ);
    double ne = functional_norm(eta, etap, x[0], x[n-1], n);
    return (ne < 1e-15) ? 0.0 : R / ne;
}
void frechet_derivative_convergence(LagrangianDensity F,
    const double *x, const double *y, const double *yp,
    const double *eta, const double *etap, int n, int n_steps, double *rates)
{
    if (!F || !x || !y || !yp || !eta || !etap || !rates || n < 2 || n_steps < 2) return;
    double *rv = (double*)malloc((size_t)n_steps * sizeof(double));
    double *es = (double*)malloc((size_t)n * 2 * sizeof(double));
    if (!rv || !es) { free(rv); free(es); return; }
    double *ep = es + n;
    for (int k = 0; k < n_steps; k++) {
        double sc = 1.0 / (double)(1 << k);
        for (int i = 0; i < n; i++) { es[i] = eta[i]*sc; ep[i] = etap[i]*sc; }
        rv[k] = frechet_derivative_test(F, x, y, yp, es, ep, n);
    }
    for (int k = 0; k < n_steps - 1; k++)
        rates[k] = (rv[k+1] < 1e-15) ? 0.0 : rv[k] / rv[k+1];
    free(rv); free(es);
}

/* ---- L4: Beltrami Identity ---- */
double beltrami_identity(LagrangianDensity F, LagrangianDensity dF_dyp,
    double x, double y, double yp)
{
    if (!F || !dF_dyp) return 0.0;
    return F(x, y, yp) - yp * dF_dyp(x, y, yp);
}
double verify_beltrami(LagrangianDensity F, LagrangianDensity dF_dyp,
    const double *x, const double *y, const double *yp, int n, double *values_out)
{
    if (!F || !dF_dyp || !x || !y || !yp || n < 2) return 0.0;
    double mn = DBL_MAX, mx = -DBL_MAX;
    for (int i = 0; i < n; i++) {
        double b = beltrami_identity(F, dF_dyp, x[i], y[i], yp[i]);
        if (values_out) values_out[i] = b;
        if (b < mn) mn = b;
        if (b > mx) mx = b;
    }
    return mx - mn;
}
double cyclic_coordinate_integral(LagrangianDensity dF_dyp,
    double x, double y, double yp)
{ if (!dF_dyp) return 0.0; return dF_dyp(x, y, yp); }

/* ---- L4: Weak/Strong Extrema ---- */
int weak_extremum_type(LagrangianDensity F, const double *x,
    const double *y, const double *yp, int n, double delta, int n_pert)
{
    if (!F || !x || !y || !yp || n < 2 || n_pert < 1) return 0;
    double J0 = evaluate_functional(F, x, y, yp, n);
    int ca = 0, cb = 0;
    double *buf = (double*)malloc((size_t)n * 4 * sizeof(double));
    if (!buf) return 0;
    double *eta = buf, *etap = buf + n, *ya = etap + n, *ypp = ya + n;
    double L = x[n-1] - x[0];
    for (int k = 0; k < n_pert; k++) {
        double ph = 2.0 * M_PI * (double)k / (double)n_pert;
        for (int i = 0; i < n; i++) {
            double xi = (x[i] - x[0]) / L;
            eta[i]  = delta * sin(ph + 3.0*xi) * xi * (1.0-xi);
            etap[i] = delta * 3.0 * cos(ph + 3.0*xi) * xi * (1.0-xi) / L;
            ya[i]  = y[i]  + eta[i];
            ypp[i] = yp[i] + etap[i];
        }
        double Jp = evaluate_functional(F, x, ya, ypp, n);
        if (Jp >= J0 - 1e-12) ca++;
        if (Jp <= J0 + 1e-12) cb++;
    }
    free(buf);
    if (ca >= n_pert) return 1;
    if (cb >= n_pert) return -1;
    return 0;
}
void compare_functionals(LagrangianDensity F, const double *x,
    const double *y1, const double *yp1, const double *y2, const double *yp2,
    int n, double *J1, double *J2, double *diff, double *rel_diff)
{
    double v1 = evaluate_functional(F, x, y1, yp1, n);
    double v2 = evaluate_functional(F, x, y2, yp2, n);
    double d = v2 - v1, rd = d / (fabs(v1) > DBL_MIN ? fabs(v1) : DBL_MIN);
    if (J1) *J1 = v1;
    if (J2) *J2 = v2;
    if (diff) *diff = d;
    if (rel_diff) *rel_diff = rd;
}

/* ---- L4: First and Second Variations ---- */
double first_variation_formula(LagrangianDensity dF_dy, LagrangianDensity dF_dyp,
    const double *x, const double *y, const double *yp,
    const double *eta, const double *etap, int n)
{
    if (!dF_dy || !dF_dyp || !x || !y || !yp || !eta || !etap || n < 2) return 0.0;
    double *ig = (double*)malloc((size_t)n * sizeof(double));
    if (!ig) return 0.0;
    for (int i = 0; i < n; i++)
        ig[i] = dF_dy(x[i],y[i],yp[i])*eta[i] + dF_dyp(x[i],y[i],yp[i])*etap[i];
    double r = trap_int(ig, x, n);
    free(ig);
    return r;
}
double second_variation_formula(LagrangianDensity d2F_dyp2,
    LagrangianDensity d2F_dy2, const double *x,
    const double *eta, const double *etap, int n)
{
    if (!d2F_dyp2 || !d2F_dy2 || !x || !eta || !etap || n < 2) return 0.0;
    double *ig = (double*)malloc((size_t)n * sizeof(double));
    if (!ig) return 0.0;
    for (int i = 0; i < n; i++) {
        double P = d2F_dyp2(x[i], 0.0, 0.0), Q = d2F_dy2(x[i], 0.0, 0.0);
        ig[i] = P * etap[i] * etap[i] + Q * eta[i] * eta[i];
    }
    double r = trap_int(ig, x, n);
    free(ig);
    return r;
}
double second_variation_jacobi_form(LagrangianDensity d2F_dyp2,
    LagrangianDensity d2F_dy2, LagrangianDensity d2F_dydyp,
    const double *x, const double *y, const double *yp,
    const double *eta, const double *etap, int n)
{
    if (!d2F_dyp2 || !d2F_dy2 || !d2F_dydyp || !x || !eta || !etap || n < 3) return 0.0;
    double dx = x[1] - x[0];
    double *Je = (double*)malloc((size_t)n * 2 * sizeof(double));
    if (!Je) return 0.0;
    double *ig = Je + n;
    for (int i = 1; i < n - 1; i++) {
        double Pip = d2F_dyp2(x[i+1], y[i+1], yp[i+1]);
        double Pim = d2F_dyp2(x[i-1], y[i-1], yp[i-1]);
        double Qi  = d2F_dy2(x[i],    y[i],   yp[i]);
        double mxd = (d2F_dydyp(x[i+1],y[i+1],yp[i+1]) -
                      d2F_dydyp(x[i-1],y[i-1],yp[i-1])) / (2.0*dx);
        double Qeff = Qi - mxd;
        double lap  = (Pip*(eta[i+1]-eta[i])/dx - Pim*(eta[i]-eta[i-1])/dx) / dx;
        (void)d2F_dyp2; /* Pi not needed separately, Pip/Pim used */
        Je[i] = -lap + Qeff * eta[i];
    }
    Je[0] = 0.0; Je[n-1] = 0.0;
    for (int i = 0; i < n; i++) ig[i] = eta[i] * Je[i];
    double r = trap_int(ig, x, n);
    free(Je);
    return r;
}

/* ---- L3: Convexity ---- */
bool is_convex_functional(LagrangianDensity F, const double *x,
    const double *y1, const double *yp1, const double *y2, const double *yp2,
    int n, int n_theta)
{
    if (!F || !x || !y1 || !yp1 || !y2 || !yp2 || n < 2 || n_theta < 2) return false;
    double J1 = evaluate_functional(F, x, y1, yp1, n);
    double J2 = evaluate_functional(F, x, y2, yp2, n);
    double *buf = (double*)malloc((size_t)n * 2 * sizeof(double));
    if (!buf) return false;
    double *ym = buf, *ypm = buf + n;
    bool cvx = true;
    for (int t = 0; t < n_theta; t++) {
        double th = (double)(t+1) / (double)(n_theta+1);
        for (int i = 0; i < n; i++) {
            ym[i]  = th * y1[i]  + (1.0-th) * y2[i];
            ypm[i] = th * yp1[i] + (1.0-th) * yp2[i];
        }
        double Jm = evaluate_functional(F, x, ym, ypm, n);
        double Jl = th * J1 + (1.0-th) * J2;
        if (Jm > Jl + 1e-12) { cvx = false; break; }
    }
    free(buf);
    return cvx;
}
bool is_quadratic_functional(LagrangianDensity F,
    const double *x, const double *y, const double *yp, int n)
{
    if (!F || !x || !y || !yp || n < 2) return false;
    double sc[5] = {0.5, 0.75, 1.0, 1.25, 1.5};
    bool q = true;
    for (int i = 0; i < n && q; i++) {
        double fv[5];
        for (int s = 0; s < 5; s++) fv[s] = F(x[i], sc[s]*y[i], sc[s]*yp[i]);
        double d1[4], d2[3], d3[2];
        for (int j = 0; j < 4; j++) d1[j] = fv[j+1] - fv[j];
        for (int j = 0; j < 3; j++) d2[j] = d1[j+1] - d1[j];
        for (int j = 0; j < 2; j++) d3[j] = d2[j+1] - d2[j];
        if (fabs(d3[0]) > 1e-8 || fabs(d3[1]) > 1e-8) q = false;
    }
    return q;
}

/* ---- L5: Utility Functions ---- */
void generate_uniform_nodes(double a, double b, int n, double *x) {
    if (!x || n < 2) return;
    double dx = (b - a) / (double)(n - 1);
    for (int i = 0; i < n; i++) x[i] = a + (double)i * dx;
}
void generate_chebyshev_nodes(double a, double b, int n, double *x) {
    if (!x || n < 2) return;
    double mid = 0.5 * (a + b), hf = 0.5 * (b - a);
    for (int i = 0; i < n; i++) {
        double th = M_PI * (double)i / (double)(n - 1);
        x[i] = mid + hf * cos(th);
    }
}
void evaluate_function_samples(double (*yf)(double), double (*ypf)(double),
    const double *x, int n, double *y, double *yp) {
    if (!x || !y || !yp || n < 1) return;
    for (int i = 0; i < n; i++) {
        y[i]  = yf  ? yf(x[i])  : 0.0;
        yp[i] = ypf ? ypf(x[i]) : 0.0;
    }
}
void initialize_functional_context(FunctionalContext *ctx,
    LagrangianDensity F, LagrangianDensity dF_dy,
    LagrangianDensity dF_dyp, LagrangianDensity d2F_dyp2,
    double a, double b) {
    if (!ctx) return;
    ctx->F = F; ctx->dF_dy = dF_dy; ctx->dF_dyp = dF_dyp;
    ctx->d2F_dyp2 = d2F_dyp2; ctx->a = a; ctx->b = b;
}
void compute_derivative_fd(double (*yf)(double), const double *x, int n, double *yp) {
    if (!yf || !x || !yp || n < 5) {
        if (yp && n >= 2) {
            double h = x[1] - x[0];
            for (int i = 0; i < n; i++) yp[i] = (yf(x[i]+h) - yf(x[i]-h)) / (2.0*h);
        }
        return;
    }
    double h = x[1] - x[0];
    for (int i = 2; i < n - 2; i++) {
        double fm2 = yf(x[i-2]), fm1 = yf(x[i-1]);
        double fp1 = yf(x[i+1]), fp2 = yf(x[i+2]);
        yp[i] = (-fp2 + 8.0*fp1 - 8.0*fm1 + fm2) / (12.0 * h);
    }
    double f0 = yf(x[0]), f1 = yf(x[1]), f2 = yf(x[2]);
    double f3 = yf(x[3]), f4 = yf(x[4]);
    yp[0] = (-25.0*f0 + 48.0*f1 - 36.0*f2 + 16.0*f3 - 3.0*f4) / (12.0 * h);
    yp[1] = (-3.0*f0  - 10.0*f1 + 18.0*f2 -  6.0*f3 +     f4) / (12.0 * h);
    double fn5 = yf(x[n-5]), fn4 = yf(x[n-4]), fn3 = yf(x[n-3]);
    double fn2 = yf(x[n-2]), fn1 = yf(x[n-1]);
    yp[n-2] = (-fn5 + 6.0*fn4 - 18.0*fn3 + 10.0*fn2 + 3.0*fn1) / (12.0 * h);
    yp[n-1] = (3.0*fn5 - 16.0*fn4 + 36.0*fn3 - 48.0*fn2 + 25.0*fn1) / (12.0 * h);
}
double linear_interpolate(const double *x, const double *y, int n, double xq) {
    if (!x || !y || n < 1) return 0.0;
    if (n == 1) return y[0];
    int lo = 0, hi = n - 1;
    if (xq <= x[lo]) return y[lo];
    if (xq >= x[hi]) return y[hi];
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (x[mid] <= xq) lo = mid; else hi = mid;
    }
    double t = (xq - x[lo]) / (x[hi] - x[lo]);
    return y[lo] + t * (y[hi] - y[lo]);
}
