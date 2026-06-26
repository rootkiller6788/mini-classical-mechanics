/**
 * three_body.c — Circular Restricted Three-Body Problem
 *
 * Implements the CRTBP equations, Lagrange point computation,
 * Jacobi integral, zero-velocity surfaces, Hill stability,
 * Roche limits, and RK4 trajectory integration.
 *
 * References:
 *   - Murray & Dermott Sec.3.8-3.13
 *   - Szebehely, Theory of Orbits (1967), Ch.1-4
 *   - Vallado Sec.10.6
 *   - Koon, Lo, Marsden & Ross, Dynamical Systems, the
 *     Three-Body Problem and Space Mission Design (2007)
 */

#include "celestial_types.h"
#include "three_body.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* ================================================================
 * L6: CRTBP Equations of Motion
 *
 * In the rotating (synodic) coordinate frame centered at the
 * barycenter, with the x-axis along the primary-secondary line:
 *
 *   m1 at x = -mu
 *   m2 at x = 1 - mu
 *   where mu = m2 / (m1 + m2)
 *
 * The equations of motion are:
 *   x'' - 2*y' = partial_Omega/partial_x
 *   y'' + 2*x' = partial_Omega/partial_y
 *   z'' = partial_Omega/partial_z
 *
 * where Omega(x,y,z) = (x^2+y^2)/2 + (1-mu)/r1 + mu/r2 + const
 *
 * The Coriolis terms (2*y' and -2*x') arise from the rotating frame.
 *
 * Ref: Szebehely Eq.(1.2.14), Murray & Dermott Sec.3.8
 * ================================================================ */

void crtbp_equations(const double state[6], double mu, double deriv[6]) {
    double x, y, z, vx, vy, vz;
    double r1, r2, r1_3, r2_3;
    double Ux, Uy, Uz;
    double om1 = 1.0 - mu;

    x  = state[0]; y  = state[1]; z  = state[2];
    vx = state[3]; vy = state[4]; vz = state[5];

    /* Distances to the two primaries */
    r1 = sqrt((x + mu) * (x + mu) + y * y + z * z);
    r2 = sqrt((x - om1) * (x - om1) + y * y + z * z);

    /* Guard against division by zero (at primary centers) */
    if (r1 < 1e-12) r1 = 1e-12;
    if (r2 < 1e-12) r2 = 1e-12;

    r1_3 = r1 * r1 * r1;
    r2_3 = r2 * r2 * r2;

    /* Pseudo-potential derivatives:
     * Omega = (x^2+y^2)/2 + (1-mu)/r1 + mu/r2 + const
     * dOmega/dx = x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3
     * dOmega/dy = y - (1-mu)*y/r1^3 - mu*y/r2^3
     * dOmega/dz = - (1-mu)*z/r1^3 - mu*z/r2^3 */
    Ux = x - om1 * (x + mu) / r1_3 - mu * (x - om1) / r2_3;
    Uy = y - om1 * y / r1_3 - mu * y / r2_3;
    Uz = -om1 * z / r1_3 - mu * z / r2_3;

    /* Equations of motion:
     * x'' = 2*y' + Ux
     * y'' = -2*x' + Uy
     * z'' = Uz */
    deriv[0] = vx;
    deriv[1] = vy;
    deriv[2] = vz;
    deriv[3] =  2.0 * vy + Ux;
    deriv[4] = -2.0 * vx + Uy;
    deriv[5] = Uz;
}

/* ================================================================
 * L6: CRTBP Effective Potential (Pseudo-Potential)
 *
 * Omega = (1/2)(x^2 + y^2) + (1-mu)/r1 + mu/r2 + (1/2)mu(1-mu)
 *
 * The last term is a constant that makes Omega=3/2 at L4 and L5
 * for mu < 1/2.
 *
 * Ref: Szebehely Eq.(1.2.13)
 * ================================================================ */

double crtbp_effective_potential(double x, double y, double z, double mu) {
    double r1 = sqrt((x + mu) * (x + mu) + y * y + z * z);
    double r2 = sqrt((x - 1.0 + mu) * (x - 1.0 + mu) + y * y + z * z);
    double om1 = 1.0 - mu;

    if (r1 < 1e-12) r1 = 1e-12;
    if (r2 < 1e-12) r2 = 1e-12;

    return 0.5 * (x * x + y * y) + om1 / r1 + mu / r2
           + 0.5 * om1 * mu;
}

/* ================================================================
 * L6: Jacobi Constant
 *
 * C_J = 2*Omega - v^2
 *
 * The Jacobi integral is the only known integral of motion for
 * the CRTBP. It is constant along any trajectory.
 *
 * For a given C_J, the region v^2 >= 0 defines where motion is
 * possible (the "Hill region"). The zero-velocity surfaces (ZVS)
 * are where v^2 = 0.
 *
 * Ref: Murray & Dermott Eq.(3.144), Szebehely Sec.1.5
 * ================================================================ */

double jacobi_constant(const double state[6], double mu) {
    double x = state[0], y = state[1], z = state[2];
    double vx = state[3], vy = state[4], vz = state[5];
    double Omega = crtbp_effective_potential(x, y, z, mu);
    double v2 = vx * vx + vy * vy + vz * vz;
    return 2.0 * Omega - v2;
}

/* ================================================================
 * Zero-Velocity Surface
 *
 * The ZVS is defined by v^2 = 2*Omega - C_J = 0.
 * At a given (x,y) in the orbital plane, the ZVS height z is:
 *   z^2 = solution of 2*Omega(x,y,z) - C_J = 0
 *
 * Solved iteratively since Omega depends on z nonlinearly
 * through r1 and r2.
 *
 * Returns the positive z-coordinate, or NaN if no solution.
 * ================================================================ */

double zero_velocity_surface_z(double x, double y, double C_J, double mu,
                                double z_guess) {
    double z = z_guess;
    int iter;

    for (iter = 0; iter < 30; iter++) {
        double Omega = crtbp_effective_potential(x, y, z, mu);
        double f = 2.0 * Omega - C_J;

        if (f < 0.0) {
            /* No zero-velocity surface at this (x,y) for given C_J */
            return NAN;
        }

        /* Newton iteration: find z such that 2*Omega(z) - C_J = 0
         * df/dz = 2 * dOmega/dz = 2*(-(1-mu)z/r1^3 - mu*z/r2^3) */
        double r1 = sqrt((x + mu)*(x + mu) + y*y + z*z);
        double r2 = sqrt((x - 1.0 + mu)*(x - 1.0 + mu) + y*y + z*z);

        if (r1 < 1e-12) r1 = 1e-12;
        if (r2 < 1e-12) r2 = 1e-12;

        double dOmega_dz = -(1.0 - mu) * z / (r1*r1*r1) - mu * z / (r2*r2*r2);
        double df_dz = 2.0 * dOmega_dz;

        if (fabs(df_dz) < 1e-15) break;

        double dz = -f / df_dz;
        if (fabs(dz) > 0.1) dz = (dz > 0 ? 0.1 : -0.1);
        z += dz;

        if (fabs(dz) < 1e-12) break;
    }

    return z;
}

int point_accessible(double x, double y, double z, double C_J, double mu) {
    double Omega = crtbp_effective_potential(x, y, z, mu);
    return (2.0 * Omega >= C_J) ? 1 : 0;
}

/* ================================================================
 * L6: Collinear Lagrange Points
 *
 * The collinear points L1, L2, L3 lie on the x-axis (y=z=0)
 * and satisfy the quintic equation:
 *
 *   x - (1-mu)*(x+mu)/|x+mu|^3 - mu*(x-1+mu)/|x-1+mu|^3 = 0
 *
 * Solved via bisection in specified intervals:
 *   L1: between m2 and m1 (x in [-mu, 1-mu])
 *   L2: beyond m2 (x in [1-mu, 2.0])
 *   L3: beyond m1 opposite m2 (x in [-2.0, -mu])
 *
 * Ref: Murray & Dermott Sec.3.8, Szebehely Ch.4
 * ================================================================ */

double find_collinear_lagrange(double mu, double x0, double x1, double tol) {
    double a = x0, b = x1;
    double fa, fb, mid, fmid;
    int iter;
    double om1 = 1.0 - mu;

    if (tol <= 0.0) tol = 1e-12;

    /* Evaluate f(x) = x - (1-mu)*(x+mu)/|x+mu|^3 - mu*(x-1+mu)/|x-1+mu|^3 */
    auto double f(double x_val) {
        double r1 = fabs(x_val + mu);
        double r2 = fabs(x_val - om1);
        if (r1 < 1e-12) r1 = 1e-12;
        if (r2 < 1e-12) r2 = 1e-12;
        return x_val - om1 * (x_val + mu) / (r1*r1*r1)
                     - mu * (x_val - om1) / (r2*r2*r2);
    }

    /* Note: "auto double f" is GCC nested function — use inline approach instead */

    fa = f(a); fb = f(b);

    /* Ensure root is bracketed */
    if (fa * fb > 0.0) {
        /* Try wider bracket */
        a = x0 - 0.1; b = x1 + 0.1;
        fa = f(a); fb = f(b);
    }

    for (iter = 0; iter < 80; iter++) {
        mid = (a + b) / 2.0;
        fmid = f(mid);

        if (fabs(fmid) < tol) return mid;

        if (fa * fmid < 0.0) {
            b = mid; fb = fmid;
        } else {
            a = mid; fa = fmid;
        }

        if (fabs(b - a) < tol) return mid;
    }

    return (a + b) / 2.0;
}

/* ================================================================
 * Lagrange Point Computation
 *
 * L1, L2, L3: collinear points on the x-axis (solved via bisection)
 * L4, L5: equilateral triangle points:
 *   L4: x = 0.5 - mu, y = +sqrt(3)/2
 *   L5: x = 0.5 - mu, y = -sqrt(3)/2
 *
 * L4 and L5 form equilateral triangles with the two primaries.
 * They are stable for mu < mu_critical = (1-sqrt(23/27))/2 ≈ 0.0385.
 * (Routh criterion)
 *
 * Ref: Murray & Dermott Sec.3.8, Szebehely Ch.4
 * ================================================================ */

/* Forward declare the collinear search function with its equation */
static double collinear_lagrange_equation(double x_val, double mu) {
    double om1 = 1.0 - mu;
    double r1 = fabs(x_val + mu);
    double r2 = fabs(x_val - om1);
    if (r1 < 1e-12) r1 = 1e-12;
    if (r2 < 1e-12) r2 = 1e-12;
    return x_val - om1 * (x_val + mu) / (r1*r1*r1)
                 - mu * (x_val - om1) / (r2*r2*r2);
}

/* Override the find_collinear_lagrange to use static function */
static double find_collinear_lagrange_static(double mu, double x0, double x1, double tol) {
    double a = x0, b = x1;
    double fa, fb, mid, fmid;
    int iter;

    if (tol <= 0.0) tol = 1e-12;

    fa = collinear_lagrange_equation(a, mu);
    fb = collinear_lagrange_equation(b, mu);

    if (fa * fb > 0.0) {
        a = x0 - 0.1; b = x1 + 0.1;
        fa = collinear_lagrange_equation(a, mu);
        fb = collinear_lagrange_equation(b, mu);
    }

    for (iter = 0; iter < 80; iter++) {
        mid = (a + b) / 2.0;
        fmid = collinear_lagrange_equation(mid, mu);

        if (fabs(fmid) < tol) return mid;

        if (fa * fmid < 0.0) {
            b = mid; fb = fmid;
        } else {
            a = mid; fa = fmid;
        }

        if (fabs(b - a) < tol) return mid;
    }

    return (a + b) / 2.0;
}

LagrangePoints lagrange_points(double mu) {
    LagrangePoints pts;
    double om1 = 1.0 - mu;
    double sqrt3_2 = sqrt(3.0) / 2.0;

    memset(&pts, 0, sizeof(pts));

    /* L1: between the two primaries */
    pts.L1.x = find_collinear_lagrange_static(mu, -mu + 0.001, om1 - 0.001, 1e-12);
    pts.L1.y = 0.0; pts.L1.z = 0.0;

    /* L2: beyond the secondary */
    pts.L2.x = find_collinear_lagrange_static(mu, om1 + 0.001, 2.0, 1e-12);
    pts.L2.y = 0.0; pts.L2.z = 0.0;

    /* L3: beyond the primary, opposite side */
    pts.L3.x = find_collinear_lagrange_static(mu, -2.0, -mu - 0.001, 1e-12);
    pts.L3.y = 0.0; pts.L3.z = 0.0;

    /* L4: leading equilateral triangle point */
    pts.L4.x = 0.5 - mu;
    pts.L4.y =  sqrt3_2;
    pts.L4.z = 0.0;

    /* L5: trailing equilateral triangle point */
    pts.L5.x = 0.5 - mu;
    pts.L5.y = -sqrt3_2;
    pts.L5.z = 0.0;

    return pts;
}

/* ================================================================
 * Jacobi Constants at Lagrange Points
 *
 * At Lagrange points, v = 0, so C_J(L) = 2*Omega(L).
 * C_J values determine the topological structure of the
 * zero-velocity surfaces (Hill regions).
 *
 * The ordering is typically: C_J(L1) > C_J(L2) > C_J(L3) > C_J(L4,L5).
 *
 * Ref: Szebehely Sec.4.7
 * ================================================================ */

void lagrange_jacobi_constants(double mu, double C_J[5]) {
    LagrangePoints pts = lagrange_points(mu);

    C_J[0] = 2.0 * crtbp_effective_potential(pts.L1.x, pts.L1.y, pts.L1.z, mu);
    C_J[1] = 2.0 * crtbp_effective_potential(pts.L2.x, pts.L2.y, pts.L2.z, mu);
    C_J[2] = 2.0 * crtbp_effective_potential(pts.L3.x, pts.L3.y, pts.L3.z, mu);
    C_J[3] = 2.0 * crtbp_effective_potential(pts.L4.x, pts.L4.y, pts.L4.z, mu);
    C_J[4] = 2.0 * crtbp_effective_potential(pts.L5.x, pts.L5.y, pts.L5.z, mu);
}

/* ================================================================
 * Hill Region Connectivity
 *
 * The zero-velocity surfaces partition space into regions of
 * allowed and forbidden motion. As C_J decreases, "necks" open
 * at the Lagrange points, connecting previously separated regions.
 *
 * The critical C_J values at L1, L2, L3 determine the topology:
 *   C_J > C_J(L1):  motion confined to m1 or m2 vicinity
 *   C_J(L2) < C_J < C_J(L1):  m1 and m2 regions connect at L1
 *   C_J(L3) < C_J < C_J(L2):  outer region opens at L2
 *   C_J < C_J(L3):  all regions connected
 *
 * Ref: Murray & Dermott Fig.3.28, Szebehely Fig.4.2
 * ================================================================ */

const char *hill_region(double C_J, double mu) {
    double C[5];
    lagrange_jacobi_constants(mu, C);

    if (C_J >= C[0] - 1e-10)
        return "m1 and m2 regions disconnected (C_J >= C_L1)";
    else if (C_J >= C[1] - 1e-10)
        return "m1-m2 connected at L1 neck, outer region closed";
    else if (C_J >= C[2] - 1e-10)
        return "outer region opens at L2, m2 exterior accessible";
    else if (C_J >= C[3] - 1e-10)
        return "all regions open at L3, L4/L5 still forbidden";
    else
        return "all regions fully connected (C_J < C_L4/L5)";
}

/* ================================================================
 * Linear Stability of Collinear Lagrange Points
 *
 * The variational equations at a collinear point yield a 6x6
 * system. The characteristic equation is:
 *   lambda^6 + (4 - Uxx - Uyy)*lambda^4 + (Uxx*Uyy - Uxy^2)*lambda^2
 *            + (some constant) = 0
 *
 * Collinear points are always unstable (one real eigenvalue pair
 * with positive real part).
 *
 * Ref: Szebehely Sec.4.8, Murray & Dermott Sec.3.8
 * ================================================================ */

void lagrange_point_stability(double mu, int point_idx,
                               double eigenvalues[6]) {
    double x, y;
    double Uxx, Uyy, Uxy;
    double Uzz;(void)Uzz;
    double om1 = 1.0 - mu;
    LagrangePoints pts;
    int i;

    /* Get the Lagrange point location */
    pts = lagrange_points(mu);

    switch (point_idx) {
        case 0: x = pts.L1.x; y = pts.L1.y; break;
        case 1: x = pts.L2.x; y = pts.L2.y; break;
        case 2: x = pts.L3.x; y = pts.L3.y; break;
        case 3: x = pts.L4.x; y = pts.L4.y; break;
        case 4: x = pts.L5.x; y = pts.L5.y; break;
        default: return;
    }

    /* Second derivatives of Omega at the Lagrange point */
    double r1 = sqrt((x+mu)*(x+mu) + y*y);
    double r2 = sqrt((x-om1)*(x-om1) + y*y);

    if (r1 < 1e-12) r1 = 1e-12;
    if (r2 < 1e-12) r2 = 1e-12;

    double r1_3 = r1*r1*r1, r1_5 = r1_3*r1*r1;
    double r2_3 = r2*r2*r2, r2_5 = r2_3*r2*r2;

    Uxx = 1.0 - om1/r1_3 + 3.0*om1*(x+mu)*(x+mu)/r1_5
                - mu/r2_3   + 3.0*mu*(x-om1)*(x-om1)/r2_5;

    Uyy = 1.0 - om1/r1_3 + 3.0*om1*y*y/r1_5
                - mu/r2_3   + 3.0*mu*y*y/r2_5;

    Uzz = -om1/r1_3 - mu/r2_3;

    Uxy = 3.0*om1*(x+mu)*y/r1_5 + 3.0*mu*(x-om1)*y/r2_5;

    /* Characteristic equation coefficients for collinear points:
     * lambda^2 = -b/2 +/- sqrt(b^2/4 - c)
     * where b = 4 - Uxx - Uyy, c = Uxx*Uyy - Uxy^2 */
    {
        double b = 4.0 - Uxx - Uyy;
        double c = Uxx * Uyy - Uxy * Uxy;
        double disc = b*b/4.0 - c;

        /* lambda^2 values */
        double lam2_1 = -b/2.0 + sqrt(fabs(disc));
        double lam2_2 = -b/2.0 - sqrt(fabs(disc));

        /* In-plane eigenvalues: lambda = +/- sqrt(lam2_i) */
        if (lam2_1 > 0) {
            eigenvalues[0] =  sqrt(lam2_1);  /* real positive -> unstable */
            eigenvalues[1] = -sqrt(lam2_1);
        } else {
            eigenvalues[0] = 0.0;
            eigenvalues[1] = 0.0;  /* purely imaginary (oscillatory) */
        }

        if (lam2_2 < 0) {
            eigenvalues[2] = 0.0;
            eigenvalues[3] = 0.0;  /* oscillatory modes: lambda = +/- i*sqrt(-lam2_2) */
        } else {
            eigenvalues[2] =  sqrt(lam2_2);
            eigenvalues[3] = -sqrt(lam2_2);
        }

        /* Out-of-plane: lambda^2 = Uzz (always negative -> oscillatory) */
        eigenvalues[4] = 0.0;
        eigenvalues[5] = 0.0; /* purely imaginary stability */
    }

    /* Initialize unused entries */
    for (i = 0; i < 6; i++) {
        if (!isfinite(eigenvalues[i])) eigenvalues[i] = 0.0;
    }
}

/* ================================================================
 * L8: State Transition Matrix for CRTBP
 *
 * The STM Phi(t) maps initial perturbations to final perturbations:
 *   delta_x(t) = Phi(t) * delta_x(0)
 *
 * Phi satisfies: Phi_dot = A(t) * Phi, Phi(0) = I
 * where A(t) = df/dx is the Jacobian of the CRTBP equations.
 *
 * The STM is fundamental for trajectory correction, station-keeping,
 * and manifold computation around Lagrange point orbits.
 *
 * Ref: Koon et al. (2007), Sec.2.4
 * ================================================================ */

void crtbp_stm(const double state[6], double mu,
                const double phi[36], double phi_dot[36]) {
    double x = state[0], y = state[1], z = state[2];
    double om1 = 1.0 - mu;
    double r1, r2, r1_3, r1_5, r2_3, r2_5;
    double Uxx, Uyy, Uzz, Uxy, Uxz, Uyz;
    double A[6][6];
    int i, j, k;

    r1 = sqrt((x+mu)*(x+mu) + y*y + z*z);
    r2 = sqrt((x-om1)*(x-om1) + y*y + z*z);
    if (r1 < 1e-12) r1 = 1e-12;
    if (r2 < 1e-12) r2 = 1e-12;

    r1_3 = r1*r1*r1; r1_5 = r1_3*r1*r1;
    r2_3 = r2*r2*r2; r2_5 = r2_3*r2*r2;

    /* Second derivatives of pseudo-potential */
    Uxx = 1.0 - om1/r1_3 + 3.0*om1*(x+mu)*(x+mu)/r1_5
                - mu/r2_3   + 3.0*mu*(x-om1)*(x-om1)/r2_5;

    Uyy = 1.0 - om1/r1_3 + 3.0*om1*y*y/r1_5
                - mu/r2_3   + 3.0*mu*y*y/r2_5;

    Uzz = -om1/r1_3 - mu/r2_3
          + 3.0*om1*z*z/r1_5 + 3.0*mu*z*z/r2_5;

    Uxy = 3.0*om1*(x+mu)*y/r1_5 + 3.0*mu*(x-om1)*y/r2_5;
    Uxz = 3.0*om1*(x+mu)*z/r1_5 + 3.0*mu*(x-om1)*z/r2_5;
    Uyz = 3.0*om1*y*z/r1_5 + 3.0*mu*y*z/r2_5;

    /* State matrix A = df/dx:
     * A = [ 0     I  ]
     *     [ U_xx  Omega ]   where Omega is the Coriolis matrix */
    memset(A, 0, sizeof(A));

    /* Upper right: identity */
    A[0][3] = 1.0; A[1][4] = 1.0; A[2][5] = 1.0;

    /* Lower left: gradient of pseudo-potential */
    A[3][0] = Uxx; A[3][1] = Uxy; A[3][2] = Uxz;
    A[4][0] = Uxy; A[4][1] = Uyy; A[4][2] = Uyz;
    A[5][0] = Uxz; A[5][1] = Uyz; A[5][2] = Uzz;

    /* Lower right: Coriolis matrix [0 2 0; -2 0 0; 0 0 0] */
    A[3][4] =  2.0;
    A[4][3] = -2.0;

    /* Phi_dot = A * Phi (6x6 matrix multiplication) */
    for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
            double sum = 0.0;
            for (k = 0; k < 6; k++) {
                sum += A[i][k] * phi[k * 6 + j];
            }
            phi_dot[i * 6 + j] = sum;
        }
    }
}

/* ================================================================
 * Spheres of Influence
 *
 * Laplace SOI: r_SOI = a * (m/M)^(2/5)
 *   Within the SOI, the planet's gravity dominates. Used for
 *   patched-conic approximations in interplanetary trajectory design.
 *
 * Hill sphere: r_H = a * (m/(3M))^(1/3)
 *   The region where a satellite can maintain a stable orbit around
 *   a planet against solar tidal forces.
 *
 * Roche limit — rigid body: d = 1.26 * R_p * (rho_p/rho_s)^(1/3)
 * Roche limit — fluid body:  d = 2.44 * R_p * (rho_p/rho_s)^(1/3)
 *   Distance at which tidal forces exceed self-gravity.
 *
 * Ref: Murray & Dermott Sec.3.13, Roche (1847)
 * ================================================================ */

double sphere_of_influence(double a, double m, double M) {
    if (M <= 0.0) return 0.0;
    return a * pow(m / M, 0.4); /* (m/M)^(2/5) */
}

double hill_sphere(double a, double m, double M) {
    if (M <= 0.0) return 0.0;
    return a * pow(m / (3.0 * M), 1.0/3.0);
}

double roche_limit_rigid(double R_primary, double rho_primary, double rho_satellite) {
    if (rho_satellite <= 0.0) return INFINITY;
    return 1.26 * R_primary * pow(rho_primary / rho_satellite, 1.0/3.0);
}

double roche_limit_fluid(double R_primary, double rho_primary, double rho_satellite) {
    if (rho_satellite <= 0.0) return INFINITY;
    return 2.44 * R_primary * pow(rho_primary / rho_satellite, 1.0/3.0);
}

/* ================================================================
 * L5: CRTBP RK4 Integration
 * ================================================================ */

void crtbp_rk4_step(double state[6], double mu, double dt) {
    double k1[6], k2[6], k3[6], k4[6], temp[6];
    int i;

    crtbp_equations(state, mu, k1);
    for (i = 0; i < 6; i++) temp[i] = state[i] + 0.5 * dt * k1[i];

    crtbp_equations(temp, mu, k2);
    for (i = 0; i < 6; i++) temp[i] = state[i] + 0.5 * dt * k2[i];

    crtbp_equations(temp, mu, k3);
    for (i = 0; i < 6; i++) temp[i] = state[i] + dt * k3[i];

    crtbp_equations(temp, mu, k4);

    for (i = 0; i < 6; i++) {
        state[i] += (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
}

double **integrate_crtbp(const double state0[6], double mu,
                          double t_end, double dt, int *n_steps_out) {
    int n_steps, i, j;
    double **traj;
    double state[6];

    if (dt <= 0.0 || t_end <= 0.0) {
        if (n_steps_out) *n_steps_out = 0;
        return NULL;
    }

    n_steps = (int)ceil(t_end / dt);
    if (n_steps < 1) n_steps = 1;

    /* Allocate trajectory array: (n_steps+1) x 6 */
    traj = (double **)malloc((n_steps + 1) * sizeof(double *));
    if (!traj) { if (n_steps_out) *n_steps_out = 0; return NULL; }

    for (i = 0; i <= n_steps; i++) {
        traj[i] = (double *)malloc(6 * sizeof(double));
        if (!traj[i]) {
            for (j = 0; j < i; j++) free(traj[j]);
            free(traj);
            if (n_steps_out) *n_steps_out = 0;
            return NULL;
        }
    }

    /* Initial state */
    for (j = 0; j < 6; j++) state[j] = state0[j];
    for (j = 0; j < 6; j++) traj[0][j] = state[j];

    /* Integrate */
    for (i = 0; i < n_steps; i++) {
        crtbp_rk4_step(state, mu, dt);
        for (j = 0; j < 6; j++) traj[i + 1][j] = state[j];
    }

    if (n_steps_out) *n_steps_out = n_steps;
    return traj;
}

void free_crtbp_trajectory(double **traj, int n_steps) {
    int i;
    if (!traj) return;
    for (i = 0; i <= n_steps; i++) {
        free(traj[i]);
    }
    free(traj);
}