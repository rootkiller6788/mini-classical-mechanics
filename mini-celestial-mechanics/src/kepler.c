/**
 * kepler.c — Kepler Equation Solvers and Orbital Element Conversions
 *
 * Implements the transcendental Kepler equation and all related
 * computations: Newton-Raphson solvers for elliptic/hyperbolic/parabolic
 * orbits, anomaly conversions, and the fundamental elements↔state
 * transformations with full 3-1-3 Euler angle rotation.
 *
 * References:
 *   - Goldstein, Poole & Safko Sec.3.5-3.7
 *   - Murray & Dermott Sec.2.4-2.8
 *   - Vallado Sec.2.2-2.5
 *   - Danby, Fundamentals of Celestial Mechanics (1988), Ch.6
 *   - Battin, An Introduction to the Mathematics and Methods of
 *     Astrodynamics (1987), Sec.4.4
 */

#include "celestial_types.h"
#include "kepler.h"
#include <math.h>
#include <float.h>
#include <string.h>

#define KEPLER_DEFAULT_TOL 1e-12
#define KEPLER_MAX_ITER     50

/* ================================================================
 * L4: Elliptic Kepler Equation Solver — M = E - e sin(E)
 *
 * Newton-Raphson with Danby (1988) initial guess.
 * E_{n+1} = E_n - (E_n - e sin(E_n) - M) / (1 - e cos(E_n))
 * ================================================================ */

double solve_kepler(double M, double e, double tol) {
    double M_red, E, dE, f, fp;
    int iter;

    if (tol <= 0.0) tol = KEPLER_DEFAULT_TOL;

    M_red = fmod(M, TWOPI);
    if (M_red < 0.0) M_red += TWOPI;

    /* Circular orbit: trivial solution */
    if (e < 1e-15) return M_red;

    /* Danby optimal initial guess */
    {
        double sM = sin(M_red);
        double denom = 1.0 - sin(M_red + e) + sM;
        if (fabs(denom) < 1e-15) denom = 1e-15;
        E = M_red + e * sM / denom;
    }

    E = clamp_val(E, M_red - e, M_red + e);

    for (iter = 0; iter < KEPLER_MAX_ITER; iter++) {
        double sinE = sin(E);
        double cosE = cos(E);
        f = E - e * sinE - M_red;
        fp = 1.0 - e * cosE;

        if (fabs(fp) < 1e-30) {
            dE = -f / (fp + 0.5 * e * sinE * (f / (fp + 1e-15)));
        } else {
            dE = -f / fp;
        }

        if (fabs(dE) > 1.0) dE = (dE > 0 ? 1.0 : -1.0);
        E += dE;
        if (fabs(dE) < tol) break;
    }
    return wrap_2pi(E);
}

/* ================================================================
 * L4: Hyperbolic Kepler Equation — M_h = e sinh(H) - H
 *
 * For e > 1. Uses asymptotic initial guess and Newton-Raphson.
 * ================================================================ */

double solve_kepler_hyperbolic(double M_h, double e, double tol) {
    double H, dH;
    int iter;

    if (tol <= 0.0) tol = KEPLER_DEFAULT_TOL;

    /* Handle negative M_h by symmetry: H(-M) = -H(M) */
    if (M_h < 0.0) return -solve_kepler_hyperbolic(-M_h, e, tol);

    if (e <= 1.0) e = 1.000001;

    /* Asymptotic initial guess */
    if (e > 1.8 && M_h > 10.0) {
        H = log(2.0 * M_h / e + 1.8);
    } else if (e < 1.2) {
        H = M_h / (e - 1.0);
    } else {
        H = asinh(M_h / e);
    }

    for (iter = 0; iter < KEPLER_MAX_ITER; iter++) {
        double sh = sinh(H);
        double ch = cosh(H);
        double f_val = e * sh - H - M_h;
        double fp_val = e * ch - 1.0;

        if (fabs(fp_val) < 1e-30) {
            dH = -f_val / (fp_val + 0.5 * e * sh * (f_val / (fp_val + 1e-15)));
        } else {
            dH = -f_val / fp_val;
        }

        if (fabs(dH) > 5.0) dH = (dH > 0 ? 5.0 : -5.0);
        H += dH;
        if (fabs(dH) < tol) break;
    }
    return H;
}

/* ================================================================
 * L4: Parabolic (Barker) Equation — exact analytic solution
 *
 * For e = 1: M_p = tan(nu/2) + (1/3)tan^3(nu/2)
 *
 * Cubic formula: Let A = (3/2)M_p, B = (A + sqrt(A^2+1))^(1/3)
 * Then tan(nu/2) = B - 1/B, nu = 2*atan(B - 1/B)
 *
 * Ref: Barker (1757), Vallado Sec.2.2.5
 * ================================================================ */

double solve_kepler_parabolic(double M_p) {
    double A, B, tan_half;
    A = 1.5 * M_p;
    B = cbrt(A + sqrt(A * A + 1.0));
    tan_half = B - 1.0 / B;
    return 2.0 * atan(tan_half);
}

/* ================================================================
 * L5: Universal Kepler Solver — auto-dispatch by eccentricity
 * ================================================================ */

double solve_kepler_universal(double x, double e, double tol) {
    if (tol <= 0.0) tol = KEPLER_DEFAULT_TOL;
    if (e < 1.0 - 1e-12)
        return solve_kepler(x, e, tol);
    else if (e > 1.0 + 1e-12)
        return solve_kepler_hyperbolic(x, e, tol);
    else
        return solve_kepler_parabolic(x);
}

/* ================================================================
 * Stumpff functions for universal variable formulation
 *
 * c(z) = sum_{k=0}^{inf} (-z)^k / (2k+2)!
 *   c(0)=1/2, c(z>0)=(1-cos(sqrt(z)))/z, c(z<0)=(cosh(sqrt(-z))-1)/(-z)
 *
 * s(z) = sum_{k=0}^{inf} (-z)^k / (2k+3)!
 *   s(0)=1/6, s(z>0)=(sqrt(z)-sin(sqrt(z)))/(z*sqrt(z)),
 *   s(z<0)=(sinh(sqrt(-z))-sqrt(-z))/((-z)*sqrt(-z))
 *
 * Ref: Stumpff (1959), Battin Sec.4.3
 * ================================================================ */

double stumpff_c(double z) {
    if (fabs(z) < 1e-6) {
        return 0.5 - z/24.0 + z*z/720.0 - z*z*z/40320.0;
    }
    if (z > 0.0) {
        double sz = sqrt(z);
        return (1.0 - cos(sz)) / z;
    } else {
        double sz = sqrt(-z);
        return (cosh(sz) - 1.0) / (-z);
    }
}

double stumpff_s(double z) {
    if (fabs(z) < 1e-6) {
        return 1.0/6.0 - z/120.0 + z*z/5040.0 - z*z*z/362880.0;
    }
    if (z > 0.0) {
        double sz = sqrt(z);
        return (sz - sin(sz)) / (z * sz);
    } else {
        double sz = sqrt(-z);
        return (sinh(sz) - sz) / ((-z) * sz);
    }
}

/* ================================================================
 * L5: Universal Variable Formulation (Battin / Goodyear)
 *
 * Solves: sqrt(mu)*t = r0*x + sigma0*x^2*C(alpha*x^2)
 *                       + r0*alpha*x^3*S(alpha*x^2)
 *
 * where x is the universal anomaly, alpha = 1/a,
 * r0 = initial radius, sigma0 = (r dot v)/sqrt(mu).
 *
 * This unified formulation handles elliptic, parabolic, and
 * hyperbolic orbits without singularity at e = 1.
 *
 * Ref: Battin Sec.4.4, Vallado Sec.2.4
 * ================================================================ */

double solve_kepler_universal_variable(double r0, double sigma0,
                                        double dt, double alpha, double tol) {
    double x, dx;
    int iter;

    if (tol <= 0.0) tol = KEPLER_DEFAULT_TOL;

    /* Initial guess */
    x = dt / r0;

    for (iter = 0; iter < KEPLER_MAX_ITER; iter++) {
        double z, C, S, x2, x3, f_val, df_val;
        z = alpha * x * x;
        C = stumpff_c(z);
        S = stumpff_s(z);
        x2 = x * x;
        x3 = x2 * x;

        /* Universal Kepler equation */
        f_val = r0 * x + sigma0 * x2 * C + r0 * alpha * x3 * S - dt;

        /* Derivative: df/dx = r0 + sigma0*x*(1 - z*S) + r0*z*C */
        df_val = r0 + sigma0 * x * (1.0 - z * S) + r0 * z * C;

        if (fabs(df_val) < 1e-30) {
            dx = -f_val / (df_val + 1e-12);
        } else {
            dx = -f_val / df_val;
        }

        if (fabs(dx) > fabs(x) * 0.5) {
            dx = (dx > 0 ? fabs(x) * 0.5 : -fabs(x) * 0.5);
        }

        x += dx;
        if (fabs(dx) < tol) break;
    }
    return x;
}

/* ================================================================
 * L2: Anomaly Conversions — Elliptic
 *
 * tan(nu/2) = sqrt((1+e)/(1-e)) * tan(E/2)
 * tan(E/2) = sqrt((1-e)/(1+e)) * tan(nu/2)
 * M = E - e*sin(E)
 *
 * All use atan2 for proper quadrant resolution.
 * Ref: Murray & Dermott Sec.2.4
 * ================================================================ */

double eccentric_to_true_anomaly(double E, double e) {
    double coeff = sqrt((1.0 + e) / (1.0 - e));
    return 2.0 * atan2(coeff * sin(E / 2.0), cos(E / 2.0));
}

double true_to_eccentric_anomaly(double nu, double e) {
    double coeff = sqrt((1.0 - e) / (1.0 + e));
    return 2.0 * atan2(coeff * sin(nu / 2.0), cos(nu / 2.0));
}

double eccentric_to_mean_anomaly(double E, double e) {
    return E - e * sin(E);
}

double mean_to_eccentric_anomaly(double M, double e, double tol) {
    return solve_kepler(M, e, tol);
}

double true_to_mean_anomaly(double nu, double e) {
    double E = true_to_eccentric_anomaly(nu, e);
    return eccentric_to_mean_anomaly(E, e);
}

double mean_to_true_anomaly(double M, double e, double tol) {
    double E = solve_kepler(M, e, tol);
    return eccentric_to_true_anomaly(E, e);
}

/* ================================================================
 * L2: Anomaly Conversions — Hyperbolic
 *
 * tan(nu/2) = sqrt((e+1)/(e-1)) * tanh(H/2)
 * tanh(H/2) = sqrt((e-1)/(e+1)) * tan(nu/2)
 * ================================================================ */

double hyperbolic_to_true_anomaly(double H, double e) {
    double coeff = sqrt((e + 1.0) / (e - 1.0));
    double tan_nu2 = coeff * tanh(H / 2.0);
    return 2.0 * atan(tan_nu2);
}

double true_to_hyperbolic_anomaly(double nu, double e) {
    double coeff = sqrt((e - 1.0) / (e + 1.0));
    double tan_nu2 = tan(nu / 2.0);
    double tanh_H2 = coeff * tan_nu2;
    if (tanh_H2 >= 1.0) tanh_H2 = 1.0 - 1e-15;
    if (tanh_H2 <= -1.0) tanh_H2 = -1.0 + 1e-15;
    return 2.0 * atanh(tanh_H2);
}

double mean_to_hyperbolic_anomaly(double M_h, double e, double tol) {
    return solve_kepler_hyperbolic(M_h, e, tol);
}

/* ================================================================
 * L4: Orbital Elements -> State Vector (Cartesian)
 *
 * Transforms from perifocal frame (PQW) to reference frame
 * using 3-1-3 Euler rotation: R = R_z(Omega) * R_x(i) * R_z(omega)
 *
 * Perifocal position: r = [r*cos(nu), r*sin(nu), 0]
 *   where r = a(1-e^2)/(1+e*cos(nu))
 * Perifocal velocity: v = sqrt(mu/p) * [-sin(nu), e+cos(nu), 0]
 *   where p = a(1-e^2)
 *
 * Ref: Murray & Dermott Eq.(2.122), Vallado Sec.2.5.1
 * ================================================================ */

void orbital_elements_to_state(const OrbitalElements *el, double mu,
                                Vector3 *r_out, Vector3 *v_out) {
    double a, e, i_val, Omega, omega, nu;
    double r_mag, p, h_mag;
    Vector3 r_pqw, v_pqw;
    Matrix33 R;

    if (!el || !r_out || !v_out) return;

    a = el->a; e = el->e; i_val = el->i;
    Omega = el->Omega; omega = el->omega; nu = el->nu;

    if (fabs(a) < 1e-30 || mu <= 0.0) {
        *r_out = vec3_zero();
        *v_out = vec3_zero();
        return;
    }

    r_mag = a * (1.0 - e * e) / (1.0 + e * cos(nu));
    p = a * (1.0 - e * e);
    h_mag = sqrt(mu * p);

    /* Perifocal position */
    r_pqw.x = r_mag * cos(nu);
    r_pqw.y = r_mag * sin(nu);
    r_pqw.z = 0.0;

    /* Perifocal velocity */
    v_pqw.x = -mu / h_mag * sin(nu);
    v_pqw.y =  mu / h_mag * (e + cos(nu));
    v_pqw.z = 0.0;

    /* 3-1-3 rotation: R = R_z(Omega) * R_x(i) * R_z(omega) */
    {
        Matrix33 Rz_om = mat33_rotz(omega);
        Matrix33 Rx_i  = mat33_rotx(i_val);
        Matrix33 Rz_Om = mat33_rotz(Omega);
        R = mat33_mul(Rz_Om, mat33_mul(Rx_i, Rz_om));
    }

    *r_out = mat33_mul_vec(R, r_pqw);
    *v_out = mat33_mul_vec(R, v_pqw);
}

/* ================================================================
 * L4: State Vector -> Orbital Elements
 *
 * Inverse transformation using angular momentum and eccentricity
 * vectors (Laplace-Runge-Lenz approach).
 *
 * Steps:
 *   1. h = r x v  (specific angular momentum)
 *   2. n = k x h  (node vector, k = [0,0,1])
 *   3. i = acos(h_z / |h|)
 *   4. Omega = atan2(n_y, n_x) with quadrant logic
 *   5. a = 1 / (2/r - v^2/mu)  (vis-viva)
 *   6. e_vec = (v x h)/mu - r_hat
 *   7. omega from n dot e_vec
 *   8. nu from e_vec dot r
 *
 * Quadrant ambiguities resolved via z-component signs.
 *
 * Ref: Vallado Algorithm 9 (Sec.2.5.2), Murray & Dermott Sec.2.8
 * ================================================================ */

void state_to_orbital_elements(Vector3 r, Vector3 v, double mu,
                                OrbitalElements *el) {
    Vector3 h_vec, n_vec, e_vec, k;
    double r_mag, v_mag, h_mag, n_mag, e_mag;
    double xi, a_val, incl, Omega, omega, nu, r_dot_v;

    if (!el) return;

    r_mag = vec3_norm(r);
    v_mag = vec3_norm(v);

    if (r_mag < 1e-30 || v_mag < 1e-30) {
        memset(el, 0, sizeof(*el));
        el->a = INFINITY; el->e = 1.0;
        return;
    }

    /* Angular momentum */
    h_vec = vec3_cross(r, v);
    h_mag = vec3_norm(h_vec);

    /* Node vector: k x h, k = [0,0,1] */
    k = vec3_new(0.0, 0.0, 1.0);
    n_vec = vec3_cross(k, h_vec);
    n_mag = vec3_norm(n_vec);

    /* Semi-major axis from vis-viva */
    xi = 2.0 / r_mag - v_mag*v_mag / mu;
    if (fabs(xi) < 1e-30) {
        a_val = 1e30; /* near-parabolic */
    } else {
        a_val = 1.0 / xi;
    }

    /* Eccentricity vector: e = (v x h)/mu - r_hat */
    {
        Vector3 v_cross_h = vec3_cross(v, h_vec);
        Vector3 r_hat = vec3_scale(1.0 / r_mag, r);
        e_vec = vec3_sub(vec3_scale(1.0 / mu, v_cross_h), r_hat);
    }
    e_mag = vec3_norm(e_vec);

    /* Inclination */
    {
        double cos_i;
        if (h_mag > 1e-30) {
            cos_i = h_vec.z / h_mag;
            cos_i = clamp_val(cos_i, -1.0, 1.0);
            incl = acos(cos_i);
        } else {
            incl = 0.0;
        }
    }

    /* RAAN */
    if (n_mag > 1e-15) {
        double cos_Om = n_vec.x / n_mag;
        cos_Om = clamp_val(cos_Om, -1.0, 1.0);
        Omega = acos(cos_Om);
        if (n_vec.y < 0.0) Omega = TWOPI - Omega;
    } else {
        Omega = 0.0; /* equatorial: RAAN undefined */
    }

    r_dot_v = vec3_dot(r, v);

    /* Argument of periapsis */
    if (e_mag > 1e-12) {
        double cos_om;
        if (n_mag > 1e-15) {
            cos_om = vec3_dot(n_vec, e_vec) / (n_mag * e_mag);
        } else {
            /* Equatorial: define omega from e_vec x-component */
            cos_om = e_vec.x / e_mag;
        }
        cos_om = clamp_val(cos_om, -1.0, 1.0);
        omega = acos(cos_om);
        if (e_vec.z < 0.0) omega = TWOPI - omega;
    } else {
        omega = 0.0; /* circular: argument of periapsis undefined */
    }

    /* True anomaly */
    if (e_mag > 1e-12) {
        double cos_nu = vec3_dot(e_vec, r) / (e_mag * r_mag);
        cos_nu = clamp_val(cos_nu, -1.0, 1.0);
        nu = acos(cos_nu);
        if (r_dot_v < 0.0) nu = TWOPI - nu;
    } else {
        /* Circular orbit: use argument of latitude */
        if (n_mag > 1e-15) {
            double cos_u = vec3_dot(n_vec, r) / (n_mag * r_mag);
            cos_u = clamp_val(cos_u, -1.0, 1.0);
            nu = acos(cos_u);
            if (r.z < 0.0) nu = TWOPI - nu;
        } else {
            nu = atan2(r.y, r.x);
            if (nu < 0.0) nu += TWOPI;
        }
    }

    el->a = a_val;
    el->e = e_mag;
    el->i = incl;
    el->Omega = Omega;
    el->omega = omega;
    el->nu = nu;
}

/* ================================================================
 * Round-trip validation: elements -> state -> elements
 * ================================================================ */

double kepler_roundtrip_error(const OrbitalElements *el, double mu) {
    Vector3 r, v;
    OrbitalElements el2;
    double err_a, err_e, err_i, err_Om, err_om, max_err;

    if (!el) return INFINITY;

    orbital_elements_to_state(el, mu, &r, &v);
    state_to_orbital_elements(r, v, mu, &el2);

    err_a  = fabs(el->a - el2.a);
    err_e  = fabs(el->e - el2.e);
    err_i  = fabs(el->i - el2.i);
    err_Om = fabs(el->Omega - el2.Omega);
    err_om = fabs(el->omega - el2.omega);

    max_err = err_a;
    if (err_e  > max_err) max_err = err_e;
    if (err_i  > max_err) max_err = err_i;
    if (err_Om > max_err) max_err = err_Om;
    if (err_om > max_err) max_err = err_om;
    return max_err;
}

/* ================================================================
 * Orbit state at a given true anomaly
 * ================================================================ */

void orbit_state_at_anomaly(double a, double e, double nu, double mu,
                             double *r_mag, double *v_mag) {
    double denom = 1.0 + e * cos(nu);
    if (fabs(denom) < 1e-30) {
        *r_mag = INFINITY;
        *v_mag = 0.0;
        return;
    }
    *r_mag = a * (1.0 - e * e) / denom;
    *v_mag = sqrt(mu * (2.0 / (*r_mag) - 1.0 / a));
}