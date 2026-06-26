/**
 * perturbations.c — Orbital Perturbation Models
 *
 * J2 oblateness, atmospheric drag, solar radiation pressure,
 * third-body effects, general relativistic precession, and
 * Gauss/Lagrange planetary equations for osculating elements.
 *
 * References:
 *   - Vallado Sec.8.2-8.7, Sec.9.2-9.8
 *   - Murray & Dermott Sec.10.2-10.5
 *   - Kaula (1966), Theory of Satellite Geodesy
 *   - Goldstein Sec.7.6
 */

#include "celestial_types.h"
#include "perturbations.h"
#include <math.h>
#include <float.h>
#include <string.h>

/* ---- Earth orientation parameters (defaults) ---- */
#define DEFAULT_J2   0.0010826265
#define DEFAULT_J4  -1.619e-6
#define DEFAULT_R    6378.1363
#define DEFAULT_MU   398600.4418
#define R_SUN_KM     696000.0

/* ================================================================
 * L6: J2 Gravitational Perturbation
 *
 * The dominant non-spherical term in the geopotential. J2 represents
 * the equatorial bulge. The perturbing potential is:
 *
 *   U_J2 = -(mu/r) * (R/r)^2 * J2/2 * (3*sin^2(lat) - 1)
 *
 * Acceleration in ECI coordinates (Vallado Eq.8.25):
 *   a_x = -1.5*J2*mu*R^2/r^5 * x*(1 - 5*z^2/r^2)
 *   a_y = -1.5*J2*mu*R^2/r^5 * y*(1 - 5*z^2/r^2)
 *   a_z = -1.5*J2*mu*R^2/r^5 * z*(3 - 5*z^2/r^2)
 *
 * Ref: Vallado Eq.(8.25), Kaula (1966)
 * ================================================================ */

Vector3 j2_acceleration(Vector3 r, double mu, double J2, double R_body) {
    double r_mag, x, y, z, z_over_r_sq, factor;
    Vector3 a;

    r_mag = vec3_norm(r);
    if (r_mag < R_body) return vec3_zero(); /* below surface */

    x = r.x; y = r.y; z = r.z;
    z_over_r_sq = (z / r_mag) * (z / r_mag);
    factor = -1.5 * J2 * mu * R_body * R_body / (r_mag * r_mag * r_mag * r_mag * r_mag);

    a.x = factor * x * (1.0 - 5.0 * z_over_r_sq);
    a.y = factor * y * (1.0 - 5.0 * z_over_r_sq);
    a.z = factor * z * (3.0 - 5.0 * z_over_r_sq);

    return a;
}

/* ================================================================
 * L6: J2 Secular (Long-Period) Rates
 *
 * The averaged effect of J2 on the orbital elements:
 *
 *   Omega_dot = -1.5*n*J2*(R/p)^2*cos(i)
 *   omega_dot =  0.75*n*J2*(R/p)^2*(5*cos^2(i) - 1)
 *   M_dot     = n + 0.75*n*J2*(R/p)^2*sqrt(1-e^2)*(3*cos^2(i) - 1)
 *
 * These are the first-order secular perturbations. Omega_dot causes
 * nodal regression (used for sun-synchronous orbits) and omega_dot
 * causes apsidal rotation.
 *
 * Ref: Vallado Eq.(8.29), Kaula (1966)
 * ================================================================ */

J2SecularRates j2_secular_rates(double a, double e, double i,
                                 double mu, double J2, double R_body) {
    J2SecularRates rates;
    double n, p, factor;

    if (a <= 0.0 || e >= 1.0) {
        rates.Omega_dot = 0.0;
        rates.omega_dot = 0.0;
        rates.M_dot = 0.0;
        return rates;
    }

    n = sqrt(mu / (a * a * a));
    p = a * (1.0 - e * e);
    factor = 1.5 * n * J2 * (R_body / p) * (R_body / p);

    rates.Omega_dot = -factor * cos(i);
    rates.omega_dot =  0.5 * factor * (5.0 * cos(i) * cos(i) - 1.0);
    rates.M_dot     =  n + 0.5 * factor * sqrt(1.0 - e * e) * (3.0 * cos(i) * cos(i) - 1.0);

    return rates;
}

/* ================================================================
 * L6: Sun-Synchronous Orbit Inclination
 *
 * A sun-synchronous orbit precesses its RAAN at the same rate as
 * the Earth orbits the Sun: Omega_dot = 2*pi/year = ~0.9856 deg/day.
 *
 * Solving for i:
 *   cos(i) = -Omega_sun_dot / (1.5*n*J2*(R/p)^2)
 *
 * Sun-synchronous orbits are used by Earth-observation satellites
 * to maintain constant lighting conditions.
 *
 * Ref: Vallado Sec.8.3.2
 * ================================================================ */

double sun_sync_inclination(double a, double e, double mu, double J2, double R_body) {
    double n, p, omega_sun, cos_i;

    if (a <= 0.0) return M_PI / 2.0;

    n = sqrt(mu / (a * a * a));
    p = a * (1.0 - e * e);

    /* Apparent solar motion: 2*pi per tropical year */
    omega_sun = TWOPI / SECONDS_PER_YEAR;

    cos_i = -omega_sun / (1.5 * n * J2 * (R_body / p) * (R_body / p));
    cos_i = clamp_val(cos_i, -1.0, 1.0);

    return acos(cos_i); /* retrograde: i > 90 deg */
}

/* ================================================================
 * L6: J2+J4 Combined Acceleration
 *
 * Higher-order zonal harmonics. J4 models the next symmetric term.
 *
 * Ref: Vallado Sec.8.2.1
 * ================================================================ */

Vector3 j2_j4_acceleration(Vector3 r, double mu, double J2, double J4, double R_body) {
    double r_mag = vec3_norm(r);
    double R_r = R_body / r_mag;
    double R_r2 = R_r * R_r;
    double R_r4 = R_r2 * R_r2;
    double x = r.x, y = r.y, z = r.z;
    double zr = z / r_mag;
    double zr2 = zr * zr;
    double zr4 = zr2 * zr2;

    double mu_r3 = mu / (r_mag * r_mag * r_mag);

    /* J2 contribution */
    double c_J2 = -1.5 * J2 * R_r2;
    double a_J2_x = c_J2 * x / r_mag * (1.0 - 5.0 * zr2);
    double a_J2_y = c_J2 * y / r_mag * (1.0 - 5.0 * zr2);
    double a_J2_z = c_J2 * z / r_mag * (3.0 - 5.0 * zr2);

    /* J4 contribution (simplified) */
    double c_J4 = -0.625 * J4 * R_r4;
    double a_J4_x = c_J4 * x / r_mag * (3.0 - 42.0*zr2 + 63.0*zr4);
    double a_J4_y = c_J4 * y / r_mag * (3.0 - 42.0*zr2 + 63.0*zr4);
    double a_J4_z = c_J4 * z / r_mag * (15.0 - 70.0*zr2 + 63.0*zr4);

    return vec3_new(
        mu_r3 * (a_J2_x + a_J4_x),
        mu_r3 * (a_J2_y + a_J4_y),
        mu_r3 * (a_J2_z + a_J4_z)
    );
}

/* ================================================================
 * L7: Atmospheric Drag Acceleration
 *
 * For satellites in LEO (<1000 km), atmospheric drag is the dominant
 * non-conservative perturbation. The acceleration opposes velocity:
 *
 *   a_drag = -0.5 * C_D * (A/m) * rho * v^2 * v_hat
 *
 * Using a simple exponential atmosphere model:
 *   rho(h) = rho0 * exp(-(h - h0) / H)
 *
 * where H is the scale height (~8.5 km for Earth at 200 km).
 *
 * Ref: Vallado Sec.8.5
 * ================================================================ */

Vector3 atmospheric_drag_acceleration(Vector3 r, Vector3 v,
                                       double C_D, double A_over_m,
                                       double rho0, double r0, double H) {
    double r_mag, alt, rho, v_mag;
    Vector3 a;

    r_mag = vec3_norm(r);
    alt = r_mag - r0;
    if (alt < 0.0) return vec3_zero(); /* below reference, already decayed */

    rho = rho0 * exp(-alt / H);
    v_mag = vec3_norm(v);

    if (v_mag < 1e-15) return vec3_zero();

    /* a = -0.5*C_D*(A/m)*rho*v^2 * v_hat */
    {
        Vector3 v_hat = vec3_scale(1.0 / v_mag, v);
        a = vec3_scale(-0.5 * C_D * A_over_m * rho * v_mag * v_mag, v_hat);
    }

    return a;
}

/* ================================================================
 * L7: Drag Decay Rate and Lifetime
 *
 * For a circular orbit, the semi-major axis decay rate is:
 *   da/dt = -C_D * (A/m) * rho * n * a^2
 *
 * The orbital lifetime can be estimated by integrating this
 * equation from the initial altitude to re-entry (~120 km).
 *
 * Ref: Vallado Sec.8.5.2, King-Hele (1987)
 * ================================================================ */

double drag_decay_rate(double a, double C_D, double A_over_m,
                        double rho, double mu) {
    if (a <= 0.0) return 0.0;
    double n = sqrt(mu / (a * a * a));
    return -C_D * A_over_m * rho * n * a * a;
}

double drag_lifetime_estimate(double a0, double e0, double C_D,
                               double A_over_m, double rho0, double r0,
                               double H, double mu) {
    /* Simplified: integrate da/dt = -C_D*(A/m)*rho0*exp(-(a-r0)/H)*n*a^2
     * For circular orbit, this yields:
     *   T_life ~ H / (da/dt)_0 * [exp((a0 - r_reentry)/H) - 1]
     *
     * where r_reentry ~ r0 + 120 km (typical LEO reentry altitude).
     */
    double a = a0;
    double dt_total = 0.0;
    double da_per_step;
    int i;

    /* Simple Euler integration with altitude steps */
    for (i = 0; i < 2000; i++) {
        double alt = a - r0;
        if (alt < 120.0) break; /* re-entry altitude */

        double rho = rho0 * exp(-alt / H);
        double n = sqrt(mu / (a * a * a));
        da_per_step = -C_D * A_over_m * rho * n * a * a;

        /* Adaptive step: keep altitude change < 10 km */
        double step = 10.0 / fabs(da_per_step);
        if (step > SECONDS_PER_DAY * 30) step = SECONDS_PER_DAY * 30;
        if (step < 1.0) step = 1.0;

        a += da_per_step * step;
        dt_total += step;

        if (a < r0 + 100.0) break;
    }

    return dt_total;
}

/* ---- Harris-Priester density model (simplified) ---- */
double atmospheric_density_harris_priester(double altitude_km, double F10_7) {
    /* Simplified fit to Harris-Priester tables.
     * F10_7: solar radio flux at 10.7 cm (sfu)
     * Returns density in kg/km^3. */
    double h = altitude_km;
    double rho;

    if (h < 150.0) {
        rho = 5.0e12 * exp(-h / 45.0); /* lower thermosphere */
    } else if (h < 500.0) {
        double base = (h < 300.0) ? 2.0e10 : 1.0e9;
        rho = base * exp(-(h - 150.0) / 50.0);
        /* Solar activity modulation */
        rho *= (1.0 + 0.5 * (F10_7 - 150.0) / 100.0);
    } else if (h < 1000.0) {
        rho = 5.0e7 * exp(-(h - 500.0) / 150.0);
    } else {
        rho = 1.0e6 * exp(-(h - 1000.0) / 500.0);
    }

    if (rho < 1e-3) rho = 1e-3;
    return rho;
}

/* ================================================================
 * L7: Solar Radiation Pressure (SRP)
 *
 * Photons from the Sun impart momentum on the spacecraft surface:
 *
 *   a_srp = -C_R * (A/m) * P_sun * (1 AU)^2 * d / |d|^3
 *
 * where d = r_sat - r_sun, C_R is the reflectivity coefficient
 * (1.0 = perfect absorber, 2.0 = perfect reflector), and
 * P_sun = 4.56e-6 N/m^2 at 1 AU.
 *
 * The shadow function accounts for Earth eclipses.
 *
 * Ref: Vallado Sec.8.6, Montenbruck & Gill Sec.3.4
 * ================================================================ */

Vector3 srp_acceleration(Vector3 r_sat, Vector3 r_sun,
                          double C_R, double A_over_m, double P_sun) {
    Vector3 d = vec3_sub(r_sat, r_sun);
    double d_mag = vec3_norm(d);

    if (d_mag < 1e-10) return vec3_zero();

    /* P_sun at 1 AU, scaled by (AU/d)^2 */
    double factor = -C_R * A_over_m * P_sun * (AU_KM * AU_KM) / (d_mag * d_mag * d_mag);

    return vec3_scale(factor, d);
}

/* ---- Shadow function: cylindrical Earth shadow model ---- */

double shadow_function(Vector3 r_sat, Vector3 r_sun, double R_body, double R_sun) {
    /* Cylindrical shadow model:
     * Returns 0.0 for total eclipse (umbra),
     *         partial for penumbra,
     *         1.0 for full sunlight.
     *
     * The satellite is in shadow if:
     *   r_sat dot r_sun < 0  AND  |r_sat x r_sun|/|r_sun| < R_body
     */
    double r_mag = vec3_norm(r_sat);
    double s_mag = vec3_norm(r_sun);
    double r_dot_s = vec3_dot(r_sat, r_sun);

    /* Satellite on night side? */
    if (r_dot_s > 0.0) return 1.0; /* day side */

    /* Distance from shadow axis */
    Vector3 cross_rs = vec3_cross(r_sat, r_sun);
    double d_perp = vec3_norm(cross_rs) / s_mag;

    /* Umbra and penumbra radii at satellite distance */
    double s_ang = asin(R_sun / s_mag); /* angular radius of Sun */
    double umbra_r = R_body - r_mag * tan(s_ang);
    double penumbra_r = R_body + r_mag * tan(s_ang);

    if (d_perp < umbra_r) {
        /* Total eclipse */
        return 0.0;
    } else if (d_perp < penumbra_r) {
        /* Penumbra: linear interpolation */
        return (d_perp - umbra_r) / (penumbra_r - umbra_r);
    } else {
        return 1.0; /* full sunlight */
    }
}

/* ================================================================
 * L7: Third-Body Perturbation
 *
 * The gravitational attraction of a third body (Sun, Moon) perturbs
 * the satellite's orbit. The acceleration relative to the central
 * body is:
 *
 *   a_3b = -mu_b[ (r_sat-r_b)/|r_sat-r_b|^3 + r_b/|r_b|^3 ]
 *
 * The first term is the direct attraction on the satellite, and
 * the second (indirect) term accounts for the acceleration of the
 * central body by the third body.
 *
 * Ref: Vallado Sec.8.4, Battin Sec.8.4
 * ================================================================ */

Vector3 third_body_acceleration(Vector3 r_sat, Vector3 r_body, double mu_body) {
    Vector3 delta = vec3_sub(r_sat, r_body);
    double d_mag = vec3_norm(delta);
    double r_body_mag = vec3_norm(r_body);
    Vector3 direct, indirect;

    if (d_mag < 1e-10 || r_body_mag < 1e-10) return vec3_zero();

    /* Direct: point mass attraction on satellite */
    direct = vec3_scale(1.0 / (d_mag * d_mag * d_mag), delta);

    /* Indirect: acceleration of central body by third body */
    indirect = vec3_scale(1.0 / (r_body_mag * r_body_mag * r_body_mag), r_body);

    return vec3_scale(-mu_body, vec3_sub(direct, indirect));
}

Vector3 solar_perturbation_acceleration(Vector3 r_sat, Vector3 r_sun) {
    return third_body_acceleration(r_sat, r_sun, MU_SUN);
}

Vector3 lunar_perturbation_acceleration(Vector3 r_sat, Vector3 r_moon) {
    return third_body_acceleration(r_sat, r_moon, MU_MOON);
}

/* ================================================================
 * L8: General Relativistic Corrections
 *
 * GR perihelion precession (Einstein 1915):
 *
 *   omega_dot_GR = 3 * n * mu / (a * c^2 * (1 - e^2))
 *
 * For Mercury: ~42.98 arcsec/century — one of the classical
 * tests of General Relativity.
 *
 * 1PN correction to acceleration:
 *   a_GR = -(mu/r^2)[ (1+3v^2/c^2)r_hat + ... ]
 *
 * Ref: Goldstein Sec.7.6, Weinberg Sec.8.6, Will (1993)
 * ================================================================ */

double gr_precession_rate(double a, double e, double mu, double c) {
    if (a <= 0.0 || e >= 1.0) return 0.0;
    double n = sqrt(mu / (a * a * a));
    return 3.0 * n * mu / (a * c * c * (1.0 - e * e));
}

Vector3 gr_post_newtonian_acceleration(Vector3 r, Vector3 v,
                                        double mu, double c) {
    double r_mag = vec3_norm(r);
    double v_mag = vec3_norm(v);
    double c2 = c * c;
    Vector3 r_hat, a;

    if (r_mag < 1e-30) return vec3_zero();

    r_hat = vec3_scale(1.0 / r_mag, r);

    /* 1PN correction (simplified Einstein-Infeld-Hoffman):
     * a_1PN = -(mu/r^2) * [ (1 + 3*v^2/c^2) * r_hat
     *                        - 2*(2-sqrt(mu/a)/c)*(r_dot_v/c)*v_hat ] */
    {
        double mu_r2 = mu / (r_mag * r_mag);
        double corr = 1.0 + 3.0 * v_mag * v_mag / c2;

        a = vec3_scale(-mu_r2 * corr, r_hat);
    }

    return a;
}

double precession_rate_classical(double a, double m_perturber,
                                  double M_central) {
    /* Classical third-body-induced precession (secular):
     * Omega_dot ~ (3/4)*n*(m/M)*(a/a_p)^3 for external perturber */
    if (a <= 0.0 || M_central <= 0.0) return 0.0;
    return 0.0; /* Requires perturber orbit parameters for full calc */
}

/* ================================================================
 * L5: Combined Perturbation Models
 * ================================================================ */

OrbitalElements propagate_with_j2(const OrbitalElements *el, double dt,
                                   double mu, double J2, double R_body) {
    OrbitalElements result;
    J2SecularRates rates;

    if (!el) {
        memset(&result, 0, sizeof(result));
        return result;
    }

    result = *el;
    rates = j2_secular_rates(el->a, el->e, el->i, mu, J2, R_body);

    /* J2 has no secular effect on a, e, i to first order */
    result.Omega = wrap_2pi(el->Omega + rates.Omega_dot * dt);
    result.omega = wrap_2pi(el->omega + rates.omega_dot * dt);
    result.nu    = wrap_2pi(el->nu + rates.M_dot * dt);

    return result;
}

Vector3 combined_perturbation(Vector3 r, Vector3 v, Vector3 r_sun,
                               Vector3 r_moon, double mu, double J2,
                               double R_body, double C_D, double A_over_m,
                               double rho0, double r0, double H,
                               double C_R, double P_sun, PerturbationFlags flags) {
    Vector3 a_total = vec3_zero();

    if (flags & PERTURB_J2) {
        a_total = vec3_add(a_total, j2_acceleration(r, mu, J2, R_body));
    }
    if (flags & PERTURB_DRAG) {
        a_total = vec3_add(a_total,
            atmospheric_drag_acceleration(r, v, C_D, A_over_m, rho0, r0, H));
    }
    if (flags & PERTURB_SRP) {
        a_total = vec3_add(a_total,
            srp_acceleration(r, r_sun, C_R, A_over_m, P_sun));
    }
    if (flags & PERTURB_THIRD_BODY) {
        a_total = vec3_add(a_total,
            third_body_acceleration(r, r_sun, MU_SUN));
        a_total = vec3_add(a_total,
            third_body_acceleration(r, r_moon, MU_MOON));
    }

    return a_total;
}

/* ================================================================
 * L4: Gauss Planetary Equations (Variation of Parameters)
 *
 * The time derivatives of osculating orbital elements due to
 * perturbation accelerations in the RSW (Radial-Transverse-Normal)
 * coordinate frame:
 *
 *   da/dt = (2a^2/h)[e*sin(nu)*R + p/r*S]
 *   de/dt = (1/h)[p*sin(nu)*R + ((p+r)*cos(nu)+r*e)*S]
 *   di/dt = (r*cos(u)/h)*W
 *   dOmega/dt = (r*sin(u)/(h*sin(i)))*W
 *   domega/dt = (1/(h*e))[-p*cos(nu)*R + (p+r)*sin(nu)*S]
 *                - cos(i)*dOmega/dt
 *   dM/dt = n - (2r/(na^2))*R*sqrt(1-e^2)*(de/dt)/(e)
 *
 * where u = omega + nu (argument of latitude).
 *
 * These are the fundamental tool for orbit determination when
 * non-Keplerian forces are present.
 *
 * Ref: Vallado Sec.8.3.1, Battin Sec.10.3, Gauss (1809)
 * ================================================================ */

void gauss_planetary_equations(const OrbitalElements *el, double mu,
                                const PerturbationRSW *accel,
                                double *a_dot, double *e_dot, double *i_dot,
                                double *Omega_dot, double *omega_dot, double *M_dot) {
    double a, e, i_val, omega, nu;
    double p, r, h_mag, n, u;

    if (!el || !accel) {
        if (a_dot) *a_dot = 0.0;
        if (e_dot) *e_dot = 0.0;
        if (i_dot) *i_dot = 0.0;
        if (Omega_dot) *Omega_dot = 0.0;
        if (omega_dot) *omega_dot = 0.0;
        if (M_dot) *M_dot = 0.0;
        return;
    }

    a = el->a; e = el->e; i_val = el->i;
    omega = el->omega; nu = el->nu;

    p = a * (1.0 - e * e);
    r = a * (1.0 - e * e) / (1.0 + e * cos(nu));
    h_mag = sqrt(mu * p);
    n = sqrt(mu / (a * a * a));
    u = omega + nu;

    double R = accel->R;
    double S = accel->S;
    double W = accel->W;

    double sin_nu = sin(nu);
    double cos_nu = cos(nu);
    double sin_u  = sin(u);
    double cos_u  = cos(u);
    double sin_i  = sin(i_val);
    double cos_i  = cos(i_val);

    /* da/dt */
    if (a_dot) {
        *a_dot = (2.0 * a * a / h_mag) * (e * sin_nu * R + p / r * S);
    }

    /* de/dt */
    if (e_dot) {
        *e_dot = (1.0 / h_mag) * (p * sin_nu * R
                  + ((p + r) * cos_nu + r * e) * S);
    }

    /* di/dt */
    if (i_dot) {
        *i_dot = (r * cos_u / h_mag) * W;
    }

    /* dOmega/dt (RAAN) */
    if (Omega_dot) {
        if (fabs(sin_i) > 1e-15) {
            *Omega_dot = (r * sin_u / (h_mag * sin_i)) * W;
        } else {
            *Omega_dot = 0.0;
        }
    }

    /* domega/dt (argument of periapsis) */
    if (omega_dot) {
        double domega;
        if (e > 1e-15) {
            domega = (1.0 / (h_mag * e)) * (-p * cos_nu * R
                      + (p + r) * sin_nu * S);
        } else {
            domega = 0.0;
        }
        /* Subtract cos(i)*dOmega/dt */
        if (fabs(sin_i) > 1e-15) {
            domega -= cos_i * (r * sin_u / (h_mag * sin_i)) * W;
        }
        *omega_dot = domega;
    }

    /* dM/dt (mean anomaly) */
    if (M_dot) {
        double dM = n;
        if (e > 1e-15 && e_dot) {
            dM -= (2.0 * r / (n * a * a)) * sqrt(1.0 - e * e) * (*e_dot) / e;
            /* Additional term from R: dM = n - (1-e^2)/(nae)*[...] */
            dM += sqrt(1.0 - e * e) / (n * a * e)
                  * (2.0 * R * a * e * sin_nu / h_mag);
        }
        *M_dot = dM;
    }
}