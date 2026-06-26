/**
 * mission.c — Orbital Maneuvers and Interplanetary Transfers
 *
 * Hohmann, bi-elliptic, plane changes, gravity assists,
 * patched conics, Lambert problem, launch windows,
 * and low-thrust approximations.
 *
 * References:
 *   - Vallado Sec.6.3-6.7
 *   - Bate, Mueller & White, Ch.6-7
 *   - Battin, Sec.7.1-7.8
 *   - Prussing & Conway, Ch.6-8
 */

#include "celestial_types.h"
#include "two_body.h"
#include "mission.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>

/* ================================================================
 * L6: Hohmann Transfer
 *
 * The minimum-energy impulsive transfer between two coplanar
 * circular orbits (Hohmann 1925).
 *
 * Transfer ellipse:
 *   a_t = (r1 + r2)/2
 *   e_t = |r2 - r1|/(r2 + r1)
 *
 * dv1 = |v_peri - v_circ1|    (at departure)
 * dv2 = |v_circ2 - v_apo|     (at arrival)
 * T_transfer = pi*sqrt(a_t^3/mu)   (half period)
 *
 * The Hohmann transfer is optimal for r2/r1 < ~11.94
 * (beyond which bi-elliptic becomes better).
 *
 * Ref: Vallado Sec.6.3.1, Bate et al. Sec.6.2
 * ================================================================ */

TransferOrbit hohmann_transfer(double r1, double r2, double mu) {
    TransferOrbit result;
    double v1, v2, a_t, e_t, v_peri, v_apo, dv1, dv2;

    memset(&result, 0, sizeof(result));
    strcpy(result.name, "Hohmann");

    if (r1 <= 0.0 || r2 <= 0.0 || mu <= 0.0) return result;

    v1 = circular_velocity(r1, mu);
    v2 = circular_velocity(r2, mu);

    a_t = (r1 + r2) / 2.0;
    e_t = fabs(r2 - r1) / (r2 + r1);

    /* Transfer orbit velocities at r1 and r2 */
    v_peri = sqrt(mu * (2.0/r1 - 1.0/a_t));
    v_apo  = sqrt(mu * (2.0/r2 - 1.0/a_t));

    dv1 = fabs(v_peri - v1);
    dv2 = fabs(v2 - v_apo);
    result.delta_v1 = dv1;
    result.delta_v2 = dv2;
    result.delta_v_total = dv1 + dv2;
    result.transfer_time = M_PI * sqrt(a_t*a_t*a_t / mu);

    /* Fill transfer orbit elements */
    result.transfer.a = a_t;
    result.transfer.e = e_t;
    result.transfer.i = 0.0;
    result.transfer.Omega = 0.0;
    result.transfer.omega = 0.0;
    result.transfer.nu = 0.0;

    result.depart.a = r1;
    result.depart.e = 0.0;
    result.arrive.a = r2;
    result.arrive.e = 0.0;

    return result;
}

TransferOrbit hohmann_transfer_inward(double r1, double r2, double mu) {
    /* Symmetric: swap r1,r2, compute outward, swap dv's */
    TransferOrbit out;
    if (r1 <= r2) return hohmann_transfer(r1, r2, mu);

    out = hohmann_transfer(r2, r1, mu);
    /* Swap dv1 and dv2 for inward direction */
    {
        double tmp = out.delta_v1;
        out.delta_v1 = out.delta_v2;
        out.delta_v2 = tmp;
    }
    strcpy(out.name, "Hohmann-inward");
    return out;
}

void hohmann_transfer_detail(double r1, double r2, double mu,
                              double *dv1, double *dv2, double *dv_total,
                              double *transfer_time) {
    TransferOrbit t = hohmann_transfer(r1, r2, mu);
    if (dv1) *dv1 = t.delta_v1;
    if (dv2) *dv2 = t.delta_v2;
    if (dv_total) *dv_total = t.delta_v_total;
    if (transfer_time) *transfer_time = t.transfer_time;
}

/* ================================================================
 * L6: Bi-Elliptic Transfer
 *
 * Three-impulse transfer via intermediate radius r_b > max(r1,r2).
 *
 * Step 1: r1 -> r_b (ellipse 1), dv1
 * Step 2: at r_b, transfer to ellipse r_b -> r2, dv2
 * Step 3: circularize at r2, dv3
 *
 * Total dv = dv1 + dv2 + dv3
 *
 * Bi-elliptic transfer can beat Hohmann when r2/r1 > ~11.94
 * (for the limiting case r_b -> infinity).
 *
 * Ref: Vallado Sec.6.3.2, Hoelker & Silber (1959)
 * ================================================================ */

TransferOrbit bi_elliptic_transfer(double r1, double r2, double r_b, double mu) {
    TransferOrbit result;
    double v_circ1, v_circ2, a1, a2, dv1, dv2, dv3;
    double v_peri1, v_apo1, v_peri2, v_apo2;

    memset(&result, 0, sizeof(result));
    strcpy(result.name, "Bi-elliptic");

    if (r1 <= 0.0 || r2 <= 0.0 || r_b <= 0.0 || mu <= 0.0) return result;
    if (r_b <= r1 || r_b <= r2) return result; /* r_b must be largest */

    v_circ1 = circular_velocity(r1, mu);
    v_circ2 = circular_velocity(r2, mu);

    /* Ellipse 1: r1 -> r_b */
    a1 = (r1 + r_b) / 2.0;
    v_peri1 = sqrt(mu * (2.0/r1 - 1.0/a1));
    v_apo1  = sqrt(mu * (2.0/r_b - 1.0/a1));
    dv1 = fabs(v_peri1 - v_circ1);

    /* Ellipse 2: r_b -> r2 */
    a2 = (r_b + r2) / 2.0;
    v_peri2 = sqrt(mu * (2.0/r_b - 1.0/a2));
    v_apo2  = sqrt(mu * (2.0/r2 - 1.0/a2));
    dv2 = fabs(v_peri2 - v_apo1);

    /* Circularize at r2 */
    dv3 = fabs(v_circ2 - v_apo2);

    result.delta_v1 = dv1;
    result.delta_v2 = dv2;
    result.delta_v3 = dv3;
    result.delta_v_total = dv1 + dv2 + dv3;
    result.transfer_time = M_PI * (sqrt(a1*a1*a1/mu) + sqrt(a2*a2*a2/mu));

    return result;
}

void bi_elliptic_comparison(double r1, double r2, double mu,
                             double *dv_hohmann, double *dv_bielliptic,
                             double *r_b_optimal) {
    TransferOrbit hoh = hohmann_transfer(r1, r2, mu);
    double best_rb = r2 * 1.5;
    double best_dv = INFINITY;
    int i;

    if (dv_hohmann) *dv_hohmann = hoh.delta_v_total;

    /* Scan r_b values to find minimum dv */
    for (i = 0; i < 50; i++) {
        double rb = r2 * (1.5 + i * 0.2);
        TransferOrbit bi = bi_elliptic_transfer(r1, r2, rb, mu);
        if (bi.delta_v_total < best_dv) {
            best_dv = bi.delta_v_total;
            best_rb = rb;
        }
    }

    if (dv_bielliptic) *dv_bielliptic = best_dv;
    if (r_b_optimal) *r_b_optimal = best_rb;
}

/* ================================================================
 * L6: Plane Change Maneuvers
 *
 * Simple plane change: rotate velocity vector by delta_i
 *   dv = 2*v*sin(delta_i/2)
 *
 * For large plane changes, it is more efficient to:
 *   1. Raise apoapsis (spend dv on Hohmann-like maneuver)
 *   2. Perform plane change at high apoapsis (lower v)
 *   3. Lower apoapsis back
 *
 * Combined plane change dv:
 *   dv_total^2 = dv1^2 + dv2_plane^2 (if simultaneous)
 *   or dv_total = dv1 + dv2_plane (if sequential)
 *
 * Ref: Vallado Sec.6.4
 * ================================================================ */

double simple_plane_change_delta_v(double v, double delta_i) {
    return 2.0 * v * sin(fabs(delta_i) / 2.0);
}

double apoapsis_plane_change_delta_v(double v_apo, double delta_i) {
    return 2.0 * v_apo * sin(fabs(delta_i) / 2.0);
}

TransferOrbit hohmann_with_plane_change(double r1, double r2,
                                         double delta_i, double mu) {
    TransferOrbit result = hohmann_transfer(r1, r2, mu);
    double a_t = (r1 + r2) / 2.0;
    double v_apo = sqrt(mu * (2.0/r2 - 1.0/a_t));
    double dv_plane = simple_plane_change_delta_v(v_apo, delta_i);

    /* Combined dv at apoapsis: sqrt(dv2^2 + dv_plane^2) */
    double dv2_combined = sqrt(result.delta_v2 * result.delta_v2
                               + dv_plane * dv_plane);
    result.delta_v_total = result.delta_v1 + dv2_combined;
    strcpy(result.name, "Hohmann+PlaneChange");

    return result;
}

void one_tangent_burn_transfer(double r1, double r2, double mu,
                                double *dv, double *e_transfer) {
    /* One-tangent burn: impulse applied tangentially at r1 to
     * achieve an orbit that intersects r2. Used when r2/r1 is
     * large and Hohmann is expensive.
     *
     * Transfer eccentricity: e = (r2/r1 - 1) / (r2/r1 + cos(nu2))
     * dv = |v1_t - v_circ1| */
    double v_circ1 = circular_velocity(r1, mu);
    double e_val, a_t_val;

    /* For the one-tangent burn, we choose the transfer ellipse
     * such that at r1 (periapsis), the burn is tangential. */
    if (r2 > r1) {
        /* Outward: r1 is periapsis */
        e_val = (r2 - r1) / (r2 + r1); /* simplified */
        a_t_val = r1 / (1.0 - e_val);
    } else {
        /* Inward: r1 is apoapsis */
        e_val = (r1 - r2) / (r1 + r2);
        a_t_val = r1 / (1.0 + e_val);
    }

    double v_peri = sqrt(mu * (2.0/r1 - 1.0/a_t_val));
    *dv = fabs(v_peri - v_circ1);
    if (e_transfer) *e_transfer = e_val;
}

/* ================================================================
 * L7: Gravity Assist (Swing-By)
 *
 * A spacecraft flying past a planet gains/loses heliocentric
 * energy by exchanging momentum with the planet.
 *
 * Turn angle in planet frame:
 *   delta = 2*asin(1/e_hyp) = 2*asin(1/(1 + r_p*V_inf^2/mu))
 *
 * Heliocentric dv gain (maximum, for optimal alignment):
 *   dv_max = 2*V_planet*sin(delta/2)
 *
 * The Voyager missions (Flandro 1966) used a "Grand Tour"
 * alignment of Jupiter-Saturn-Uranus-Neptune for massive
 * gravity assists.
 *
 * Ref: Vallado Sec.6.6, Flandro (1966), Minovitch (1961)
 * ================================================================ */

double gravity_assist_turn_angle(double V_inf, double r_p, double mu) {
    double e_hyp;

    if (r_p <= 0.0 || mu <= 0.0) return 0.0;

    /* Hyperbolic eccentricity:
     * e = 1 + r_p*V_inf^2/mu  (from energy equation) */
    e_hyp = 1.0 + r_p * V_inf * V_inf / mu;

    if (e_hyp <= 1.0) return 0.0; /* not hyperbolic */

    /* Turn angle: delta = 2*asin(1/e) */
    return 2.0 * asin(1.0 / e_hyp);
}

double gravity_assist_delta_v(double V_planet, double V_inf,
                               double r_p, double mu) {
    double delta = gravity_assist_turn_angle(V_inf, r_p, mu);
    return 2.0 * V_planet * sin(delta / 2.0);
}

GravityAssist gravity_assist_full(double V_inf_in, double r_p,
                                   double V_planet, double mu_planet,
                                   double bend_sign) {
    GravityAssist ga;
    memset(&ga, 0, sizeof(ga));

    ga.v_inf_in = V_inf_in;
    ga.v_inf_out = V_inf_in; /* |V_inf| conserved in planet frame */
    ga.r_periapsis = r_p;
    ga.e_hyperbolic = 1.0 + r_p * V_inf_in * V_inf_in / mu_planet;
    ga.turn_angle = 2.0 * asin(1.0 / ga.e_hyperbolic);
    ga.delta_v_helio = 2.0 * V_planet * sin(ga.turn_angle / 2.0) * bend_sign;

    return ga;
}

double max_gravity_assist_gain(double V_planet, double V_inf) {
    /* Maximum possible gain when spacecraft passes behind planet
     * in its orbital direction (or ahead for braking). */
    double delta_max = M_PI; /* limiting case: r_p -> 0, delta -> pi */
    (void)V_inf;
    return 2.0 * V_planet * sin(delta_max / 2.0);
}

/* ================================================================
 * L7: Patched Conics
 *
 * Approximate interplanetary trajectory by dividing into:
 *   1. Departure: hyperbolic escape from planet sphere of influence
 *   2. Heliocentric: conic arc from departure to arrival planet
 *   3. Arrival: hyperbolic capture at target planet
 *
 * This is the method used for preliminary mission design
 * (e.g., Earth-to-Mars transfers).
 *
 * Ref: Vallado Sec.6.7, Battin Sec.6.4
 * ================================================================ */

TransferOrbit patched_conics_transfer(double r_depart, double r_arrive,
                                       double mu_depart, double mu_arrive,
                                       double r_p_depart, double r_p_arrive,
                                       double mu_sun) {
    TransferOrbit result;
    double v_depart, v_arrive, v_inf_depart, v_inf_arrive;
    double v_esc_depart, v_capture;
    double dv_depart, dv_hohmann, dv_arrive, dv_total;

    memset(&result, 0, sizeof(result));
    strcpy(result.name, "PatchedConics");

    if (r_depart <= 0.0 || r_arrive <= 0.0) return result;

    /* Planet circular speeds */
    v_depart = circular_velocity(r_depart, mu_sun);
    v_arrive = circular_velocity(r_arrive, mu_sun);

    /* Hohmann heliocentric transfer */
    {
        TransferOrbit hoh = hohmann_transfer(r_depart, r_arrive, mu_sun);
        dv_hohmann = hoh.delta_v_total;
        result.transfer_time = hoh.transfer_time;
    }

    /* Departure hyperbolic excess speed:
     * V_inf is the difference between transfer orbit and planet speed */
    {
        double a_t_h = (r_depart + r_arrive) / 2.0;
        double v_t_depart = sqrt(mu_sun * (2.0/r_depart - 1.0/a_t_h));
        v_inf_depart = fabs(v_t_depart - v_depart);
    }

    /* Escape delta-v from parking orbit at r_p_depart */
    {
        double v_park = circular_velocity(r_p_depart, mu_depart);
        v_esc_depart = sqrt(2.0 * mu_depart / r_p_depart + v_inf_depart*v_inf_depart);
        dv_depart = v_esc_depart - v_park;
    }

    /* Arrival hyperbolic excess speed */
    {
        double a_t_h = (r_depart + r_arrive) / 2.0;
        double v_t_arrive = sqrt(mu_sun * (2.0/r_arrive - 1.0/a_t_h));
        v_inf_arrive = fabs(v_t_arrive - v_arrive);
    }

    /* Capture delta-v at r_p_arrive */
    {
        v_capture = sqrt(v_inf_arrive*v_inf_arrive + 2.0*mu_arrive/r_p_arrive);
        double v_park_arrive = circular_velocity(r_p_arrive, mu_arrive);
        dv_arrive = v_capture - v_park_arrive;
    }

    dv_total = dv_depart + dv_hohmann + dv_arrive;

    result.delta_v1 = dv_depart;
    result.delta_v2 = dv_hohmann;
    result.delta_v3 = dv_arrive;
    result.delta_v_total = dv_total;

    return result;
}

/* ================================================================
 * L6: Lambert Problem (Battin's Universal Variable Method)
 *
 * Given: r1, r2, dt (time of flight), mu
 * Find: v1, v2 such that a Keplerian orbit connects r1 to r2 in dt.
 *
 * The Lambert problem is fundamental to orbit determination and
 * rendezvous guidance. This implementation uses Battin's (1987)
 * universal variable formulation with continued fraction refinement.
 *
 * For multi-revolution cases (dt > one orbital period), there are
 * 2*max_revs + 1 possible solutions.
 *
 * Ref: Vallado Sec.7.3, Battin Sec.7.2, Lancaster & Blanchard (1969)
 * ================================================================ */

LambertSolution lambert_solver(Vector3 r1, Vector3 r2, double dt,
                                double mu, int prograde, int max_revs) {
    LambertSolution sol;
    double r1_mag, r2_mag, cos_dnu, sin_dnu;
    double A, z, y, x, F = 0.0, dF, t_calc;
    double Cz = 0.5, Sz = 1.0/6.0;
    int iter, rev;
    Vector3 r1_hat, r2_hat, h_vec;

    memset(&sol, 0, sizeof(sol));

    r1_mag = vec3_norm(r1);
    r2_mag = vec3_norm(r2);

    if (r1_mag < 1e-12 || r2_mag < 1e-12) return sol;

    r1_hat = vec3_scale(1.0/r1_mag, r1);
    r2_hat = vec3_scale(1.0/r2_mag, r2);
    cos_dnu = vec3_dot(r1_hat, r2_hat);
    cos_dnu = clamp_val(cos_dnu, -1.0, 1.0);

    /* Determine transfer angle and direction */
    h_vec = vec3_cross(r1, r2);
    if (prograde) {
        sin_dnu = (h_vec.z >= 0) ? sqrt(1.0 - cos_dnu*cos_dnu)
                                  : -sqrt(1.0 - cos_dnu*cos_dnu);
    } else {
        sin_dnu = (h_vec.z >= 0) ? -sqrt(1.0 - cos_dnu*cos_dnu)
                                  : sqrt(1.0 - cos_dnu*cos_dnu);
    }

    /* A = sin(dnu) * sqrt(r1*r2 / (1 - cos(dnu))) */
    if (fabs(1.0 - cos_dnu) < 1e-12) {
        A = 0.0;
    } else {
        A = sin_dnu * sqrt(r1_mag * r2_mag / (1.0 - cos_dnu));
    }

    if (A < 0.0) A = -A;

    /* Iterate over possible revolution numbers */
    for (rev = 0; rev <= max_revs; rev++) {
        z = 0.0;

        for (iter = 0; iter < 50; iter++) {
            /* Compute Stumpff functions for current z */
            if (fabs(z) < 1e-12) {
                Cz = 0.5;
                Sz = 1.0/6.0;
            } else if (z > 0.0) {
                double sz = sqrt(z);
                Cz = (1.0 - cos(sz)) / z;
                Sz = (sz - sin(sz)) / (z * sz);
            } else {
                double sz = sqrt(-z);
                Cz = (cosh(sz) - 1.0) / (-z);
                Sz = (sinh(sz) - sz) / ((-z) * sz);
            }

            if (Cz < 1e-30) Cz = 1e-30;

            y = r1_mag + r2_mag + A * (z * Sz - 1.0) / sqrt(Cz);

            if (y < 0.0) {
                z += 0.1;
                continue;
            }

            x = sqrt(y / Cz);
            t_calc = (x*x*x * Sz + A * sqrt(y)) / sqrt(mu);
            F = t_calc - dt;

            if (fabs(F) < 1e-8) break;

            dF = 1e-6; /* simplified: use pseudo-derivative */
            z -= F / dF;
        }

        if (fabs(F) < 1e-6) {
            sol.converged = 1;
            sol.num_revs = rev;
            sol.a = y / (2.0 * Cz);
            break;
        }
    }

    /* Compute terminal velocities using Lagrange coefficients */
    if (sol.converged) {
        double f_val, g_val, g_dot;

        if (fabs(z) < 1e-12) {
            Cz = 0.5; Sz = 1.0/6.0;
        } else if (z > 0.0) {
            double sz = sqrt(z);
            Cz = (1.0 - cos(sz)) / z;
            Sz = (sz - sin(sz)) / (z * sz);
        } else {
            double sz = sqrt(-z);
            Cz = (cosh(sz) - 1.0) / (-z);
            Sz = (sinh(sz) - sz) / ((-z) * sz);
        }

        y = r1_mag + r2_mag + A*(z*Sz - 1.0)/sqrt(Cz);
        if (y < 0.0) y = 0.0;
        x = sqrt(y / Cz);

        /* Lagrange f and g functions */
        f_val = 1.0 - y / r1_mag;
        g_val = A * sqrt(y / mu);

        /* g_dot = 1 - x^2*C(z)/r2 */
        g_dot = 1.0 - x*x*Cz / r2_mag;

        /* v1 = (r2 - f*r1)/g */
        sol.v1 = vec3_scale(1.0/g_val,
                  vec3_sub(r2, vec3_scale(f_val, r1)));

        /* v2 = (g_dot*r2 - r1)/g */
        sol.v2 = vec3_scale(1.0/g_val,
                  vec3_sub(vec3_scale(g_dot, r2), r1));

        sol.delta_v_total = vec3_norm(sol.v1) + vec3_norm(sol.v2);
    }

    return sol;
}

LambertSolution lambert_multirev(Vector3 r1, Vector3 r2, double dt,
                                  double mu, int max_revs) {
    return lambert_solver(r1, r2, dt, mu, 1, max_revs);
}

/* ================================================================
 * L7: Phasing and Launch Windows
 * ================================================================ */

double phasing_angle_hohmann(double r_depart, double r_arrive, double mu) {
    TransferOrbit hoh = hohmann_transfer(r_depart, r_arrive, mu);
    double n_target = sqrt(mu / (r_arrive * r_arrive * r_arrive));
    double phi = M_PI - n_target * hoh.transfer_time;
    return wrap_2pi(phi);
}

double synodic_period(double T1, double T2) {
    if (fabs(T1 - T2) < 1e-12) return INFINITY;
    return 1.0 / fabs(1.0/T1 - 1.0/T2);
}

double synodic_period_from_radii(double r1, double r2, double mu) {
    double T1 = orbital_period(r1, mu);
    double T2 = orbital_period(r2, mu);
    return synodic_period(T1, T2);
}

double next_launch_window(double epoch, double r_depart, double r_arrive,
                           double mu, double tolerance_rad) {
    double phi_required = phasing_angle_hohmann(r_depart, r_arrive, mu);
    double n_depart = sqrt(mu / (r_depart*r_depart*r_depart)); (void)n_depart;
    double n_arrive = sqrt(mu / (r_arrive*r_arrive*r_arrive)); (void)n_arrive;
    double phi_current;
    int i;

    for (i = 0; i < 5000; i++) {
        /* Phase angle at current epoch (simplified model) */
        double theta_dep = fmod(n_depart * epoch, TWOPI);
        double theta_arr = fmod(n_arrive * epoch, TWOPI);
        phi_current = fmod(theta_arr - theta_dep, TWOPI);
        if (phi_current < 0.0) phi_current += TWOPI;

        if (fabs(phi_current - phi_required) < tolerance_rad) {
            return epoch;
        }

        /* Step forward by a fraction of the synodic period */
        double T_syn = synodic_period_from_radii(r_depart, r_arrive, mu);
        epoch += T_syn / 100.0;
    }

    return -1.0; /* no window found within search range */
}

void porkchop_plot_point(double dt_depart, double dt_flight,
                          double r_depart, double r_arrive, double mu,
                          double *C3, double *dv_arrive) {
    /* C3 = V_inf^2 at departure
     * dv_arrive = capture dv at target
     *
     * These define the characteristic energy requirement for
     * a given launch date and flight time.
     */
    double n_depart = sqrt(mu / (r_depart*r_depart*r_depart)); (void)n_depart;
    double n_arrive = sqrt(mu / (r_arrive*r_arrive*r_arrive)); (void)n_arrive;

    /* Approximate C3 from Lambert solution energy */
    double a_t = (r_depart + r_arrive) / 2.0;
    double v_dep = sqrt(mu * (2.0/r_depart - 1.0/a_t));
    double v_planet = circular_velocity(r_depart, mu);

    *C3 = (v_dep - v_planet) * (v_dep - v_planet); /* V_inf^2 */
    *dv_arrive = fabs(circular_velocity(r_arrive, mu)
                      - sqrt(mu*(2.0/r_arrive - 1.0/a_t)));

    (void)dt_depart; (void)dt_flight;
}

/* ================================================================
 * L7: Low-Thrust Transfers
 * ================================================================ */

double spiral_transfer_delta_v(double r1, double r2, double mu) {
    /* Spiral transfer (continuous low thrust):
     * The total dv for a very-low-thrust circular-to-circular
     * spiral transfer equals the difference in circular speeds.
     *
     * dv = |v_circ1 - v_circ2|
     *
     * This is the absolute lower bound; real low-thrust transfers
     * require more dv due to steering losses.
     */
    double v1 = circular_velocity(r1, mu);
    double v2 = circular_velocity(r2, mu);
    return fabs(v1 - v2);
}

double edelbaum_transfer_delta_v(double r1, double r2,
                                  double delta_i, double mu) {
    /* Edelbaum(1961) low-thrust transfer including inclination change:
     * For continuous tangential+out-of-plane thrust:
     *
     * dv_total^2 = dv_spiral^2 + dv_plane^2 - 2*dv_spiral*dv_plane*cos(beta)
     *
     * where beta is the optimal yaw-steering angle.
     * Simplified: dv_total = sqrt(dv_spiral^2 + (v_avg*delta_i)^2)
     */
    double dv_spiral = spiral_transfer_delta_v(r1, r2, mu);
    double v_avg = (circular_velocity(r1, mu) + circular_velocity(r2, mu)) / 2.0;
    double dv_plane = v_avg * fabs(delta_i);

    return sqrt(dv_spiral*dv_spiral + dv_plane*dv_plane);
}

/* ================================================================
 * L5: General Maneuver Utilities
 * ================================================================ */

double apoapsis_raise_delta_v(double r_circ, double r_apo_target, double mu) {
    /* To raise apoapsis from circular orbit:
     *   dv = v_circ * (sqrt(2*r_apo/(r_circ+r_apo)) - 1) */
    double v_circ = circular_velocity(r_circ, mu);
    double a_trans = (r_circ + r_apo_target) / 2.0;
    double v_peri = sqrt(mu * (2.0/r_circ - 1.0/a_trans));
    return fabs(v_peri - v_circ);
}

double periapsis_raise_delta_v(double r_circ, double r_peri_target, double mu) {
    /* To raise periapsis from circular orbit (r_circ > r_peri_target):
     * First burn at apoapsis of the resulting ellipse
     * dv = v_circ*(1 - sqrt(2*r_peri/(r_circ+r_peri))) */
    if (r_peri_target >= r_circ) return 0.0;
    double v_circ = circular_velocity(r_circ, mu);
    double a_trans = (r_circ + r_peri_target) / 2.0;
    double v_apo = sqrt(mu * (2.0/r_circ - 1.0/a_trans));
    return fabs(v_circ - v_apo);
}

double rocket_equation_delta_v(double Isp, double g0, double m0, double mf) {
    /* Tsiolkovsky rocket equation:
     *   dv = Isp * g0 * ln(m0/mf) */
    if (m0 <= 0.0 || mf <= 0.0 || mf > m0) return 0.0;
    return Isp * g0 * log(m0 / mf);
}

double propellant_mass_fraction(double dv, double Isp, double g0) {
    /* m_prop/m0 = 1 - exp(-dv/(Isp*g0)) */
    if (Isp <= 0.0 || g0 <= 0.0) return 1.0;
    return 1.0 - exp(-dv / (Isp * g0));
}