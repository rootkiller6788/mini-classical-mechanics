/**
 * analysis.c - Trajectory analysis and physics diagnostics
 * Each function implements an independent analysis knowledge point.
 *
 * References:
 *   Goldstein Ch.10 (action-angle variables)
 *   Strogatz (1994) Nonlinear Dynamics and Chaos
 *   Thijssen (2007) Computational Physics Ch.8
 *   Murray & Dermott (2000) Solar System Dynamics
 */
#include "analysis.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 * L5: Phase space analysis
 * ================================================================ */

void extract_phase_portrait_1d(const Trajectory *traj, int component,
                               double mass, PhasePoint points[]) {
    for (int i = 0; i < traj->n_points; i++) {
        double q, v;
        switch (component) {
            case 0: q = traj->positions[i].x; v = traj->velocities[i].x; break;
            case 1: q = traj->positions[i].y; v = traj->velocities[i].y; break;
            default: q = traj->positions[i].z; v = traj->velocities[i].z; break;
        }
        points[i].q = q;
        points[i].p = mass * v;
    }
}

/* Area of closed curve in phase space via trapezoidal rule
 * A = (1/2) * |sum(q_i * Delta(p_i) - p_i * Delta(q_i))|
 */
double phase_space_area(const PhasePoint points[], int n) {
    if (n < 3) return 0.0;
    double area = 0.0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        area += points[i].q * points[j].p - points[j].q * points[i].p;
    }
    return 0.5 * fabs(area);
}

/* ================================================================
 * L5: Poincare section
 *
 * Knowledge point: For driven systems with period T_drive,
 * recording the state at t = n*T_drive produces a Poincare map.
 * In a regular (integrable) system, points lie on invariant tori.
 * In a chaotic system, points scatter stochastically.
 * ================================================================ */

int compute_poincare_section(const Trajectory *traj, double drive_period,
                             int component, double mass,
                             PoincarePoint section[], int max_points) {
    int count = 0;
    for (int i = 0; i < traj->n_points && count < max_points; i++) {
        double phase = fmod(traj->ts[i], drive_period);
        if (phase < 1e-6 || fabs(phase - drive_period) < 1e-6) {
            double q, v;
            switch (component) {
                case 0: q = traj->positions[i].x; v = traj->velocities[i].x; break;
                case 1: q = traj->positions[i].y; v = traj->velocities[i].y; break;
                default: q = traj->positions[i].z; v = traj->velocities[i].z; break;
            }
            section[count].q = q;
            section[count].p = mass * v;
            count++;
        }
    }
    return count;
}

/* ================================================================
 * L5: Lyapunov exponent estimation
 *
 * Knowledge point: The largest Lyapunov exponent measures the
 * exponential divergence rate of nearby trajectories.
 *
 * lambda = lim_{t->inf} (1/t) * ln(|delta(t)| / |delta(0)|)
 *
 * Positive lambda: chaos (exponential sensitivity)
 * Zero lambda: regular motion (linear divergence at most)
 * Negative lambda: stable fixed point (convergence)
 *
 * Course: Strogatz Ch.10, Goldstein Ch.11.7
 * ================================================================ */

double estimate_lyapunov_exponent(const Trajectory *traj1,
                                  const Trajectory *traj2) {
    int n = (traj1->n_points < traj2->n_points) ? traj1->n_points : traj2->n_points;
    if (n < 2) return 0.0;

    /* Initial separation */
    Vec3 d0 = vec3_sub(traj1->positions[0], traj2->positions[0]);
    double dist0 = vec3_norm(d0);
    if (dist0 < 1e-15) return 0.0;

    /* Log-averaged divergence rate */
    double sum_lambda = 0.0;
    int valid_points = 0;
    for (int i = 1; i < n; i++) {
        double t = traj1->ts[i];
        if (t < 1e-10) continue;
        Vec3 di = vec3_sub(traj1->positions[i], traj2->positions[i]);
        double dist = vec3_norm(di);
        if (dist < 1e-15) continue;
        sum_lambda += log(dist / dist0) / t;
        valid_points++;
    }
    if (valid_points == 0) return 0.0;
    return sum_lambda / valid_points;
}

/* ================================================================
 * L5: Frequency analysis
 * ================================================================ */

/* Zero-crossing: count sign changes, divide by 2 for periods */
double dominant_frequency_zero_crossing(const double signal[], int n, double dt) {
    if (n < 3 || dt < 1e-15) return 0.0;
    int crossings = 0;
    for (int i = 1; i < n; i++) {
        if (signal[i-1] * signal[i] < 0.0) crossings++;
    }
    double period = (n - 1) * dt / (0.5 * crossings + 0.5);
    return 1.0 / period;
}

/* Autocorrelation method: first peak in ACF gives dominant period */
double dominant_frequency_autocorrelation(const double signal[], int n, double dt) {
    if (n < 10 || dt < 1e-15) return 0.0;
    /* Compute mean */
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += signal[i];
    mean /= n;
    /* Compute ACF for lags up to n/2 */
    int max_lag = n / 2;
    double *acf = (double *)malloc(max_lag * sizeof(double));
    for (int lag = 0; lag < max_lag; lag++) {
        double num = 0.0, den = 0.0;
        for (int i = 0; i < n - lag; i++) {
            double da = signal[i] - mean;
            double db = signal[i + lag] - mean;
            num += da * db;
            if (lag == 0) den += da * da;
        }
        acf[lag] = (den > 1e-15) ? num / den : 0.0;
    }
    /* Find first peak after zero crossing */
    int peak_lag = 1;
    for (int lag = 2; lag < max_lag - 1; lag++) {
        if (acf[lag] > acf[lag-1] && acf[lag] > acf[lag+1] && acf[lag] > 0.3) {
            peak_lag = lag;
            break;
        }
    }
    free(acf);
    if (peak_lag < 1) peak_lag = 1;
    return 1.0 / (peak_lag * dt);
}

/* ================================================================
 * L5: Action-angle variables
 * ================================================================ */

/* Action: J = (1/2pi) * integral(p dq) for 1D periodic motion
 * Numerically: sum over closed phase-space orbit */
double compute_action_variable(const PhasePoint points[], int n) {
    /* Trapezoidal integration: integral p dq */
    double action = 0.0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        double dq = points[j].q - points[i].q;
        double p_avg = 0.5 * (points[i].p + points[j].p);
        action += p_avg * dq;
    }
    return fabs(action) / (2.0 * M_PI);
}

/* Check if J is conserved (adiabatic invariant) */
double adiabatic_invariant_check(const PhasePoint points[], int n) {
    double J0 = compute_action_variable(points, n);
    /* Normalized deviation of J */
    double sum_dev = 0.0;
    for (int i = 0; i < n; i++) {
        double J_i = points[i].q * points[i].p;  /* simplified: local p*q */
        sum_dev += fabs(J_i - J0) / (fabs(J0) + 1e-15);
    }
    return sum_dev / n;
}

/* ================================================================
 * L6: Orbital elements from state vectors
 *
 * Knowledge point: Keplerian orbital elements uniquely define a
 * conic section orbit. Conversion from Cartesian state vectors
 * (r, v) to Kepler elements is fundamental to celestial mechanics.
 *
 * Algorithm: compute specific angular momentum h = r x v,
 * eccentricity vector e = (v x h)/mu - r_hat,
 * from which all orbital elements follow.
 *
 * Course: Goldstein Ch.3.7, Murray & Dermott Ch.2
 * ================================================================ */

int orbital_elements_from_state(Vec3 r, Vec3 v, double mu,
                                OrbitalElements *elem) {
    Vec3 h = vec3_cross(r, v);
    double h_mag = vec3_norm(h);
    if (h_mag < 1e-15) return 0;  /* radial infall */

    Vec3 r_hat = vec3_normalize(r);
    double r_mag = vec3_norm(r);
    double v2 = vec3_norm2(v);

    /* Eccentricity vector: e = (v x h)/mu - r_hat */
    Vec3 e_vec = vec3_sub(vec3_div(vec3_cross(v, h), mu), r_hat);
    double e = vec3_norm(e_vec);

    /* Semi-major axis from vis-viva: v^2/2 - mu/r = -mu/(2a) */
    double specific_energy = 0.5 * v2 - mu / r_mag;
    double a;
    if (fabs(specific_energy) < 1e-15) {
        a = 1e300;  /* parabolic */
    } else {
        a = -mu / (2.0 * specific_energy);
    }

    /* Inclination */
    double i = acos(h.z / h_mag);

    /* Node vector: n = (0,0,1) x h */
    Vec3 z_hat = vec3_make(0, 0, 1);
    Vec3 n_vec = vec3_cross(z_hat, h);
    double n_mag = vec3_norm(n_vec);

    /* Longitude of ascending node */
    double Omega;
    if (n_mag > 1e-15) {
        Omega = acos(n_vec.x / n_mag);
        if (n_vec.y < 0) Omega = 2.0 * M_PI - Omega;
    } else {
        Omega = 0.0;
    }

    /* Argument of periapsis */
    double omega;
    if (e > 1e-10 && n_mag > 1e-15) {
        omega = acos(vec3_dot(n_vec, e_vec) / (n_mag * e));
        if (e_vec.z < 0) omega = 2.0 * M_PI - omega;
    } else if (e > 1e-10) {
        omega = acos(e_vec.x / e);
        if (e_vec.y < 0) omega = 2.0 * M_PI - omega;
    } else {
        omega = 0.0;
    }

    /* True anomaly */
    double f;
    if (e > 1e-10) {
        f = acos(vec3_dot(e_vec, r) / (e * r_mag));
        if (vec3_dot(r, v) < 0) f = 2.0 * M_PI - f;
    } else {
        /* Circular: use argument of latitude */
        if (n_mag > 1e-15) {
            f = acos(vec3_dot(n_vec, r) / (n_mag * r_mag));
            if (r.z < 0) f = 2.0 * M_PI - f;
        } else {
            f = acos(r.x / r_mag);
            if (r.y < 0) f = 2.0 * M_PI - f;
        }
    }

    elem->semi_major_axis = a;
    elem->eccentricity = e;
    elem->inclination = i;
    elem->arg_periapsis = omega;
    elem->longitude_ascending = Omega;
    elem->true_anomaly = f;
    return 1;
}

double orbital_period(double semi_major_axis, double mu) {
    if (semi_major_axis <= 0.0 || mu < 1e-15) return 0.0;
    return 2.0 * M_PI * sqrt(semi_major_axis * semi_major_axis * semi_major_axis / mu);
}

double orbital_energy_from_elements(double semi_major_axis, double mu) {
    if (fabs(semi_major_axis) < 1e-15) return 0.0;
    return -mu / (2.0 * semi_major_axis);
}

/* ================================================================
 * L7: Statistical analysis of trajectories
 * ================================================================ */

double mean_position_offset(const Trajectory *traj, Vec3 reference) {
    double sum = 0.0;
    for (int i = 0; i < traj->n_points; i++) {
        sum += vec3_norm(vec3_sub(traj->positions[i], reference));
    }
    return sum / traj->n_points;
}

double rms_velocity(const Trajectory *traj) {
    double sum_v2 = 0.0;
    for (int i = 0; i < traj->n_points; i++) {
        sum_v2 += vec3_norm2(traj->velocities[i]);
    }
    return sqrt(sum_v2 / traj->n_points);
}

double max_deviation(const Trajectory *traj, Vec3 reference) {
    double max_d = 0.0;
    for (int i = 0; i < traj->n_points; i++) {
        double d = vec3_norm(vec3_sub(traj->positions[i], reference));
        if (d > max_d) max_d = d;
    }
    return max_d;
}

/* Velocity autocorrelation function: C(tau) = <v(0).v(tau)> / <v^2>
 * Knowledge point: ACF decay time characterizes memory of the system.
 * Fast decay: chaotic/diffusive; slow/oscillatory decay: regular.
 */
void velocity_autocorrelation(const Trajectory *traj, double lag_times[],
                              double acf[], int n_lags) {
    int n = traj->n_points;
    if (n < 2) return;
    /* Mean square velocity */
    double v2_mean = 0.0;
    for (int i = 0; i < n; i++) v2_mean += vec3_norm2(traj->velocities[i]);
    v2_mean /= n;
    if (v2_mean < 1e-15) return;
    /* Compute ACF for each lag */
    for (int k = 0; k < n_lags && k < n; k++) {
        double sum = 0.0;
        int count = 0;
        for (int i = 0; i < n - k; i++) {
            sum += vec3_dot(traj->velocities[i], traj->velocities[i + k]);
            count++;
        }
        acf[k] = (count > 0) ? sum / (count * v2_mean) : 0.0;
        lag_times[k] = k * (traj->ts[n-1] - traj->ts[0]) / (n - 1);
    }
}

/* ================================================================
 * L7: Collision detection
 * ================================================================ */

/* Detect when two particles come within collision_radius */
int detect_collision_time(const Trajectory *traj1, const Trajectory *traj2,
                          double collision_radius, double *t_collision) {
    int n = (traj1->n_points < traj2->n_points) ? traj1->n_points : traj2->n_points;
    for (int i = 0; i < n; i++) {
        double d = vec3_norm(vec3_sub(traj1->positions[i], traj2->positions[i]));
        if (d <= collision_radius) {
            *t_collision = traj1->ts[i];
            return 1;
        }
    }
    return 0;
}

/* Count close approaches below threshold distance */
int count_close_approaches(const Trajectory *traj1, const Trajectory *traj2,
                           double threshold, double *times, int max_events) {
    int count = 0;
    int n = (traj1->n_points < traj2->n_points) ? traj1->n_points : traj2->n_points;
    for (int i = 0; i < n && count < max_events; i++) {
        double d = vec3_norm(vec3_sub(traj1->positions[i], traj2->positions[i]));
        if (d < threshold) {
            times[count] = traj1->ts[i];
            count++;
            /* Skip subsequent points in same event */
            while (i + 1 < n &&
                   vec3_norm(vec3_sub(traj1->positions[i+1], traj2->positions[i+1])) < threshold)
                i++;
        }
    }
    return count;
}
