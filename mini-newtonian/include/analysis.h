/**
 * analysis.h - Trajectory analysis and physics diagnostics (L5/L7)
 *
 * Post-processing tools for trajectory data:
 * phase space analysis, Lyapunov characterization, spectral analysis,
 * conservation law verification, statistical measures.
 *
 * Course: MIT 8.012 + Computational Physics (Thijssen 2007) Ch.8
 */
#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "vec3.h"
#include "integrators.h"

/* ==== L5: Phase space analysis ==== */
/* Compute phase space trajectory (position, momentum) pairs */
typedef struct { double q; double p; } PhasePoint;
void extract_phase_portrait_1d(const Trajectory *traj, int component,
                               double mass, PhasePoint points[]);
double phase_space_area(const PhasePoint points[], int n);

/* ==== L5: Poincare section ==== */
/* Surface of section for detecting chaos in driven systems */
typedef struct { double q; double p; } PoincarePoint;
int compute_poincare_section(const Trajectory *traj, double drive_period,
                             int component, double mass,
                             PoincarePoint section[], int max_points);

/* ==== L5: Lyapunov exponent estimation ==== */
/* Measure sensitivity to initial conditions
 * lambda = (1/t) * ln(|delta(t)| / |delta(0)|)
 * lambda > 0: chaotic; lambda = 0: regular/periodic; lambda < 0: stable fixed point
 */
double estimate_lyapunov_exponent(const Trajectory *traj1,
                                  const Trajectory *traj2);

/* ==== L5: Frequency analysis (FFT-based) ==== */
/* Compute dominant frequency via zero-crossing or FFT */
double dominant_frequency_zero_crossing(const double signal[], int n, double dt);
double dominant_frequency_autocorrelation(const double signal[], int n, double dt);

/* ==== L5: Action-angle variables (1D) ==== */
/* Action: J = (1/2pi) * integral(p dq) over one period */
double compute_action_variable(const PhasePoint points[], int n);
double adiabatic_invariant_check(const PhasePoint points[], int n);

/* ==== L6: Orbital elements from state vectors ==== */
typedef struct {
    double semi_major_axis;   /* a */
    double eccentricity;      /* e */
    double inclination;       /* i (radians) */
    double arg_periapsis;     /* omega (radians) */
    double longitude_ascending; /* Omega (radians) */
    double true_anomaly;      /* f (radians) */
} OrbitalElements;

int orbital_elements_from_state(Vec3 r, Vec3 v, double mu,
                                OrbitalElements *elem);
double orbital_period(double semi_major_axis, double mu);
double orbital_energy_from_elements(double semi_major_axis, double mu);

/* ==== L7: Statistical analysis of trajectories ==== */
double mean_position_offset(const Trajectory *traj, Vec3 reference);
double rms_velocity(const Trajectory *traj);
double max_deviation(const Trajectory *traj, Vec3 reference);
void velocity_autocorrelation(const Trajectory *traj, double lag_times[],
                              double acf[], int n_lags);

/* ==== L7: Collision detection ==== */
int detect_collision_time(const Trajectory *traj1, const Trajectory *traj2,
                          double collision_radius, double *t_collision);
int count_close_approaches(const Trajectory *traj1, const Trajectory *traj2,
                           double threshold, double *times, int max_events);

#endif /* ANALYSIS_H */
