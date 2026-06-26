/**
 * mission.h — Orbital Maneuvers and Interplanetary Transfer
 *
 * Hohmann transfer, bi-elliptic transfer, plane changes,
 * gravity assists (swing-by), patched conics, Lambert problem,
 * phasing, launch windows, and low-thrust approximations.
 *
 * References:
 *   - Vallado, Fundamentals of Astrodynamics, Sec.6.3-6.7
 *   - Bate, Mueller & White, Fundamentals of Astrodynamics, Ch.6-7
 *   - Battin, An Introduction to the Mathematics and Methods of
 *     Astrodynamics, Sec.7.1-7.8
 *   - Prussing & Conway, Orbital Mechanics, Ch.6-8
 *
 * Knowledge: L6 Canonical Systems (Hohmann, Lambert),
 *   L7 Applications (gravity assist, patched conics)
 */

#ifndef MISSION_H
#define MISSION_H

#include "celestial_types.h"

/* ---- L6: Hohmann Transfer ---- */

/** Hohmann transfer between two coplanar circular orbits.
 *  r1 < r2 for outward transfer. Returns TransferOrbit with dv budget. */
TransferOrbit hohmann_transfer(double r1, double r2, double mu);

/** Hohmann transfer inward (r1 > r2). Symmetric to outward case. */
TransferOrbit hohmann_transfer_inward(double r1, double r2, double mu);

/** Optimal Hohmann transfer: finds the transfer orbit and dv costs.
 *  Returns delta_v_total. Also provides dv1, dv2, transfer time. */
void hohmann_transfer_detail(double r1, double r2, double mu,
                              double *dv1, double *dv2, double *dv_total,
                              double *transfer_time);

/* ---- L6: Bi-Elliptic Transfer ---- */

/** Bi-elliptic transfer via intermediate radius r_b.
 *  r_b > max(r1, r2) for potential dv savings over Hohmann
 *  when r2/r1 > ~11.94. */
TransferOrbit bi_elliptic_transfer(double r1, double r2, double r_b, double mu);

/** Compare bi-elliptic vs Hohmann dv costs */
void bi_elliptic_comparison(double r1, double r2, double mu,
                             double *dv_hohmann, double *dv_bielliptic,
                             double *r_b_optimal);

/* ---- L6: Plane Change Maneuvers ---- */

/** Simple plane change at constant speed:
 *  dv = 2*v*sin(delta_i/2) */
double simple_plane_change_delta_v(double v, double delta_i);

/** Plane change at apoapsis (more efficient for large changes) */
double apoapsis_plane_change_delta_v(double v_apo, double delta_i);

/** Combined Hohmann transfer with plane change at apoapsis */
TransferOrbit hohmann_with_plane_change(double r1, double r2,
                                         double delta_i, double mu);

/** One-tangent burn transfer (single-impulse transfer to change orbit) */
void one_tangent_burn_transfer(double r1, double r2, double mu,
                                double *dv, double *e_transfer);

/* ---- L7: Gravity Assist (Swing-By) ---- */

/** Gravity assist turn angle for hyperbolic flyby:
 *  delta = 2*asin(1 / (1 + r_p*V_inf^2/mu)) */
double gravity_assist_turn_angle(double V_inf, double r_p, double mu);

/** Maximum heliocentric dv gain from gravity assist:
 *  dv_max = 2*V_planet*sin(delta/2) */
double gravity_assist_delta_v(double V_planet, double V_inf, double r_p, double mu);

/** Full gravity assist computation with incoming/outgoing vectors */
GravityAssist gravity_assist_full(double V_inf_in, double r_p,
                                   double V_planet, double mu_planet,
                                   double bend_sign);

/** Maximum gravity assist gain (best-case alignment) */
double max_gravity_assist_gain(double V_planet, double V_inf);

/* ---- L7: Patched Conics ---- */

/** Patched conics interplanetary transfer:
 *  Computes total dv for escape + heliocentric transfer + capture. */
TransferOrbit patched_conics_transfer(double r_depart, double r_arrive,
                                       double mu_depart, double mu_arrive,
                                       double r_p_depart, double r_p_arrive,
                                       double mu_sun);

/* ---- L6: Lambert Problem ---- */

/** Lambert problem: find the orbit connecting r1 to r2 in time dt.
 *  Universal variable (Battin) formulation.
 *  Returns LambertSolution with v1, v2, a, e. */
LambertSolution lambert_solver(Vector3 r1, Vector3 r2, double dt,
                                double mu, int prograde, int max_revs);

/** Multi-revolution Lambert solver (for 0 to max_revs complete orbits) */
LambertSolution lambert_multirev(Vector3 r1, Vector3 r2, double dt,
                                  double mu, int max_revs);

/* ---- L7: Phasing and Launch Windows ---- */

/** Phasing angle for Hohmann transfer (target planet must lead by phi):
 *  phi = pi - n_target * T_transfer */
double phasing_angle_hohmann(double r_depart, double r_arrive, double mu);

/** Synodic period: time between same relative configurations
 *  1/T_syn = |1/T1 - 1/T2| */
double synodic_period(double T1, double T2);
double synodic_period_from_radii(double r1, double r2, double mu);

/** Launch window: next epoch when phasing angle is within tolerance */
double next_launch_window(double epoch, double r_depart, double r_arrive,
                           double mu, double tolerance_rad);

/** Porkchop plot data: C3 and dv for range of departure/arrival dates */
void porkchop_plot_point(double dt_depart, double dt_flight,
                          double r_depart, double r_arrive, double mu,
                          double *C3, double *dv_arrive);

/* ---- L7: Low-Thrust Transfers ---- */

/** Spiral transfer dv estimate (low-thrust continuous):
 *  dv_spiral = |v_circ1 - v_circ2| (approximate for very low thrust) */
double spiral_transfer_delta_v(double r1, double r2, double mu);

/** Edelbaum low-thrust plane change + altitude change dv */
double edelbaum_transfer_delta_v(double r1, double r2, double delta_i, double mu);

/* ---- L5: General Maneuver Utilities ---- */

/** Delta-v for raising/lowering apoapsis from circular orbit */
double apoapsis_raise_delta_v(double r_circ, double r_apo_target, double mu);

/** Delta-v for raising/lowering periapsis from circular orbit */
double periapsis_raise_delta_v(double r_circ, double r_peri_target, double mu);

/** Total mission delta-v from ideal rocket equation:
 *  dv = Isp * g0 * ln(m0/mf) */
double rocket_equation_delta_v(double Isp, double g0, double m0, double mf);

/** Propellant mass fraction from delta-v:
 *  m_prop/m0 = 1 - exp(-dv/(Isp*g0)) */
double propellant_mass_fraction(double dv, double Isp, double g0);

#endif /* MISSION_H */