/**
 * integrators.h - ODE numerical integrators for Newtonian mechanics
 * (L5 Computational Methods / L6 Canonical Systems)
 *
 * Solves d^2r/dt^2 = a(r,v,t) = F(r,v,t)/m
 * using explicit, symplectic, and adaptive methods.
 *
 * Course mapping:
 *   MIT 8.012 Ch.7 (Numerical methods)
 *   Hairer, Lubich, Wanner (2006) Geometric Numerical Integration
 *   Butcher (2008) Numerical Methods for ODEs
 *
 * Key theorems:
 *   - Euler: O(dt) local truncation error, non-symplectic
 *   - RK4: O(dt^4) local truncation error, non-symplectic
 *   - Velocity Verlet: O(dt^2) global error, SYMPLECTIC (preserves phase space)
 *   - Leapfrog: O(dt^2), symplectic, offset velocity and position grids
 *   - RK45 (Dormand-Prince): adaptive 5(4) embedded pair, error control
 */

#ifndef INTEGRATORS_H
#define INTEGRATORS_H

#include "vec3.h"

/* ===== L1: Particle state for integration ===== */
typedef struct {
    double t;
    Vec3  r;
    Vec3  v;
} ParticleState;

ParticleState particle_state_make(double t, Vec3 r, Vec3 v);

/* ===== L5: First-order ODE integrators (dy/dt = f(t,y)) ===== */
/* Signature: void (*step_fn)(f, t, y, n, dt, y_out) */
/* n = dimension of y vector */

void euler_step(void (*f)(double t, const double y[], double dydt[]),
                double t, const double y[], int n, double dt, double y_out[]);

void rk2_midpoint_step(void (*f)(double t, const double y[], double dydt[]),
                       double t, const double y[], int n, double dt, double y_out[]);

void rk4_classical_step(void (*f)(double t, const double y[], double dydt[]),
                        double t, const double y[], int n, double dt, double y_out[]);

void rk4_step_with_stages(void (*f)(double t, const double y[], double dydt[]),
                          double t, const double y[], int n, double dt,
                          double y_out[], double k1[], double k2[],
                          double k3[], double k4[]);

/* ===== L5: Adaptive RK45 (Dormand-Prince 5(4) embedded pair) ===== */
/* Returns: 0=accepted, 1=rejected; err_out = error estimate; dt_sug = suggested next dt */
int rk45_dormand_prince_step(
    void (*f)(double t, const double y[], double dydt[]),
    double t, const double y[], int n, double dt, double tol,
    double y_out[], double *err_out, double *dt_suggested);

/* ===== L5: Second-order ODE integrators for d^2r/dt^2 = a(r,v,t) ===== */
/* Each step updates particle state in-place or returns new state */

ParticleState euler_cromer_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt);

ParticleState velocity_verlet_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt);

/* Leapfrog: staggered velocity updates */
void leapfrog_init(Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
                   Vec3 r, Vec3 v, double t, double dt,
                   Vec3 *r_out, Vec3 *v_half_out);

void leapfrog_step(Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
                   Vec3 r, Vec3 v_half, double t, double dt,
                   Vec3 *r_new, Vec3 *v_half_new);

ParticleState rk4_second_order_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt);

/* ===== L5: Symplectic integrators (preserve phase space structure) ===== */
/* Knowledge points:
 *   - Separable Hamiltonian H = T(p) + V(q)
 *   - Symplectic Euler: explicit in p, implicit in q (or vice versa)
 *   - Yoshida 4th-order: composition of leapfrog steps
 */
ParticleState symplectic_euler_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt);

ParticleState yoshida4_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt);

/* ===== L5: General solver loops ===== */
typedef struct {
    double *ts;
    Vec3   *positions;
    Vec3   *velocities;
    int     n_points;
    int     capacity;
} Trajectory;

Trajectory *trajectory_create(int max_points);
void trajectory_destroy(Trajectory *traj);
int trajectory_record(Trajectory *traj, double t, Vec3 r, Vec3 v);

/* Fixed-step solve for a single particle */
Trajectory *solve_fixed_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    Vec3 r0, Vec3 v0, double t_end, double dt,
    int method,  /* 0=Euler-Cromer, 1=Verlet, 2=RK4 */
    int record_every);

/* Adaptive-step solve using RK45 */
Trajectory *solve_adaptive(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    Vec3 r0, Vec3 v0, double t_end,
    double dt_init, double tol, double dt_min, double dt_max,
    int *n_accepted, int *n_rejected);

/* ===== L5: N-body Verlet integrator ===== */
/* Knowledge point: symplectic integration for N-body gravitational systems */
typedef struct {
    double  t;
    double *masses;
    Vec3   *positions;
    Vec3   *velocities;
    int     n_bodies;
} NBodySystem;

NBodySystem *nbody_system_create(int n);
void nbody_system_destroy(NBodySystem *sys);
void nbody_verlet_step(
    void (*accel_func)(const NBodySystem *sys, Vec3 accel[]),
    NBodySystem *sys, double dt);

/* ===== L5: Backward Euler (implicit) for stiff systems ===== */
/* Knowledge point: implicit integration for stiff ODEs (e.g., large damping) */
ParticleState backward_euler_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt, int max_iter, double tol);

#endif /* INTEGRATORS_H */
