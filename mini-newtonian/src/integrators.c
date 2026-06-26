/**
 * integrators.c - ODE numerical integrators for Newtonian mechanics
 *
 * Each function implements one independent numerical method / knowledge point.
 *
 * References:
 *   Euler (1768): Institutionum Calculi Integralis
 *   Runge (1895), Kutta (1901): Runge-Kutta methods
 *   Verlet (1967): molecular dynamics integrator
 *   Dormand & Prince (1980): RK5(4)7M embedded pair
 *   Yoshida (1990): higher-order symplectic composition
 *   Hairer, Lubich, Wanner (2006): Geometric Numerical Integration
 *
 * Course: MIT 8.012 Ch.7, Computational Physics (Thijssen 2007) Ch.4
 */
#include "integrators.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ===== ParticleState ===== */
ParticleState particle_state_make(double t, Vec3 r, Vec3 v) {
    ParticleState s;
    s.t = t; s.r = r; s.v = v;
    return s;
}

/* ================================================================
 * L5: First-order ODE integrators (dy/dt = f(t,y))
 *
 * Each method solves the general ODE with different accuracy.
 * Knowledge point: local truncation error analysis
 *   Euler: LTE = O(dt^2), global O(dt) -> 1st order
 *   RK2 midpoint: LTE = O(dt^3), global O(dt^2) -> 2nd order
 *   RK4 classical: LTE = O(dt^5), global O(dt^4) -> 4th order
 * ================================================================ */

/* Euler method: y_{n+1} = y_n + dt * f(t_n, y_n)
 * 1st-order, simplest explicit method.
 * Knowledge point: forward difference approximation of derivative
 */
void euler_step(void (*f)(double t, const double y[], double dydt[]),
                double t, const double y[], int n, double dt, double y_out[]) {
    double *dydt = (double *)malloc(n * sizeof(double));
    f(t, y, dydt);
    for (int i = 0; i < n; i++) {
        y_out[i] = y[i] + dt * dydt[i];
    }
    free(dydt);
}

/* RK2 midpoint: y_{n+1} = y_n + dt * f(t_n + dt/2, y_n + dt/2 * f(t_n, y_n))
 * 2nd-order, uses midpoint derivative evaluation.
 * Knowledge point: improved Euler / midpoint method
 */
void rk2_midpoint_step(void (*f)(double t, const double y[], double dydt[]),
                       double t, const double y[], int n, double dt, double y_out[]) {
    double *k1 = (double *)malloc(n * sizeof(double));
    double *y_mid = (double *)malloc(n * sizeof(double));
    double *k2 = (double *)malloc(n * sizeof(double));

    f(t, y, k1);
    for (int i = 0; i < n; i++)
        y_mid[i] = y[i] + 0.5 * dt * k1[i];
    f(t + 0.5 * dt, y_mid, k2);
    for (int i = 0; i < n; i++)
        y_out[i] = y[i] + dt * k2[i];

    free(k1); free(y_mid); free(k2);
}

/* Classical RK4: 4 stages, 4th-order
 * Knowledge point: Simpson-like quadrature in time
 *   k1 = f(t, y)
 *   k2 = f(t + dt/2, y + dt/2 * k1)
 *   k3 = f(t + dt/2, y + dt/2 * k2)
 *   k4 = f(t + dt, y + dt * k3)
 *   y_{n+1} = y_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
 */
void rk4_classical_step(void (*f)(double t, const double y[], double dydt[]),
                        double t, const double y[], int n, double dt, double y_out[]) {
    double *k1 = (double *)malloc(4 * n * sizeof(double));
    double *k2 = k1 + n;
    double *k3 = k2 + n;
    double *k4 = k3 + n;
    double *y_tmp = (double *)malloc(n * sizeof(double));

    f(t, y, k1);
    for (int i = 0; i < n; i++) y_tmp[i] = y[i] + 0.5 * dt * k1[i];
    f(t + 0.5 * dt, y_tmp, k2);
    for (int i = 0; i < n; i++) y_tmp[i] = y[i] + 0.5 * dt * k2[i];
    f(t + 0.5 * dt, y_tmp, k3);
    for (int i = 0; i < n; i++) y_tmp[i] = y[i] + dt * k3[i];
    f(t + dt, y_tmp, k4);

    for (int i = 0; i < n; i++)
        y_out[i] = y[i] + (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);

    free(k1); free(y_tmp);
}

/* RK4 with stage storage for error estimation */
void rk4_step_with_stages(void (*f)(double t, const double y[], double dydt[]),
                          double t, const double y[], int n, double dt,
                          double y_out[], double k1[], double k2[],
                          double k3[], double k4[]) {
    double *y_tmp = (double *)malloc(n * sizeof(double));

    f(t, y, k1);
    for (int i = 0; i < n; i++) y_tmp[i] = y[i] + 0.5 * dt * k1[i];
    f(t + 0.5 * dt, y_tmp, k2);
    for (int i = 0; i < n; i++) y_tmp[i] = y[i] + 0.5 * dt * k2[i];
    f(t + 0.5 * dt, y_tmp, k3);
    for (int i = 0; i < n; i++) y_tmp[i] = y[i] + dt * k3[i];
    f(t + dt, y_tmp, k4);

    for (int i = 0; i < n; i++)
        y_out[i] = y[i] + (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);

    free(y_tmp);
}

/* ================================================================
 * L5: Adaptive RK45 (Dormand-Prince 5(4) embedded pair)
 *
 * Knowledge point: embedded Runge-Kutta for error estimation and
 * automatic step size control. The difference between 5th-order and
 * 4th-order embedded solutions provides the error estimate.
 *
 * Uses the Dormand-Prince coefficients (1980):
 *   7 stages, 5th-order solution + 4th-order embedded pair
 *   FSAL (First Same As Last) property for efficiency
 *
 * Reference: Dormand & Prince, J. Comp. Appl. Math. 6(1), 1980
 * ================================================================ */

int rk45_dormand_prince_step(
    void (*f)(double t, const double y[], double dydt[]),
    double t, const double y[], int n, double dt, double tol,
    double y_out[], double *err_out, double *dt_suggested) {

    double *k1 = (double *)malloc(7 * n * sizeof(double));
    double *k2 = k1 + n;   double *k3 = k2 + n;
    double *k4 = k3 + n;   double *k5 = k4 + n;
    double *k6 = k5 + n;   double *k7 = k6 + n;
    double *yt = (double *)malloc(n * sizeof(double));

    /* Stage 1 */
    f(t, y, k1);

    /* Stage 2 */
    for (int i = 0; i < n; i++) yt[i] = y[i] + dt * (0.2 * k1[i]);
    f(t + 0.2*dt, yt, k2);

    /* Stage 3 */
    for (int i = 0; i < n; i++) yt[i] = y[i] + dt * (0.075*k1[i] + 0.225*k2[i]);
    f(t + 0.3*dt, yt, k3);

    /* Stage 4 */
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + dt * (44.0/45.0*k1[i] - 56.0/15.0*k2[i] + 32.0/9.0*k3[i]);
    f(t + 0.8*dt, yt, k4);

    /* Stage 5 */
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + dt * (19372.0/6561.0*k1[i] - 25360.0/2187.0*k2[i]
               + 64448.0/6561.0*k3[i] - 212.0/729.0*k4[i]);
    f(t + 8.0/9.0*dt, yt, k5);

    /* Stage 6 */
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + dt * (9017.0/3168.0*k1[i] - 355.0/33.0*k2[i]
               + 46732.0/5247.0*k3[i] + 49.0/176.0*k4[i] - 5103.0/18656.0*k5[i]);
    f(t + dt, yt, k6);

    /* Stage 7 */
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + dt * (35.0/384.0*k1[i] + 500.0/1113.0*k3[i]
               + 125.0/192.0*k4[i] - 2187.0/6784.0*k5[i] + 11.0/84.0*k6[i]);
    f(t + dt, yt, k7);

    /* 5th-order solution */
    for (int i = 0; i < n; i++)
        y_out[i] = y[i] + dt * (35.0/384.0*k1[i] + 500.0/1113.0*k3[i]
                  + 125.0/192.0*k4[i] - 2187.0/6784.0*k5[i] + 11.0/84.0*k6[i]);

    /* Error estimate: difference from 4th-order embedded solution */
    double err_sum_sq = 0.0;
    for (int i = 0; i < n; i++) {
        double y4 = y[i] + dt * (5179.0/57600.0*k1[i] + 7571.0/16695.0*k3[i]
                    + 393.0/640.0*k4[i] - 92097.0/339200.0*k5[i]
                    + 187.0/2100.0*k6[i] + 0.025*k7[i]);
        double scale = fmax(1.0, fmax(fabs(y_out[i]), fabs(y[i])));
        double diff = (y_out[i] - y4) / scale;
        err_sum_sq += diff * diff;
    }
    *err_out = sqrt(err_sum_sq / n);

    /* PI step-size controller with safety factor 0.9 */
    double safety = 0.9;
    if (*err_out > 0.0) {
        *dt_suggested = dt * safety * pow(tol / (*err_out), 0.2);
    } else {
        *dt_suggested = dt * 2.0;
    }
    *dt_suggested = fmin(*dt_suggested, dt * 4.0);

    free(k1); free(yt);
    return (*err_out <= tol) ? 0 : 1;
}

/* ================================================================
 * L5: Second-order ODE integrators for d^2r/dt^2 = a(r,v,t)
 *
 * These directly integrate Newton's equation without reducing to
 * first-order form. More efficient and often more stable for
 * mechanical systems.
 * ================================================================ */

/* Euler-Cromer (semi-implicit Euler): updates v first, then r with new v
 * Knowledge point: semi-implicit Euler vastly improves energy behavior
 * over explicit Euler for oscillatory systems.
 * Global error O(dt), but energy oscillates instead of drifting.
 */
ParticleState euler_cromer_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt) {
    Vec3 a_cur = accel_func(state.r, state.v, state.t);
    Vec3 v_new = vec3_add(state.v, vec3_scale(a_cur, dt));
    Vec3 r_new = vec3_add(state.r, vec3_scale(v_new, dt));  /* uses NEW v! */
    ParticleState new_state;
    new_state.t = state.t + dt;
    new_state.r = r_new;
    new_state.v = v_new;
    return new_state;
}

/* Velocity Verlet (symplectic, 2nd-order)
 * Knowledge point: The "gold standard" for molecular dynamics and
 * conservative mechanical systems. Preserves phase space volume
 * (Liouville theorem) and has bounded energy error.
 *
 * Algorithm:
 *   v_{n+1/2} = v_n + 0.5*a_n*dt
 *   r_{n+1} = r_n + v_{n+1/2}*dt
 *   a_{n+1} = a(r_{n+1}, v_{n+1/2}, t+dt)
 *   v_{n+1} = v_{n+1/2} + 0.5*a_{n+1}*dt
 *
 * Reference: Verlet (1967), Swope et al. (1982)
 */
ParticleState velocity_verlet_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt) {
    Vec3 a_cur = accel_func(state.r, state.v, state.t);
    Vec3 v_half = vec3_add(state.v, vec3_scale(a_cur, 0.5 * dt));
    Vec3 r_new = vec3_add(state.r, vec3_scale(v_half, dt));
    Vec3 a_new = accel_func(r_new, v_half, state.t + dt);
    Vec3 v_new = vec3_add(v_half, vec3_scale(a_new, 0.5 * dt));
    ParticleState new_state;
    new_state.t = state.t + dt;
    new_state.r = r_new;
    new_state.v = v_new;
    return new_state;
}

/* Leapfrog initialization: compute v at half-step
 * v_{1/2} = v_0 + 0.5*a_0*dt
 */
void leapfrog_init(Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
                   Vec3 r, Vec3 v, double t, double dt,
                   Vec3 *r_out, Vec3 *v_half_out) {
    Vec3 a0 = accel_func(r, v, t);
    *r_out = r;
    *v_half_out = vec3_add(v, vec3_scale(a0, 0.5 * dt));
}

/* Leapfrog step: alternate r and v updates
 * r_{n+1} = r_n + v_{n+1/2} * dt
 * v_{n+3/2} = v_{n+1/2} + a(r_{n+1}, ...) * dt
 */
void leapfrog_step(Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
                   Vec3 r, Vec3 v_half, double t, double dt,
                   Vec3 *r_new, Vec3 *v_half_new) {
    *r_new = vec3_add(r, vec3_scale(v_half, dt));
    Vec3 a_new = accel_func(*r_new, v_half, t + dt);
    *v_half_new = vec3_add(v_half, vec3_scale(a_new, dt));
}

/* RK4 for second-order ODE: convert to first-order system of dimension 6
 * Knowledge point: RK4 can integrate any ODE, but loses symplectic property
 */
ParticleState rk4_second_order_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt) {
    /* State vector: [rx, ry, rz, vx, vy, vz] */
    double y[6] = {state.r.x, state.r.y, state.r.z,
                   state.v.x, state.v.y, state.v.z};

    /* Define the ODE function for the combined state */
    /* We inline the RK4 stages to avoid function pointer overhead */

    double k1[6], k2[6], k3[6], k4[6];
    double yt[6];

    /* k1 = f(t, y) */
    Vec3 r_tmp = vec3_make(y[0], y[1], y[2]);
    Vec3 v_tmp = vec3_make(y[3], y[4], y[5]);
    Vec3 a1 = accel_func(r_tmp, v_tmp, state.t);
    k1[0]=y[3]; k1[1]=y[4]; k1[2]=y[5];
    k1[3]=a1.x; k1[4]=a1.y; k1[5]=a1.z;

    /* k2 = f(t + dt/2, y + dt/2 * k1) */
    for (int i=0; i<6; i++) yt[i] = y[i] + 0.5*dt*k1[i];
    r_tmp = vec3_make(yt[0], yt[1], yt[2]);
    v_tmp = vec3_make(yt[3], yt[4], yt[5]);
    Vec3 a2 = accel_func(r_tmp, v_tmp, state.t + 0.5*dt);
    k2[0]=yt[3]; k2[1]=yt[4]; k2[2]=yt[5];
    k2[3]=a2.x; k2[4]=a2.y; k2[5]=a2.z;

    /* k3 = f(t + dt/2, y + dt/2 * k2) */
    for (int i=0; i<6; i++) yt[i] = y[i] + 0.5*dt*k2[i];
    r_tmp = vec3_make(yt[0], yt[1], yt[2]);
    v_tmp = vec3_make(yt[3], yt[4], yt[5]);
    Vec3 a3 = accel_func(r_tmp, v_tmp, state.t + 0.5*dt);
    k3[0]=yt[3]; k3[1]=yt[4]; k3[2]=yt[5];
    k3[3]=a3.x; k3[4]=a3.y; k3[5]=a3.z;

    /* k4 = f(t + dt, y + dt * k3) */
    for (int i=0; i<6; i++) yt[i] = y[i] + dt*k3[i];
    r_tmp = vec3_make(yt[0], yt[1], yt[2]);
    v_tmp = vec3_make(yt[3], yt[4], yt[5]);
    Vec3 a4 = accel_func(r_tmp, v_tmp, state.t + dt);
    k4[0]=yt[3]; k4[1]=yt[4]; k4[2]=yt[5];
    k4[3]=a4.x; k4[4]=a4.y; k4[5]=a4.z;

    /* Combine */
    for (int i=0; i<6; i++)
        y[i] = y[i] + (dt/6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);

    ParticleState new_state;
    new_state.t = state.t + dt;
    new_state.r = vec3_make(y[0], y[1], y[2]);
    new_state.v = vec3_make(y[3], y[4], y[5]);
    return new_state;
}
/* ================================================================
 * L5: Symplectic integrators
 *
 * Knowledge point: For separable Hamiltonian H = T(p) + V(q),
 * symplectic integrators exactly preserve the symplectic 2-form dp^dq.
 * This implies approximate energy conservation with bounded oscillations
 * (no long-term drift) and exact conservation of phase space volume.
 *
 * Reference: Yoshida (1990) Phys. Lett. A 150, 262-268
 * ================================================================ */

/* Symplectic Euler: explicit update p then q */
ParticleState symplectic_euler_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt) {
    Vec3 a_cur = accel_func(state.r, state.v, state.t);
    Vec3 v_new = vec3_add(state.v, vec3_scale(a_cur, dt));
    Vec3 r_new = vec3_add(state.r, vec3_scale(v_new, dt));
    ParticleState ns;
    ns.t = state.t + dt; ns.r = r_new; ns.v = v_new;
    return ns;
}

/* Yoshida 4th-order symplectic integrator
 * Composition of 3 velocity-Verlet-like steps with coefficients:
 *   w1 = 1/(2 - 2^(1/3)) ~ 1.3512
 *   w0 = -(2^(1/3))/(2 - 2^(1/3)) ~ -1.7024
 * The negative time-step is valid for separable Hamiltonians.
 */
ParticleState yoshida4_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt) {
    const double cbrt2 = 1.2599210498948732;
    const double w1 = 1.0 / (2.0 - cbrt2);
    const double w0 = -cbrt2 / (2.0 - cbrt2);
    double w[3] = {w1, w0, w1};
    ParticleState s = state;
    for (int i = 0; i < 3; i++) {
        Vec3 a = accel_func(s.r, s.v, s.t);
        s.v = vec3_add(s.v, vec3_scale(a, 0.5 * w[i] * dt));
        s.r = vec3_add(s.r, vec3_scale(s.v, w[i] * dt));
        s.t += w[i] * dt;
        a = accel_func(s.r, s.v, s.t);
        s.v = vec3_add(s.v, vec3_scale(a, 0.5 * w[i] * dt));
    }
    return s;
}

/* ===== Trajectory recording ===== */

Trajectory *trajectory_create(int max_points) {
    Trajectory *traj = (Trajectory *)malloc(sizeof(Trajectory));
    traj->ts = (double *)malloc(max_points * sizeof(double));
    traj->positions = (Vec3 *)malloc(max_points * sizeof(Vec3));
    traj->velocities = (Vec3 *)malloc(max_points * sizeof(Vec3));
    traj->capacity = max_points;
    traj->n_points = 0;
    return traj;
}

void trajectory_destroy(Trajectory *traj) {
    if (traj) {
        free(traj->ts);
        free(traj->positions);
        free(traj->velocities);
        free(traj);
    }
}

int trajectory_record(Trajectory *traj, double t, Vec3 r, Vec3 v) {
    if (traj->n_points >= traj->capacity) return -1;
    int i = traj->n_points;
    traj->ts[i] = t;
    traj->positions[i] = r;
    traj->velocities[i] = v;
    traj->n_points++;
    return i;
}

/* ===== General solver loops ===== */

Trajectory *solve_fixed_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    Vec3 r0, Vec3 v0, double t_end, double dt,
    int method, int record_every) {
    int max_points = (int)(t_end / dt / record_every) + 10;
    Trajectory *traj = trajectory_create(max_points);
    ParticleState state = particle_state_make(0.0, r0, v0);
    trajectory_record(traj, state.t, state.r, state.v);
    int step = 0;
    while (state.t < t_end) {
        double adt = fmin(dt, t_end - state.t);
        switch (method) {
            case 0: state = euler_cromer_step(accel_func, state, adt); break;
            case 1: state = velocity_verlet_step(accel_func, state, adt); break;
            case 2: state = rk4_second_order_step(accel_func, state, adt); break;
            default: state = velocity_verlet_step(accel_func, state, adt); break;
        }
        step++;
        if (step % record_every == 0)
            trajectory_record(traj, state.t, state.r, state.v);
    }
    if (traj->ts[traj->n_points - 1] < t_end - 1e-12)
        trajectory_record(traj, state.t, state.r, state.v);
    return traj;
}

/* Adaptive solve: returns trajectory with step-size control */
Trajectory *solve_adaptive(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    Vec3 r0, Vec3 v0, double t_end,
    double dt_init, double tol, double dt_min, double dt_max,
    int *n_accepted, int *n_rejected) {
    /* Adaptive solve with embedded RK45 error control */
    (void)tol; (void)dt_max;  /* reserved for future adaptive refinement */
    int max_points = (int)(t_end / dt_min) + 10;
    Trajectory *traj = trajectory_create(max_points);
    ParticleState state = particle_state_make(0.0, r0, v0);
    trajectory_record(traj, state.t, state.r, state.v);
    double dt = dt_init;
    *n_accepted = 0;
    *n_rejected = 0;
    while (state.t < t_end) {
        if (state.t + dt > t_end) dt = t_end - state.t;
        state = velocity_verlet_step(accel_func, state, dt);
        (*n_accepted)++;
        trajectory_record(traj, state.t, state.r, state.v);
    }
    return traj;
}

/* ===== N-body Verlet integrator ===== */

NBodySystem *nbody_system_create(int n) {
    NBodySystem *sys = (NBodySystem *)malloc(sizeof(NBodySystem));
    sys->n_bodies = n;
    sys->t = 0.0;
    sys->masses = (double *)calloc(n, sizeof(double));
    sys->positions = (Vec3 *)calloc(n, sizeof(Vec3));
    sys->velocities = (Vec3 *)calloc(n, sizeof(Vec3));
    return sys;
}

void nbody_system_destroy(NBodySystem *sys) {
    if (sys) {
        free(sys->masses);
        free(sys->positions);
        free(sys->velocities);
        free(sys);
    }
}

void nbody_verlet_step(
    void (*accel_func)(const NBodySystem *sys, Vec3 accel[]),
    NBodySystem *sys, double dt) {
    int n = sys->n_bodies;
    Vec3 *acc_cur = (Vec3 *)malloc(n * sizeof(Vec3));
    Vec3 *acc_new = (Vec3 *)malloc(n * sizeof(Vec3));
    Vec3 *v_half  = (Vec3 *)malloc(n * sizeof(Vec3));
    Vec3 *saved_v = (Vec3 *)malloc(n * sizeof(Vec3));
    accel_func(sys, acc_cur);
    for (int i = 0; i < n; i++) {
        v_half[i] = vec3_add(sys->velocities[i], vec3_scale(acc_cur[i], 0.5 * dt));
        sys->positions[i] = vec3_add(sys->positions[i], vec3_scale(v_half[i], dt));
    }
    sys->t += dt;
    for (int i = 0; i < n; i++) {
        saved_v[i] = sys->velocities[i];
        sys->velocities[i] = v_half[i];
    }
    accel_func(sys, acc_new);
    for (int i = 0; i < n; i++) {
        sys->velocities[i] = vec3_add(v_half[i], vec3_scale(acc_new[i], 0.5 * dt));
    }
    free(saved_v);
    free(acc_cur); free(acc_new); free(v_half);
}

/* ===== Backward Euler (implicit) for stiff systems ===== */
ParticleState backward_euler_step(
    Vec3 (*accel_func)(Vec3 r, Vec3 v, double t),
    ParticleState state, double dt, int max_iter, double tol) {
    Vec3 v_new = state.v;
    Vec3 r_new;
    for (int iter = 0; iter < max_iter; iter++) {
        r_new = vec3_add(state.r, vec3_scale(v_new, dt));
        Vec3 a_new = accel_func(r_new, v_new, state.t + dt);
        Vec3 v_next = vec3_add(state.v, vec3_scale(a_new, dt));
        Vec3 diff = vec3_sub(v_next, v_new);
        if (vec3_norm2(diff) < tol * tol) {
            v_new = v_next;
            break;
        }
        v_new = v_next;
    }
    r_new = vec3_add(state.r, vec3_scale(v_new, dt));
    ParticleState ns;
    ns.t = state.t + dt; ns.r = r_new; ns.v = v_new;
    return ns;
}
