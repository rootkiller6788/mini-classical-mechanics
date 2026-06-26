/**
 * rigid_body_energy.c — Rotational Energy, Angular Momentum, Stability
 *
 * Goldstein §5.1-5.2, §5.6: Rotational kinetic energy, angular momentum vector,
 *                            Tennis racket theorem (intermediate axis instability).
 * Landau §37: Stability of rigid body rotation.
 *
 * The two fundamental integrals of torque-free motion are:
 *   1. Rotational kinetic energy: T = ½ Σ I_i ω_i²
 *   2. Angular momentum magnitude: L² = Σ (I_i ω_i)²
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/rigid_body_types.h"
#include "../include/rigid_body_energy.h"
#include "../include/rigid_body_euler.h"
#include "../include/rigid_body_inertia.h"

/* ============================================================================
 * Rotational Kinetic Energy
 * ============================================================================ */

double rotational_ke(const InertiaTensor *I, const double omega[3]) {
    /* T = ½ ωᵀ · I · ω  (Goldstein 5.9)
     * Full quadratic form: Ixx ωx² + Iyy ωy² + Izz ωz² + 2Ixy ωxωy + 2Ixz ωxωz + 2Iyz ωyωz */
    double wx = omega[0], wy = omega[1], wz = omega[2];
    return 0.5 * (I->Ixx*wx*wx + I->Iyy*wy*wy + I->Izz*wz*wz
                + 2.0*I->Ixy*wx*wy + 2.0*I->Ixz*wx*wz + 2.0*I->Iyz*wy*wz);
}

double rotational_ke_principal(const double I_principal[3], const double omega[3]) {
    return 0.5 * (I_principal[0]*omega[0]*omega[0]
                + I_principal[1]*omega[1]*omega[1]
                + I_principal[2]*omega[2]*omega[2]);
}

double total_kinetic_energy(double M, vec3 v_cm, const InertiaTensor *I_cm,
                            const double omega[3], double *T_trans, double *T_rot) {
    /* König's theorem: T = T_cm + T_rot_about_cm */
    if (T_trans) *T_trans = 0.5 * M * vec3_norm2(v_cm);
    if (T_rot) *T_rot = rotational_ke(I_cm, omega);
    double T_t = 0.5 * M * vec3_norm2(v_cm);
    double T_r = rotational_ke(I_cm, omega);
    return T_t + T_r;
}

/* ============================================================================
 * Angular Momentum
 * ============================================================================ */

void angular_momentum_from_inertia(const InertiaTensor *I, const double omega[3], double L[3]) {
    double wx = omega[0], wy = omega[1], wz = omega[2];
    L[0] = I->Ixx*wx + I->Ixy*wy + I->Ixz*wz;
    L[1] = I->Ixy*wx + I->Iyy*wy + I->Iyz*wz;
    L[2] = I->Ixz*wx + I->Iyz*wy + I->Izz*wz;
}

void angular_momentum_principal(const double I_principal[3], const double omega[3], double L[3]) {
    L[0] = I_principal[0] * omega[0];
    L[1] = I_principal[1] * omega[1];
    L[2] = I_principal[2] * omega[2];
}

void total_angular_momentum(vec3 r_cm, double M, vec3 v_cm,
                           const InertiaTensor *I_cm, const double omega[3],
                           double L_total[3], double L_orbital[3], double L_spin[3]) {
    /* L_orbital = r_cm × (M v_cm) */
    vec3 Mv = vec3_scale(M, v_cm);
    vec3 Lorb = vec3_cross(r_cm, Mv);
    L_orbital[0] = Lorb.x; L_orbital[1] = Lorb.y; L_orbital[2] = Lorb.z;

    /* L_spin = I_cm · ω */
    angular_momentum_from_inertia(I_cm, omega, L_spin);

    L_total[0] = L_orbital[0] + L_spin[0];
    L_total[1] = L_orbital[1] + L_spin[1];
    L_total[2] = L_orbital[2] + L_spin[2];
}

void angular_velocity_from_momentum(const double I_principal[3], const double L[3], double omega[3]) {
    /* For principal frame: ω_i = L_i / I_i */
    omega[0] = (I_principal[0] > 1e-15) ? L[0] / I_principal[0] : 0.0;
    omega[1] = (I_principal[1] > 1e-15) ? L[1] / I_principal[1] : 0.0;
    omega[2] = (I_principal[2] > 1e-15) ? L[2] / I_principal[2] : 0.0;
}

/* ============================================================================
 * Intermediate Axis Theorem — Stability Analysis
 *
 * Landau §37: Linearized analysis of Euler equations.
 * Let ω = Ω e_k + δω where Ω is the spin magnitude about principal axis k.
 *   For k=1 (max inertial): δω oscillates → stable
 *   For k=2 (mid inertial): δω grows exponentially → unstable (Dzhanibekov effect)
 *   For k=3 (min inertial): δω oscillates → stable
 *
 * Growth rate: λ² = Ω² (I_k - I_j)(I_i - I_k) / (I_i I_j)
 * where i,j are the two non-k axes.
 * ============================================================================ */

StabilityResult axis_stability_analysis(const double I_principal[3], int axis,
                                        double omega_mag, double *growth_rate) {
    if (axis < 0 || axis > 2) return STABILITY_DEGENERATE;

    int i = (axis == 0) ? 1 : 0;        /* one of the other two axes */
    int j = (axis == 0) ? 2 : (axis == 1) ? 2 : 1;

    double Ii = I_principal[i], Ij = I_principal[j], Ik = I_principal[axis];

    /* Check for degeneracy (symmetric case) */
    if (fabs(Ik - Ii) < 1e-14 || fabs(Ik - Ij) < 1e-14) {
        *growth_rate = 0.0;
        return STABILITY_DEGENERATE;
    }

    /* λ² = Ω² (I_k - I_j)(I_i - I_k) / (I_i I_j) */
    double lambda_sq = omega_mag * omega_mag * (Ik - Ij) * (Ii - Ik) / (Ii * Ij);

    if (lambda_sq > 1e-15) {
        *growth_rate = sqrt(lambda_sq);
        return STABILITY_UNSTABLE;
    } else {
        *growth_rate = 0.0;
        return STABILITY_STABLE;
    }
}

void full_stability_analysis(const double I_principal[3], StabilityResult results[3]) {
    for (int axis = 0; axis < 3; axis++) {
        double gr;
        results[axis] = axis_stability_analysis(I_principal, axis, 1.0, &gr);
    }
}

double flipping_period_estimate(const double I_principal[3], double omega_mag,
                                double delta_omega0) {
    double growth_rate;
    StabilityResult stab = axis_stability_analysis(I_principal, 1, omega_mag, &growth_rate);
    if (stab != STABILITY_UNSTABLE || growth_rate < 1e-15) {
        return INFINITY;
    }
    return log(1.0 / delta_omega0) / growth_rate;
}

void dzhanibekov_simulation(const double I_principal[3], const double omega0[3],
                            double t_end, double dt,
                            int *n_steps_out, double *traj_out,
                            double *flip_times, int *n_flips_out) {
    int n_steps;
    int max_steps = (int)ceil(t_end / dt) + 1;
    simulate_free_rigid_body(I_principal, omega0, t_end, dt, INTEGRATOR_RK4,
                             &n_steps, traj_out);

    /* Detect flips: sign changes of the dominant (intermediate axis) component */
    int n_flips = 0;
    double prev_sign = (traj_out[1] >= 0) ? 1.0 : -1.0; /* start with ω_y */

    for (int i = 1; i < n_steps && n_flips < 20; i++) {
        double cur_sign = (traj_out[i*3 + 1] >= 0) ? 1.0 : -1.0;
        if (cur_sign * prev_sign < 0) {
            flip_times[n_flips] = i * dt;
            n_flips++;
        }
        prev_sign = cur_sign;
    }

    *n_steps_out = n_steps;
    *n_flips_out = n_flips;
}

/* ============================================================================
 * Energy-Momentum Geometry
 * ============================================================================ */

void energy_extrema_for_L2(const double I_principal[3], double L2,
                           double *T_min, double *T_max) {
    /* T = L²/(2I_eff). For fixed L², T is minimized when rotation is about the
     * axis of largest inertia (I_max → smallest ω → smallest T for given L²).
     * T_min = L² / (2 I_max), T_max = L² / (2 I_min) */
    double I_max = I_principal[0]; /* sorted descending: I1 ≥ I2 ≥ I3 */
    double I_min = I_principal[2];
    *T_min = L2 / (2.0 * I_max);
    *T_max = L2 / (2.0 * I_min);
}

int separatrix_energy(const double I_principal[3], double L2, double *T_sep) {
    /* Separatrix corresponds to rotation about intermediate axis: T = L²/(2 I₂) */
    double T_min, T_max;
    energy_extrema_for_L2(I_principal, L2, &T_min, &T_max);
    *T_sep = L2 / (2.0 * I_principal[1]);
    if (*T_sep < T_min || *T_sep > T_max) return -1;
    return 0;
}

void energy_on_momentum_sphere(const double I_principal[3], double L2,
                               int n_theta, int n_phi,
                               double *T_map, double *omega1, double *omega2, double *omega3) {
    for (int it = 0; it < n_theta; it++) {
        double theta = M_PI * it / (n_theta - 1);
        double st = sin(theta), ct = cos(theta);
        for (int ip = 0; ip < n_phi; ip++) {
            double phi = M_2PI * ip / n_phi;
            double sp = sin(phi), cp = cos(phi);

            /* Unit direction on sphere */
            double n1 = st * cp, n2 = st * sp, n3 = ct;

            /* Scale to angular momentum sphere: I₁²ω₁² + I₂²ω₂² + I₃²ω₃² = L² */
            double denom_sq = (I_principal[0]*n1)*(I_principal[0]*n1)
                            + (I_principal[1]*n2)*(I_principal[1]*n2)
                            + (I_principal[2]*n3)*(I_principal[2]*n3);
            double scale = sqrt(L2) / sqrt(fmax(denom_sq, 1e-30));

            double w1 = n1 * scale, w2 = n2 * scale, w3 = n3 * scale;
            int idx = it * n_phi + ip;
            if (omega1) omega1[idx] = w1;
            if (omega2) omega2[idx] = w2;
            if (omega3) omega3[idx] = w3;
            if (T_map) T_map[idx] = rotational_ke_principal(I_principal, (double[]){w1,w2,w3});
        }
    }
}

/* ============================================================================
 * Torque, Power, and Impulse
 * ============================================================================ */

double torque_power(const double torque[3], const double omega[3]) {
    return torque[0]*omega[0] + torque[1]*omega[1] + torque[2]*omega[2];
}

void torque_impulse(const double torque[3], double dt, double dL[3]) {
    dL[0] = torque[0] * dt;
    dL[1] = torque[1] * dt;
    dL[2] = torque[2] * dt;
}

double torque_work(int n_steps, const double *times, const double *torques, const double *omegas) {
    /* Trapezoidal rule: ∫ N·ω dt */
    double work = 0.0;
    for (int i = 0; i < n_steps - 1; i++) {
        double dt = times[i+1] - times[i];
        double P_i = torques[i*3]*omegas[i*3] + torques[i*3+1]*omegas[i*3+1] + torques[i*3+2]*omegas[i*3+2];
        double P_ip1 = torques[(i+1)*3]*omegas[(i+1)*3] + torques[(i+1)*3+1]*omegas[(i+1)*3+1] + torques[(i+1)*3+2]*omegas[(i+1)*3+2];
        work += 0.5 * dt * (P_i + P_ip1);
    }
    return work;
}

void torque_angular_impulse_total(int n_steps, const double *torques, double dt, double delta_L[3]) {
    delta_L[0] = delta_L[1] = delta_L[2] = 0.0;
    for (int i = 0; i < n_steps; i++) {
        delta_L[0] += torques[i*3] * dt;
        delta_L[1] += torques[i*3+1] * dt;
        delta_L[2] += torques[i*3+2] * dt;
    }
}
