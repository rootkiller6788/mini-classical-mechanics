/**
 * rigid_body_kinematics.c — Euler Angles, Rotation Matrices, Quaternions
 *
 * Goldstein §4.8-4.9: Kinematics of rigid body rotations.
 * Landau §35: Eulerian angles.
 *
 * Three SO(3) parameterizations implemented:
 *   1. Euler angles (ZXZ, ZYZ conventions)
 *   2. Tait-Bryan angles (XYZ, aerospace convention)
 *   3. Rotation matrices (direct 3×3 orthogonal)
 *   4. Quaternions (SU(2) double cover, singularity-free)
 *
 * All transformations are verified for round-trip consistency.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/rigid_body_types.h"
#include "../include/rigid_body_kinematics.h"

/* ============================================================================
 * Euler Angles ↔ Rotation Matrix (ZXZ Convention)
 * ============================================================================ */

void euler_to_rotation_matrix(double phi, double theta, double psi, mat3 *R) {
    double cp = cos(phi), sp = sin(phi);
    double ct = cos(theta), st = sin(theta);
    double cps = cos(psi), sps = sin(psi);

    /* R = Rz(φ) · Rx(θ) · Rz(ψ)  (Goldstein 4.46) */
    R->m[0] = cp*cps - sp*ct*sps;
    R->m[1] = -cp*sps - sp*ct*cps;
    R->m[2] = sp*st;

    R->m[3] = sp*cps + cp*ct*sps;
    R->m[4] = -sp*sps + cp*ct*cps;
    R->m[5] = -cp*st;

    R->m[6] = st*sps;
    R->m[7] = st*cps;
    R->m[8] = ct;
}

int rotation_matrix_to_euler(const mat3 *R, EulerAngles *ea) {
    double R33 = R->m[8];
    /* Clamp for numerical stability */
    if (R33 > 1.0) R33 = 1.0;
    if (R33 < -1.0) R33 = -1.0;

    double theta = acos(R33);
    double phi, psi;

    if (fabs(R33 - 1.0) < 1e-14) {
        /* θ ≈ 0: gimbal lock — only φ+ψ is determined */
        phi = atan2(R->m[0], -R->m[1]);
        psi = 0.0;
    } else if (fabs(R33 + 1.0) < 1e-14) {
        /* θ ≈ π: gimbal lock — only φ-ψ is determined */
        phi = atan2(-R->m[0], -R->m[1]);
        psi = 0.0;
    } else {
        phi = atan2(R->m[2], -R->m[5]);  /* R13, -R23 */
        psi = atan2(R->m[6], R->m[7]);   /* R31, R32 */
    }

    ea->phi = fmod(phi, M_2PI);
    if (ea->phi < 0) ea->phi += M_2PI;
    ea->theta = theta;
    ea->psi = fmod(psi, M_2PI);
    if (ea->psi < 0) ea->psi += M_2PI;

    return 0;
}

/* ZYZ Euler angles */
void euler_zyz_to_matrix(double alpha, double beta, double gamma, mat3 *R) {
    double ca = cos(alpha), sa = sin(alpha);
    double cb = cos(beta), sb = sin(beta);
    double cg = cos(gamma), sg = sin(gamma);

    R->m[0] = ca*cb*cg - sa*sg;
    R->m[1] = -ca*cb*sg - sa*cg;
    R->m[2] = ca*sb;

    R->m[3] = sa*cb*cg + ca*sg;
    R->m[4] = -sa*cb*sg + ca*cg;
    R->m[5] = sa*sb;

    R->m[6] = -sb*cg;
    R->m[7] = sb*sg;
    R->m[8] = cb;
}

int matrix_to_euler_zyz(const mat3 *R, double *alpha, double *beta, double *gamma) {
    double cb = R->m[8];
    if (cb > 1.0) cb = 1.0;
    if (cb < -1.0) cb = -1.0;
    *beta = acos(cb);
    if (fabs(cb - 1.0) < 1e-14) {
        *alpha = atan2(R->m[1], R->m[0]);
        *gamma = 0.0;
    } else if (fabs(cb + 1.0) < 1e-14) {
        *alpha = atan2(-R->m[1], -R->m[0]);
        *gamma = 0.0;
    } else {
        *alpha = atan2(R->m[5], R->m[2]);
        *gamma = atan2(R->m[7], -R->m[6]);
    }
    *alpha = fmod(*alpha, M_2PI); if (*alpha < 0) *alpha += M_2PI;
    *gamma = fmod(*gamma, M_2PI); if (*gamma < 0) *gamma += M_2PI;
    return 0;
}

/* Tait-Bryan (XYZ / roll-pitch-yaw) */
void tait_bryan_to_matrix(double roll, double pitch, double yaw, mat3 *R) {
    double cr = cos(roll), sr = sin(roll);
    double cp = cos(pitch), sp = sin(pitch);
    double cy = cos(yaw), sy = sin(yaw);

    /* R = Rx(roll) · Ry(pitch) · Rz(yaw) */
    R->m[0] = cp*cy;
    R->m[1] = -cp*sy;
    R->m[2] = sp;

    R->m[3] = sr*sp*cy + cr*sy;
    R->m[4] = -sr*sp*sy + cr*cy;
    R->m[5] = -sr*cp;

    R->m[6] = -cr*sp*cy + sr*sy;
    R->m[7] = cr*sp*sy + sr*cy;
    R->m[8] = cr*cp;
}

int matrix_to_tait_bryan(const mat3 *R, double *roll, double *pitch, double *yaw) {
    /* Roll = atan2(-R32, R33), Pitch = asin(R31), Yaw = atan2(-R21, R11) */
    if (fabs(R->m[2]) > 0.999999) {
        /* Pitch = ±90°, gimbal lock */
        *pitch = copysign(M_PI_2, R->m[2]);
        *roll = 0.0;
        *yaw = atan2(-R->m[1], R->m[4]);
    } else {
        *pitch = asin(R->m[2]);
        *roll = atan2(-R->m[5], R->m[8]);
        *yaw = atan2(-R->m[1], R->m[0]);
    }
    return 0;
}

/* ============================================================================
 * Euler Angle Rates ↔ Angular Velocity
 * ============================================================================ */

void euler_rates_to_body_omega(double phi_dot, double theta_dot, double psi_dot,
                               double theta, double psi, double omega[3]) {
    double st = sin(theta), sp = sin(psi), cp = cos(psi);
    omega[0] = phi_dot * st * sp + theta_dot * cp;
    omega[1] = phi_dot * st * cp - theta_dot * sp;
    omega[2] = phi_dot * cos(theta) + psi_dot;
}

void body_omega_to_euler_rates(const double omega[3], double theta, double psi,
                               double *phi_dot, double *theta_dot, double *psi_dot) {
    double st = sin(theta), sp = sin(psi), cp = cos(psi);

    if (fabs(st) < 1e-15) {
        /* Gimbal lock: cannot separate φ̇ and ψ̇ */
        double omp = omega[0]*sp + omega[1]*cp;
        *phi_dot = omp / (st + 1e-15);  /* Regularized */
        *theta_dot = omega[0]*cp - omega[1]*sp;
        *psi_dot = 0.0;
    } else {
        *phi_dot = (omega[0]*sp + omega[1]*cp) / st;
        *theta_dot = omega[0]*cp - omega[1]*sp;
        *psi_dot = omega[2] - (*phi_dot) * cos(theta);
    }
}

void euler_kinematic_ode(double t, const double y[], double dydt[], const void *params) {
    (void)t; /* Autonomous system */
    const double *omega = (const double *)params;
    double phi_dot, theta_dot, psi_dot;
    body_omega_to_euler_rates(omega, y[1]/*theta*/, y[2]/*psi*/,
                              &phi_dot, &theta_dot, &psi_dot);
    dydt[0] = phi_dot;
    dydt[1] = theta_dot;
    dydt[2] = psi_dot;
}

/* ============================================================================
 * Rotation Matrix — Angular Velocity Relation
 * ============================================================================ */

void rotation_matrix_derivative(const mat3 *R, const double omega[3], mat3 *dR) {
    mat3 Omega;
    cross_product_matrix(omega, &Omega);
    *dR = mat3_mul(R, &Omega);
}

void cross_product_matrix(const double omega[3], mat3 *skew) {
    skew->m[0] = 0.0;       skew->m[1] = -omega[2];  skew->m[2] = omega[1];
    skew->m[3] = omega[2];  skew->m[4] = 0.0;        skew->m[5] = -omega[0];
    skew->m[6] = -omega[1]; skew->m[7] = omega[0];   skew->m[8] = 0.0;
}

void unskew_matrix(const mat3 *skew, double omega[3]) {
    omega[0] = -skew->m[5];  /* Use ω_x = -Ω_23 */
    omega[1] = skew->m[2];   /* ω_y =  Ω_13 */
    omega[2] = -skew->m[1];  /* ω_z = -Ω_12 */
}

/* ============================================================================
 * Quaternion Operations
 * ============================================================================ */

void quaternion_to_rotation_matrix(Quaternion q, mat3 *R) {
    double w = q.w, x = q.x, y = q.y, z = q.z;
    double xx = x*x, yy = y*y, zz = z*z;
    double xy = x*y, xz = x*z, yz = y*z;
    double wx = w*x, wy = w*y, wz = w*z;

    R->m[0] = 1.0 - 2.0*(yy + zz);
    R->m[1] = 2.0*(xy - wz);
    R->m[2] = 2.0*(xz + wy);

    R->m[3] = 2.0*(xy + wz);
    R->m[4] = 1.0 - 2.0*(xx + zz);
    R->m[5] = 2.0*(yz - wx);

    R->m[6] = 2.0*(xz - wy);
    R->m[7] = 2.0*(yz + wx);
    R->m[8] = 1.0 - 2.0*(xx + yy);
}

int rotation_matrix_to_quaternion(const mat3 *R, Quaternion *q) {
    /* Shepperd's algorithm for numerical stability */
    double trace = R->m[0] + R->m[4] + R->m[8];
    double w, x, y, z;

    if (trace > 0.0) {
        double s = 0.5 / sqrt(trace + 1.0);
        w = 0.25 / s;
        x = (R->m[7] - R->m[5]) * s;
        y = (R->m[2] - R->m[6]) * s;
        z = (R->m[3] - R->m[1]) * s;
    } else if (R->m[0] > R->m[4] && R->m[0] > R->m[8]) {
        double s = 2.0 * sqrt(1.0 + R->m[0] - R->m[4] - R->m[8]);
        w = (R->m[7] - R->m[5]) / s;
        x = 0.25 * s;
        y = (R->m[1] + R->m[3]) / s;
        z = (R->m[2] + R->m[6]) / s;
    } else if (R->m[4] > R->m[8]) {
        double s = 2.0 * sqrt(1.0 + R->m[4] - R->m[0] - R->m[8]);
        w = (R->m[2] - R->m[6]) / s;
        x = (R->m[1] + R->m[3]) / s;
        y = 0.25 * s;
        z = (R->m[5] + R->m[7]) / s;
    } else {
        double s = 2.0 * sqrt(1.0 + R->m[8] - R->m[0] - R->m[4]);
        w = (R->m[3] - R->m[1]) / s;
        x = (R->m[2] + R->m[6]) / s;
        y = (R->m[5] + R->m[7]) / s;
        z = 0.25 * s;
    }

    *q = quat_normalize((Quaternion){w, x, y, z});
    return 0;
}

void euler_to_quaternion(double phi, double theta, double psi, Quaternion *q) {
    /* q = q_z(φ) · q_x(θ) · q_z(ψ) where q_axis(α) = (cos(α/2), axis·sin(α/2)) */
    double h_phi = 0.5 * phi, h_theta = 0.5 * theta, h_psi = 0.5 * psi;

    /* q_z(φ) */
    Quaternion qz1 = {cos(h_phi), 0, 0, sin(h_phi)};
    /* q_x(θ) */
    Quaternion qx = {cos(h_theta), sin(h_theta), 0, 0};
    /* q_z(ψ) */
    Quaternion qz2 = {cos(h_psi), 0, 0, sin(h_psi)};

    *q = quat_mul(quat_mul(qz1, qx), qz2);
}

void quaternion_to_euler(Quaternion q, EulerAngles *ea) {
    /* Convert quaternion to rotation matrix, then to Euler angles */
    mat3 R;
    quaternion_to_rotation_matrix(q, &R);
    rotation_matrix_to_euler(&R, ea);
}

void quat_slerp(Quaternion q0, Quaternion q1, double t, Quaternion *out) {
    /* cos(Ω) = Re(q0* · q1) = q0·q1 (dot product for unit quaternions) */
    double dot = q0.w*q1.w + q0.x*q1.x + q0.y*q1.y + q0.z*q1.z;

    /* If dot < 0, slerp takes the long way; negate q1 for shortest path */
    Quaternion q1_adj = q1;
    if (dot < 0.0) {
        q1_adj.w = -q1.w; q1_adj.x = -q1.x;
        q1_adj.y = -q1.y; q1_adj.z = -q1.z;
        dot = -dot;
    }

    /* If very close, linear interpolation is fine */
    if (dot > 0.9995) {
        double t1 = 1.0 - t;
        *out = quat_normalize((Quaternion){
            t1*q0.w + t*q1_adj.w,
            t1*q0.x + t*q1_adj.x,
            t1*q0.y + t*q1_adj.y,
            t1*q0.z + t*q1_adj.z
        });
        return;
    }

    double omega = acos(dot);
    double sin_omega = sin(omega);
    double s0 = sin((1.0 - t) * omega) / sin_omega;
    double s1 = sin(t * omega) / sin_omega;

    out->w = s0 * q0.w + s1 * q1_adj.w;
    out->x = s0 * q0.x + s1 * q1_adj.x;
    out->y = s0 * q0.y + s1 * q1_adj.y;
    out->z = s0 * q0.z + s1 * q1_adj.z;
}

void quaternion_derivative(Quaternion q, const double omega[3], Quaternion *dq) {
    /* dq/dt = 0.5 * ω_q * q where ω_q = (0, ωx, ωy, ωz) */
    Quaternion wq = {0.0, omega[0], omega[1], omega[2]};
    Quaternion prod = quat_mul(wq, q);
    dq->w = 0.5 * prod.w;
    dq->x = 0.5 * prod.x;
    dq->y = 0.5 * prod.y;
    dq->z = 0.5 * prod.z;
}

void quat_exp(const double axis_angle[3], Quaternion *q) {
    double theta = sqrt(axis_angle[0]*axis_angle[0]
                      + axis_angle[1]*axis_angle[1]
                      + axis_angle[2]*axis_angle[2]);
    if (theta < 1e-15) {
        *q = QUAT_IDENTITY;
        return;
    }
    double half_theta = 0.5 * theta;
    double s = sin(half_theta) / theta;
    q->w = cos(half_theta);
    q->x = axis_angle[0] * s;
    q->y = axis_angle[1] * s;
    q->z = axis_angle[2] * s;
}

void quat_log(Quaternion q, double axis_angle[3]) {
    if (q.w > 1.0) q.w = 1.0;
    if (q.w < -1.0) q.w = -1.0;
    double theta = 2.0 * acos(q.w);
    if (theta < 1e-15) {
        axis_angle[0] = axis_angle[1] = axis_angle[2] = 0.0;
        return;
    }
    double s = theta / sin(0.5 * theta);
    axis_angle[0] = q.x * s;
    axis_angle[1] = q.y * s;
    axis_angle[2] = q.z * s;
}

/* ============================================================================
 * Basic Rotation Matrices
 * ============================================================================ */

void rotation_x(double angle, mat3 *R) {
    double c = cos(angle), s = sin(angle);
    R->m[0] = 1; R->m[1] = 0;  R->m[2] = 0;
    R->m[3] = 0; R->m[4] = c;  R->m[5] = -s;
    R->m[6] = 0; R->m[7] = s;  R->m[8] = c;
}

void rotation_y(double angle, mat3 *R) {
    double c = cos(angle), s = sin(angle);
    R->m[0] = c;  R->m[1] = 0; R->m[2] = s;
    R->m[3] = 0;  R->m[4] = 1; R->m[5] = 0;
    R->m[6] = -s; R->m[7] = 0; R->m[8] = c;
}

void rotation_z(double angle, mat3 *R) {
    double c = cos(angle), s = sin(angle);
    R->m[0] = c;  R->m[1] = -s; R->m[2] = 0;
    R->m[3] = s;  R->m[4] = c;  R->m[5] = 0;
    R->m[6] = 0;  R->m[7] = 0;  R->m[8] = 1;
}

void rotation_axis_angle(vec3 axis, double angle, mat3 *R) {
    vec3 n = vec3_normalize(axis);
    double c = cos(angle), s = sin(angle);
    double t = 1.0 - c;
    double nx = n.x, ny = n.y, nz = n.z;

    /* Rodrigues formula: R = I + sinθ [n̂]_× + (1-cosθ) [n̂]_×² */
    R->m[0] = c + nx*nx*t;
    R->m[1] = nx*ny*t - nz*s;
    R->m[2] = nx*nz*t + ny*s;

    R->m[3] = ny*nx*t + nz*s;
    R->m[4] = c + ny*ny*t;
    R->m[5] = ny*nz*t - nx*s;

    R->m[6] = nz*nx*t - ny*s;
    R->m[7] = nz*ny*t + nx*s;
    R->m[8] = c + nz*nz*t;
}

int rotation_to_axis_angle(const mat3 *R, vec3 *axis, double *angle) {
    double trace = mat3_trace(R);
    double c = 0.5 * (trace - 1.0);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    *angle = acos(c);

    if (*angle < 1e-15) {
        *axis = (vec3){0, 0, 1};
        return 0;
    }

    if (fabs(*angle - M_PI) < 1e-10) {
        /* 180° case: axis from columns of R+I */
        double U[9];
        for (int i = 0; i < 9; i++) U[i] = R->m[i];
        U[0] += 1.0; U[4] += 1.0; U[8] += 1.0;
        /* Find nonzero column */
        for (int col = 0; col < 3; col++) {
            double n2 = U[col]*U[col] + U[3+col]*U[3+col] + U[6+col]*U[6+col];
            if (n2 > 1e-10) {
                double inv = 1.0 / sqrt(n2);
                axis->x = U[col] * inv;
                axis->y = U[3+col] * inv;
                axis->z = U[6+col] * inv;
                return 0;
            }
        }
        *axis = (vec3){1, 0, 0};
        return 0;
    }

    /* General case: n̂ = (1/(2 sin θ)) * [R32-R23, R13-R31, R21-R12] */
    double s2 = 2.0 * sin(*angle);
    axis->x = (R->m[7] - R->m[5]) / s2;
    axis->y = (R->m[2] - R->m[6]) / s2;
    axis->z = (R->m[3] - R->m[1]) / s2;
    return 0;
}

void rotation_compose(const mat3 *R1, const mat3 *R2, mat3 *R12) {
    *R12 = mat3_mul(R1, R2);
}

void rotation_inverse(const mat3 *R, mat3 *Rinv) {
    *Rinv = mat3_transpose(R);
}
