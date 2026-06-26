/**
 * rigid_body_types.c — Implementation of core type operations
 *
 * Reference: Goldstein §4-5, Landau §32-35
 *
 * This file implements constructors, display functions, and basic operations
 * for the foundational data types of rigid body mechanics.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/rigid_body_types.h"

/* ============================================================================
 * 3×3 Matrix Operations
 * ============================================================================ */

mat3 mat3_mul(const mat3 *A, const mat3 *B) {
    mat3 C;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) {
                sum += mat3_get(A, i, k) * mat3_get(B, k, j);
            }
            mat3_set(&C, i, j, sum);
        }
    }
    return C;
}

double mat3_frobenius_norm(const mat3 *M) {
    double sum = 0.0;
    for (int i = 0; i < 9; i++) {
        sum += M->m[i] * M->m[i];
    }
    return sqrt(sum);
}

int mat3_is_orthogonal(const mat3 *M, double tol) {
    mat3 MT = mat3_transpose(M);
    mat3 prod = mat3_mul(M, &MT);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (fabs(mat3_get(&prod, i, j) - expected) > tol) return 0;
        }
    }
    return 1;
}

double mat3_det(const mat3 *M) {
    /* Sarrus' rule for 3×3 determinant */
    return M->m[0] * (M->m[4] * M->m[8] - M->m[5] * M->m[7])
         - M->m[1] * (M->m[3] * M->m[8] - M->m[5] * M->m[6])
         + M->m[2] * (M->m[3] * M->m[7] - M->m[4] * M->m[6]);
}

/* ============================================================================
 * InertiaTensor ↔ mat3 Conversion
 * ============================================================================ */

void inertia_to_mat3(const InertiaTensor *I, mat3 *out) {
    out->m[0] = I->Ixx; out->m[1] = I->Ixy; out->m[2] = I->Ixz;
    out->m[3] = I->Ixy; out->m[4] = I->Iyy; out->m[5] = I->Iyz;
    out->m[6] = I->Ixz; out->m[7] = I->Iyz; out->m[8] = I->Izz;
}

InertiaTensor mat3_to_inertia(const mat3 *M) {
    /* Symmetrize: average off-diagonal elements */
    return (InertiaTensor){
        M->m[0],
        M->m[4],
        M->m[8],
        0.5 * (M->m[1] + M->m[3]),
        0.5 * (M->m[2] + M->m[6]),
        0.5 * (M->m[5] + M->m[7])
    };
}

/* ============================================================================
 * Euler Angle Normalization
 * ============================================================================ */

void euler_normalize(EulerAngles *ea) {
    /* phi ∈ [0, 2π) */
    ea->phi = fmod(ea->phi, M_2PI);
    if (ea->phi < 0) ea->phi += M_2PI;
    /* theta ∈ [0, π] — clamp */
    ea->theta = fmod(ea->theta, M_2PI);
    if (ea->theta < 0) ea->theta += M_2PI;
    if (ea->theta > M_PI) {
        /* Map (π, 2π) back: θ→2π-θ, φ→φ+π, ψ→ψ+π (equivalent orientation) */
        ea->theta = M_2PI - ea->theta;
        ea->phi = fmod(ea->phi + M_PI, M_2PI);
        ea->psi = fmod(ea->psi + M_PI, M_2PI);
    }
    /* psi ∈ [0, 2π) */
    ea->psi = fmod(ea->psi, M_2PI);
    if (ea->psi < 0) ea->psi += M_2PI;
}

/* ============================================================================
 * RigidBodyState Constructor
 * ============================================================================ */

RigidBodyState rb_state_new(vec3 omega, EulerAngles ea) {
    return (RigidBodyState){0.0, omega, ea, VEC3_ZERO, 0.0};
}

/* ============================================================================
 * Quaternion Operations
 * ============================================================================ */

Quaternion quat_mul(Quaternion q1, Quaternion q2) {
    /* Hamilton product: (w1 + v1)(w2 + v2) = w1w2 - v1·v2 + w1v2 + w2v1 + v1×v2 */
    return (Quaternion){
        q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z,
        q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
        q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
        q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w
    };
}

Quaternion quat_normalize(Quaternion q) {
    double n = quat_norm(q);
    if (n < 1e-15) return QUAT_IDENTITY;
    return (Quaternion){q.w / n, q.x / n, q.y / n, q.z / n};
}

vec3 quat_rotate_vec(Quaternion q, vec3 v) {
    /* v' = q * v * q* where v is treated as pure quaternion (0, vx, vy, vz) */
    Quaternion qv = {0.0, v.x, v.y, v.z};
    Quaternion qc = quat_conj(q);
    Quaternion result = quat_mul(quat_mul(q, qv), qc);
    return (vec3){result.x, result.y, result.z};
}

/* ============================================================================
 * Display / Debug Functions
 * ============================================================================ */

void inertia_print(const InertiaTensor *I) {
    printf("InertiaTensor:\n");
    printf("  Ixx=%.6f  Iyy=%.6f  Izz=%.6f\n", I->Ixx, I->Iyy, I->Izz);
    printf("  Ixy=%.6f  Ixz=%.6f  Iyz=%.6f\n", I->Ixy, I->Ixz, I->Iyz);
}

void euler_print(const EulerAngles *ea) {
    printf("EulerAngles(phi=%.4f, theta=%.4f, psi=%.4f)\n",
           ea->phi, ea->theta, ea->psi);
}

void mat3_print(const mat3 *M) {
    printf("[%.6f  %.6f  %.6f]\n", M->m[0], M->m[1], M->m[2]);
    printf("[%.6f  %.6f  %.6f]\n", M->m[3], M->m[4], M->m[5]);
    printf("[%.6f  %.6f  %.6f]\n", M->m[6], M->m[7], M->m[8]);
}

void vec3_print(vec3 v) {
    printf("(%.6f, %.6f, %.6f)\n", v.x, v.y, v.z);
}

void quat_print(Quaternion q) {
    printf("Quaternion(w=%.6f, x=%.6f, y=%.6f, z=%.6f)\n", q.w, q.x, q.y, q.z);
}

void principal_axes_print(const PrincipalAxes *pa) {
    printf("PrincipalAxes:\n");
    printf("  moments: [%.6f, %.6f, %.6f]\n",
           pa->moments[0], pa->moments[1], pa->moments[2]);
    printf("  axes:\n");
    mat3_print(&pa->axes);
}

void rigid_body_state_print(const RigidBodyState *s) {
    printf("RigidBodyState @ t=%.4f:\n", s->t);
    printf("  omega: "); vec3_print(s->omega);
    printf("  euler: "); euler_print(&s->euler);
    printf("  L:     "); vec3_print(s->L);
    printf("  T:     %.6f\n", s->T);
}
