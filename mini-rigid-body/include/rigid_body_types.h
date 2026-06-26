/**
 * rigid_body_types.h — Core Data Types for Rigid Body Mechanics
 *
 * Reference: Goldstein, Poole & Safko "Classical Mechanics" (3rd Ed.) Ch.4-5
 *            Landau & Lifshitz "Mechanics" Vol.1 Ch.6
 *            MIT 8.012 — Physics I: Classical Mechanics
 *
 * All types are designed for zero external dependency beyond standard C math.
 */

#ifndef RIGID_BODY_TYPES_H
#define RIGID_BODY_TYPES_H

#include <stddef.h>
#include <math.h>

/* Mathematical constants (if not defined by math.h) */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#ifndef M_2PI
#define M_2PI 6.28318530717958647692
#endif

/* ============================================================================
 * 3D Vector Operations (inline helpers)
 * ============================================================================ */

/** 3-vector with named components for clarity in physics contexts */
typedef struct {
    double x, y, z;
} vec3;

/** Zero vector constant */
static const vec3 VEC3_ZERO = {0.0, 0.0, 0.0};

/** Vector addition: c = a + b */
static inline vec3 vec3_add(vec3 a, vec3 b) {
    return (vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

/** Vector subtraction: c = a - b */
static inline vec3 vec3_sub(vec3 a, vec3 b) {
    return (vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

/** Scalar multiplication: c = s * a */
static inline vec3 vec3_scale(double s, vec3 a) {
    return (vec3){s * a.x, s * a.y, s * a.z};
}

/** Dot product (Euclidean inner product): a · b = a_x b_x + a_y b_y + a_z b_z */
static inline double vec3_dot(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/** Cross product: c = a × b (right-hand rule) */
static inline vec3 vec3_cross(vec3 a, vec3 b) {
    return (vec3){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

/** Euclidean norm (L2): |v| = √(v·v) */
static inline double vec3_norm(vec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/** Normalize to unit vector (safe against zero) */
static inline vec3 vec3_normalize(vec3 v) {
    double n = vec3_norm(v);
    if (n < 1e-15) return VEC3_ZERO;
    return vec3_scale(1.0 / n, v);
}

/** Vector squared norm: |v|² */
static inline double vec3_norm2(vec3 v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}


/* ============================================================================
 * 3x3 Matrix Type
 * ============================================================================ */

/** 3x3 matrix, row-major storage: m[0..2]=row0, m[3..5]=row1, m[6..8]=row2 */
typedef struct {
    double m[9];
} mat3;

/** Identity matrix */
static const mat3 MAT3_IDENTITY = {{1,0,0, 0,1,0, 0,0,1}};

/** Zero matrix */
static const mat3 MAT3_ZERO = {{0,0,0, 0,0,0, 0,0,0}};

/** Access element at (row, col), both 0-indexed */
static inline double mat3_get(const mat3 *M, int row, int col) {
    return M->m[row * 3 + col];
}

/** Set element at (row, col) */
static inline void mat3_set(mat3 *M, int row, int col, double val) {
    M->m[row * 3 + col] = val;
}

/** Matrix-vector multiply: y = M * x */
static inline vec3 mat3_mul_vec(const mat3 *M, vec3 x) {
    return (vec3){
        M->m[0]*x.x + M->m[1]*x.y + M->m[2]*x.z,
        M->m[3]*x.x + M->m[4]*x.y + M->m[5]*x.z,
        M->m[6]*x.x + M->m[7]*x.y + M->m[8]*x.z
    };
}

/** Matrix transpose: B = A^T */
static inline mat3 mat3_transpose(const mat3 *A) {
    return (mat3){{
        A->m[0], A->m[3], A->m[6],
        A->m[1], A->m[4], A->m[7],
        A->m[2], A->m[5], A->m[8]
    }};
}

/** Matrix multiply: C = A * B */
mat3 mat3_mul(const mat3 *A, const mat3 *B);

/** Matrix trace: Tr(M) = Σ M_ii */
static inline double mat3_trace(const mat3 *M) {
    return M->m[0] + M->m[4] + M->m[8];
}

/** Frobenius norm: ||M||_F = sqrt(Σ M_ij²) */
double mat3_frobenius_norm(const mat3 *M);

/** Check if matrix is orthogonal: M * M^T ≈ I */
int mat3_is_orthogonal(const mat3 *M, double tol);

/** Determinant of 3x3 matrix */
double mat3_det(const mat3 *M);


/* ============================================================================
 * InertiaTensor — 3×3 symmetric inertia tensor
 *
 * Goldstein (5.3): The moment of inertia tensor I is a real symmetric 3×3 matrix.
 * Diagonal: I_xx, I_yy, I_zz (moments of inertia)
 * Off-diagonal: I_xy, I_xz, I_yz (products of inertia, with negative sign)
 *
 * Only 6 independent components due to symmetry.
 * ============================================================================ */
typedef struct {
    double Ixx, Iyy, Izz;   /**< Diagonal: moments of inertia */
    double Ixy, Ixz, Iyz;   /**< Off-diagonal: products of inertia */
} InertiaTensor;

/** Zero inertia tensor */
static const InertiaTensor INERTIA_ZERO = {0,0,0, 0,0,0};

/** Convert InertiaTensor to 3×3 matrix form */
void inertia_to_mat3(const InertiaTensor *I, mat3 *out);

/** Create InertiaTensor from 3×3 matrix (symmetrized) */
InertiaTensor mat3_to_inertia(const mat3 *M);


/* ============================================================================
 * EulerAngles — ZXZ convention (Goldstein §4.8)
 *
 * φ (phi):   precession angle  ∈ [0, 2π)
 * θ (theta): nutation angle    ∈ [0, π]
 * ψ (psi):   spin angle        ∈ [0, 2π)
 *
 * Active rotation: R(φ,θ,ψ) = R_z(φ) · R_x(θ) · R_z(ψ)
 * ============================================================================ */
typedef struct {
    double phi;      /**< Precession angle (radians) */
    double theta;    /**< Nutation angle (radians) */
    double psi;      /**< Spin angle (radians) */
} EulerAngles;

/** Default: all zero angles */
static const EulerAngles EULER_ZERO = {0.0, 0.0, 0.0};

/** Normalize angles to canonical ranges */
void euler_normalize(EulerAngles *ea);


/* ============================================================================
 * PrincipalAxes — Diagonalization of inertia tensor
 *
 * Goldstein (5.4)-(5.6): The principal axes are the eigenvectors of I.
 * moments[0] ≥ moments[1] ≥ moments[2] (sorted descending).
 * axes columns are the corresponding eigenvectors (orthonormal).
 * ============================================================================ */
typedef struct {
    double moments[3];   /**< Principal moments [I1, I2, I3] */
    mat3   axes;         /**< 3×3 orthogonal matrix, columns = principal directions */
} PrincipalAxes;


/* ============================================================================
 * RigidBodyState — Complete kinematic/dynamic state
 *
 * Goldstein Ch.5: Full description of a rigid body's rotational motion.
 * - omega contains body-frame components [ω_x, ω_y, ω_z]
 * - euler angles describe orientation relative to space frame
 * - L is angular momentum in space frame
 * - T is rotational kinetic energy
 * ============================================================================ */
typedef struct {
    double      t;          /**< Time */
    vec3        omega;      /**< Angular velocity (body frame) */
    EulerAngles euler;      /**< Orientation (space frame) */
    vec3        L;          /**< Angular momentum (space frame) */
    double      T;          /**< Rotational kinetic energy */
} RigidBodyState;

/** Initialize state from omega and euler angles (t=0, L and T computed later) */
RigidBodyState rb_state_new(vec3 omega, EulerAngles ea);


/* ============================================================================
 * Quaternion — Hamilton convention (q_w + q_x i + q_y j + q_z k)
 *
 * Unit quaternions represent rotations in SU(2), double-covering SO(3).
 * ============================================================================ */
typedef struct {
    double w, x, y, z;   /**< w: scalar part; x,y,z: vector part */
} Quaternion;

/** Identity quaternion (no rotation) */
static const Quaternion QUAT_IDENTITY = {1.0, 0.0, 0.0, 0.0};

/** Quaternion norm squared: |q|² = w² + x² + y² + z² */
static inline double quat_norm2(Quaternion q) {
    return q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
}

/** Quaternion norm: |q| */
static inline double quat_norm(Quaternion q) {
    return sqrt(quat_norm2(q));
}

/** Check if quaternion is unit (to within tolerance) */
static inline int quat_is_unit(Quaternion q, double tol) {
    return fabs(quat_norm2(q) - 1.0) < tol;
}

/** Quaternion conjugate: q* = (w, -x, -y, -z) */
static inline Quaternion quat_conj(Quaternion q) {
    return (Quaternion){q.w, -q.x, -q.y, -q.z};
}

/** Quaternion multiply: q1 * q2 (Hamilton product) */
Quaternion quat_mul(Quaternion q1, Quaternion q2);

/** Normalize quaternion to unit */
Quaternion quat_normalize(Quaternion q);

/** Rotate a vector by unit quaternion: v' = q * v * q* */
vec3 quat_rotate_vec(Quaternion q, vec3 v);


/* ============================================================================
 * ODE Derivative Function Pointer Type
 * ============================================================================ */

/**
 * Generic ODE right-hand-side: dy/dt = f(t, y, params)
 * Used for rigid body equation integrators.
 */
typedef void (*ode_derivative_fn)(
    double t, const double y[], double dydt[], const void *params
);

#endif /* RIGID_BODY_TYPES_H */
