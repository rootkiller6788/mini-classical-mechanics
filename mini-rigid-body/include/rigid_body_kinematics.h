/**
 * rigid_body_kinematics.h — Rigid Body Kinematics: Euler Angles, Rotation Matrices, Quaternions
 *
 * Reference:
 *   Goldstein §4.8-4.9: Euler angles, rotation matrices, angular velocity
 *   Landau §35: Eulerian angles
 *   MIT 8.012 Lecture 21-22: Rotational kinematics
 *
 * Three parameterizations of SO(3) are provided:
 *   1. Euler angles (ZXZ convention, Goldstein 4.46)
 *   2. Rotation matrices (3×3 orthogonal, det=+1)
 *   3. Quaternions (SU(2) double cover, singularity-free)
 *
 * Key relationships:
 *   - Euler → Matrix:  R(φ,θ,ψ) = R_z(φ) · R_x(θ) · R_z(ψ)
 *   - Matrix → Euler:  inverse trigonometric extraction
 *   - ω (body) ↔ Euler rates: via kinematic mapping
 *   - dR/dt = R · [ω]_×   (so(3) Lie algebra)
 */

#ifndef RIGID_BODY_KINEMATICS_H
#define RIGID_BODY_KINEMATICS_H

#include "rigid_body_types.h"

/* ============================================================================
 * Euler Angles ↔ Rotation Matrix (ZXZ Convention)
 * ============================================================================ */

/**
 * Convert Euler angles to 3×3 rotation matrix.
 *
 * ZXZ convention (Goldstein 4.46):
 *   R(φ,θ,ψ) = R_z(φ) · R_x(θ) · R_z(ψ)
 *
 * Matrix elements:
 *   R₁₁ =  cosφ cosψ - sinφ cosθ sinψ
 *   R₁₂ = -cosφ sinψ - sinφ cosθ cosψ
 *   R₁₃ =  sinφ sinθ
 *   R₂₁ =  sinφ cosψ + cosφ cosθ sinψ
 *   R₂₂ = -sinφ sinψ + cosφ cosθ cosψ
 *   R₂₃ = -cosφ sinθ
 *   R₃₁ =  sinθ sinψ
 *   R₃₂ =  sinθ cosψ
 *   R₃₃ =  cosθ
 *
 * @param phi    Precession angle
 * @param theta  Nutation angle
 * @param psi    Spin angle
 * @param R      Output 3×3 rotation matrix (orthogonal, det=+1)
 */
void euler_to_rotation_matrix(double phi, double theta, double psi, mat3 *R);

/**
 * Convert 3×3 rotation matrix to Euler angles.
 *
 * Goldstein (4.47): Extract (φ, θ, ψ) from matrix elements.
 * Singularities at θ = 0 and θ = π (gimbal lock).
 *
 * @param R   Input rotation matrix (must be orthogonal)
 * @param ea  Output Euler angles
 * @return    0 on success, -1 if matrix is not a valid rotation
 */
int rotation_matrix_to_euler(const mat3 *R, EulerAngles *ea);

/**
 * Alternative Euler angle conventions:
 * - ZYZ convention (common in quantum mechanics)
 * - XYZ convention (Tait-Bryan / aerospace)
 */

/** ZYZ Euler angles → rotation matrix: R_z(α) · R_y(β) · R_z(γ) */
void euler_zyz_to_matrix(double alpha, double beta, double gamma, mat3 *R);

/** Rotation matrix → ZYZ Euler angles */
int matrix_to_euler_zyz(const mat3 *R, double *alpha, double *beta, double *gamma);

/** XYZ Tait-Bryan (roll-pitch-yaw) → rotation matrix: R_x(roll) · R_y(pitch) · R_z(yaw) */
void tait_bryan_to_matrix(double roll, double pitch, double yaw, mat3 *R);

/** Rotation matrix → Tait-Bryan angles */
int matrix_to_tait_bryan(const mat3 *R, double *roll, double *pitch, double *yaw);


/* ============================================================================
 * Euler Angle Rates ↔ Angular Velocity
 * ============================================================================ */

/**
 * Convert Euler angle rates to body-frame angular velocity.
 *
 * Goldstein (4.86-4.88):
 *   ω_x =  φ̇ sinθ sinψ + θ̇ cosψ
 *   ω_y =  φ̇ sinθ cosψ - θ̇ sinψ
 *   ω_z =  φ̇ cosθ + ψ̇
 *
 * @param phi_dot    Rate of precession
 * @param theta_dot  Rate of nutation
 * @param psi_dot    Rate of spin
 * @param theta      Current nutation angle
 * @param psi        Current spin angle
 * @param omega      Output body-frame angular velocity [ωx, ωy, ωz]
 */
void euler_rates_to_body_omega(double phi_dot, double theta_dot, double psi_dot,
                               double theta, double psi, double omega[3]);

/**
 * Convert body-frame angular velocity to Euler angle rates.
 *
 * Inverse of the above (Goldstein 4.87 inverted):
 *   φ̇ = (ω_x sinψ + ω_y cosψ) / sinθ
 *   θ̇ =  ω_x cosψ - ω_y sinψ
 *   ψ̇ =  ω_z - φ̇ cosθ
 *
 * @param omega     Body-frame angular velocity
 * @param theta     Current nutation angle
 * @param psi       Current spin angle
 * @param phi_dot   Output: precession rate
 * @param theta_dot Output: nutation rate
 * @param psi_dot   Output: spin rate
 *
 * Note: Singular at θ≈0 (gimbal lock). Uses regularization near θ=0.
 */
void body_omega_to_euler_rates(const double omega[3], double theta, double psi,
                               double *phi_dot, double *theta_dot, double *psi_dot);

/**
 * Compute Euler angle ODE right-hand side for use with integrators.
 * dy/dt where y = (φ, θ, ψ).
 *
 * @param t        Time (unused, system is autonomous)
 * @param y        Current state [φ, θ, ψ]
 * @param dydt     Output [dφ/dt, dθ/dt, dψ/dt]
 * @param params   Pointer to angular velocity array ω[3]
 */
void euler_kinematic_ode(double t, const double y[], double dydt[], const void *params);


/* ============================================================================
 * Rotation Matrix — Angular Velocity Differential Relation
 * ============================================================================ */

/**
 * Compute time derivative of rotation matrix: dR/dt = R · [ω]_×
 *
 * This is the kinematic equation for SO(3):
 *   dR/dt = R · Ω   where Ω is the skew-symmetric matrix of ω.
 *
 * @param R    Current rotation matrix
 * @param omega Body-frame angular velocity
 * @param dR   Output: time derivative of R (3×3 matrix)
 */
void rotation_matrix_derivative(const mat3 *R, const double omega[3], mat3 *dR);

/**
 * Build the skew-symmetric (cross-product) matrix of a vector.
 *
 * [ω]_× = [[0, -ωz, ωy],
 *           [ωz, 0, -ωx],
 *           [-ωy, ωx, 0]]
 *
 * This matrix is in so(3), the Lie algebra of SO(3).
 * For any vector v: [ω]_× · v = ω × v
 *
 * @param omega  Input vector
 * @param skew   Output 3×3 skew-symmetric matrix
 */
void cross_product_matrix(const double omega[3], mat3 *skew);

/**
 * Recover vector from skew-symmetric matrix (inverse of cross_product_matrix).
 * @param skew  Skew-symmetric matrix
 * @param omega Output vector
 */
void unskew_matrix(const mat3 *skew, double omega[3]);


/* ============================================================================
 * Quaternion Operations (SU(2) Parameterization of Rotations)
 * ============================================================================ */

/**
 * Convert unit quaternion to rotation matrix.
 *
 * For q = (w, x, y, z) with |q|=1:
 *   R = [[1-2(y²+z²), 2(xy-wz),   2(xz+wy)],
 *        [2(xy+wz),   1-2(x²+z²), 2(yz-wx)],
 *        [2(xz-wy),   2(yz+wx),   1-2(x²+y²)]]
 */
void quaternion_to_rotation_matrix(Quaternion q, mat3 *R);

/**
 * Convert rotation matrix to unit quaternion.
 *
 * Shepperd's algorithm (numerically stable).
 * @param R  Rotation matrix
 * @param q  Output quaternion (unit)
 * @return   0 on success
 */
int rotation_matrix_to_quaternion(const mat3 *R, Quaternion *q);

/**
 * Convert Euler angles (ZXZ) to quaternion.
 * q = R_z(φ) · R_x(θ) · R_z(ψ) expressed as quaternion product.
 */
void euler_to_quaternion(double phi, double theta, double psi, Quaternion *q);

/**
 * Convert quaternion to Euler angles (ZXZ).
 */
void quaternion_to_euler(Quaternion q, EulerAngles *ea);

/**
 * Spherical linear interpolation (SLERP) between two unit quaternions.
 *
 * For rotations q0 and q1, the shortest-path interpolation at fraction t∈[0,1]:
 *   SLERP(q0, q1, t) = q0 · (q0*·q1)^t
 *                    = (sin((1-t)Ω) q0 + sin(tΩ) q1) / sin(Ω)
 * where Ω = acos(Re(q0* · q1)).
 *
 * @param q0  Start quaternion (unit)
 * @param q1  End quaternion (unit)
 * @param t   Interpolation parameter [0, 1]
 * @param out Output interpolated quaternion
 */
void quat_slerp(Quaternion q0, Quaternion q1, double t, Quaternion *out);

/**
 * Quaternion angular velocity: dq/dt = ½ ω_q · q
 * where ω_q = (0, ωx, ωy, ωz) is the pure quaternion of angular velocity.
 *
 * @param q       Current quaternion
 * @param omega   Body-frame angular velocity
 * @param dq      Output: time derivative dq/dt
 */
void quaternion_derivative(Quaternion q, const double omega[3], Quaternion *dq);

/**
 * Exponential map from so(3) to SU(2): q = exp(½ θ n̂)
 * where θ n̂ is the rotation vector (angle * unit axis).
 *
 * @param axis_angle  Rotation vector (θ * n̂)
 * @param q           Output unit quaternion
 */
void quat_exp(const double axis_angle[3], Quaternion *q);

/**
 * Logarithm map from SU(2) to so(3): extracts rotation vector from quaternion.
 *
 * @param q          Unit quaternion
 * @param axis_angle Output: rotation vector θ n̂
 */
void quat_log(Quaternion q, double axis_angle[3]);


/* ============================================================================
 * Basic Rotation Matrices
 * ============================================================================ */

/** Rotation about x-axis by angle (radians) */
void rotation_x(double angle, mat3 *R);

/** Rotation about y-axis by angle (radians) */
void rotation_y(double angle, mat3 *R);

/** Rotation about z-axis by angle (radians) */
void rotation_z(double angle, mat3 *R);

/**
 * Rodrigues' rotation formula: rotation about arbitrary axis by angle.
 *
 * R(n̂, θ) = I + sinθ [n̂]_× + (1-cosθ) [n̂]_×²
 *
 * @param axis   Rotation axis (will be normalized)
 * @param angle  Rotation angle (radians)
 * @param R      Output rotation matrix
 */
void rotation_axis_angle(vec3 axis, double angle, mat3 *R);

/**
 * Extract axis and angle from rotation matrix (inverse of Rodrigues).
 *
 * @param R      Rotation matrix
 * @param axis   Output: rotation axis (unit vector)
 * @param angle  Output: rotation angle (in [0, π])
 * @return       0 on success, -1 if trace(R) is outside [-1, 3]
 */
int rotation_to_axis_angle(const mat3 *R, vec3 *axis, double *angle);

/**
 * Compose two rotations: R12 = R1 · R2 (apply R2 first, then R1).
 */
void rotation_compose(const mat3 *R1, const mat3 *R2, mat3 *R12);

/**
 * Inverse rotation: R⁻¹ = Rᵀ (for orthogonal matrices).
 */
void rotation_inverse(const mat3 *R, mat3 *Rinv);

#endif /* RIGID_BODY_KINEMATICS_H */
