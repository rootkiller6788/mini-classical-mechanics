/**
 * rigid_body_inertia.h — Inertia Tensor: Computation, Principal Axes, Standard Shapes
 *
 * Reference:
 *   Goldstein §5.1-5.3: Moment of inertia tensor, principal axes, parallel axis theorem
 *   Landau §32: Rigid body inertia tensor
 *   MIT 8.012 Lecture 24-25: Moment of inertia
 *
 * The inertia tensor I encodes mass distribution. For a rigid body with N particles:
 *   I_ij = Σ m_k (r_k² δ_ij - r_{k,i} r_{k,j})     (Goldstein 5.3)
 *
 * Key theorems:
 *   Parallel Axis (Steiner): I_new = I_cm + M (d² 𝟙 - d⊗d)   (Goldstein 5.18)
 *   Principal Axis Theorem:  I is diagonalizable with orthonormal eigenvectors
 */

#ifndef RIGID_BODY_INERTIA_H
#define RIGID_BODY_INERTIA_H

#include "rigid_body_types.h"

/* ============================================================================
 * Inertia Tensor from Discrete Particles
 * ============================================================================ */

/**
 * Compute inertia tensor from a set of point masses relative to their center of mass.
 *
 * @param n_particles  Number of particles
 * @param masses       Array of masses [n_particles]
 * @param positions    Array of 3D positions relative to CM [n_particles]
 * @param out          Output inertia tensor (in principal-axis-aligned or arbitrary frame)
 *
 * Formula (Goldstein 5.3):
 *   I_xx = Σ m(y²+z²), I_yy = Σ m(x²+z²), I_zz = Σ m(x²+y²)
 *   I_xy = -Σ mxy,     I_xz = -Σ mxz,     I_yz = -Σ myz
 *
 * Complexity: O(N) in number of particles
 */
void inertia_tensor_from_particles(
    int n_particles, const double masses[], const vec3 positions[],
    InertiaTensor *out
);

/**
 * Compute inertia tensor for a continuous body by numerical integration over a grid.
 *
 * For bodies defined by a density function ρ(r):
 *   I_ij = ∫_V ρ(r) (r² δ_ij - r_i r_j) d³r
 *
 * @param density_func  ρ(x,y,z) — density at point
 * @param xlim,ylim,zlim  Integration bounds [{xmin,xmax}, {ymin,ymax}, {zmin,zmax}]
 * @param n_grid        Number of grid points per dimension
 * @param out           Output inertia tensor
 *
 * Uses midpoint rule quadrature. Complexity: O(n_grid³)
 */
void inertia_tensor_numerical(
    double (*density_func)(double x, double y, double z),
    const double xlim[2], const double ylim[2], const double zlim[2],
    int n_grid, InertiaTensor *out
);


/* ============================================================================
 * Parallel Axis Theorem (Steiner's Theorem)
 *
 * Goldstein (5.18): If I_cm is the inertia tensor about the center of mass,
 * and d is the vector from CM to a parallel axis through point P,
 * then: I_P = I_cm + M (d² 𝟙 - d⊗d)
 * ============================================================================ */

/**
 * Translate inertia tensor from CM to a parallel point.
 * @param I_cm  Inertia tensor about center of mass
 * @param M     Total mass
 * @param d     Displacement vector from CM to new reference point
 * @param out   Inertia tensor about new reference point
 */
void parallel_axis_theorem(const InertiaTensor *I_cm, double M, vec3 d, InertiaTensor *out);

/**
 * Inverse: translate inertia tensor from point P back to CM.
 * @param I_p   Inertia tensor about point P
 * @param M     Total mass
 * @param d     Displacement from CM to P (same direction)
 * @param out   Inertia tensor about CM
 */
void parallel_axis_inverse(const InertiaTensor *I_p, double M, vec3 d, InertiaTensor *out);


/* ============================================================================
 * Principal Axes Decomposition
 *
 * Goldstein (5.4)-(5.6): The inertia tensor I can be diagonalized:
 *   I_principal = R^T · I · R = diag(I1, I2, I3)
 * where R is an orthogonal matrix whose columns are the principal axes (eigenvectors).
 * ============================================================================ */

/**
 * Compute principal axes via Jacobi eigenvalue algorithm (robust for 3×3 symmetric).
 * @param I   Input inertia tensor
 * @param out Output principal axes (moments sorted descending: I1 ≥ I2 ≥ I3)
 */
void principal_axes_decompose(const InertiaTensor *I, PrincipalAxes *out);

/**
 * Compute principal moments only (without axes), sorted descending.
 * @param I       Input inertia tensor
 * @param moments Output array of 3 moments [I1, I2, I3]
 */
void principal_moments_only(const InertiaTensor *I, double moments[3]);

/**
 * Check if inertia tensor is already diagonal (i.e., expressed in principal axes).
 * @param I   Inertia tensor to check
 * @param tol Tolerance for off-diagonal elements
 * @return 1 if diagonal, 0 otherwise
 */
int is_inertia_diagonal(const InertiaTensor *I, double tol);


/* ============================================================================
 * Inertia Tensor of Standard Homogeneous Shapes (CM frame)
 *
 * All formulas from Goldstein Table 5.1 / standard integration results.
 * M = mass, all results in CM frame.
 * ============================================================================ */

/**
 * Solid sphere of radius R: I_xx = I_yy = I_zz = (2/5) M R²
 * All products of inertia zero by symmetry.
 */
InertiaTensor sphere_inertia(double M, double R);

/**
 * Spherical shell of radius R: I_xx = I_yy = I_zz = (2/3) M R²
 */
InertiaTensor spherical_shell_inertia(double M, double R);

/**
 * Solid cylinder, symmetry axis = z:
 *   I_zz = (1/2) M R²        (about symmetry axis)
 *   I_xx = I_yy = (1/12) M (3R² + H²)  (about transverse axes through CM)
 */
InertiaTensor cylinder_inertia(double M, double R, double H);

/**
 * Hollow cylinder (thin-walled tube), symmetry axis = z:
 *   I_zz = M R²
 *   I_xx = I_yy = (1/2) M R² + (1/12) M H²
 */
InertiaTensor hollow_cylinder_inertia(double M, double R, double H);

/**
 * Rectangular cuboid (a × b × c along x, y, z):
 *   I_xx = (1/12) M (b² + c²)
 *   I_yy = (1/12) M (a² + c²)
 *   I_zz = (1/12) M (a² + b²)
 */
InertiaTensor cuboid_inertia(double M, double a, double b, double c);

/**
 * Thin uniform rod of length L (along z-axis through CM):
 *   I_xx = I_yy = (1/12) M L²
 *   I_zz = 0 (idealized, neglecting thickness)
 */
InertiaTensor rod_inertia(double M, double L);

/**
 * Thin circular disk of radius R (in xy-plane, z = symmetry axis):
 *   I_zz = (1/2) M R²       (about symmetry axis)
 *   I_xx = I_yy = (1/4) M R² (diametral axes through CM)
 */
InertiaTensor disk_inertia(double M, double R);

/**
 * Solid ellipsoid with semi-axes (a, b, c) along (x, y, z):
 *   I_xx = (1/5) M (b² + c²)
 *   I_yy = (1/5) M (a² + c²)
 *   I_zz = (1/5) M (a² + b²)
 */
InertiaTensor ellipsoid_inertia(double M, double a, double b, double c);

/**
 * Solid cone, symmetry axis = z, base radius R, height H, apex at z=0, base at z=H:
 * CM at z = 3H/4 from apex.
 * About CM:
 *   I_zz = (3/10) M R²
 *   I_xx = I_yy = (3/20) M R² + (3/80) M H²
 */
InertiaTensor solid_cone_inertia(double M, double R, double H);

/**
 * Solid torus, major radius R, minor (tube) radius r:
 * Symmetry axis = z.
 *   I_zz = M (R² + (3/4) r²)
 *   I_xx = I_yy = (1/2) M R² + (5/8) M r²
 */
InertiaTensor torus_inertia(double M, double R, double r);

/**
 * Composite body: sum inertia tensors about a common reference point.
 * For a body composed of parts with known inertias (all about same point):
 *   I_total = Σ I_k
 */
void inertia_compose(int n_parts, const InertiaTensor parts[], InertiaTensor *out);


/* ============================================================================
 * Inertia Ellipsoid (Poinsot Construction)
 *
 * Goldstein (5.11): The inertia ellipsoid is defined by ρᵀ · I · ρ = 1.
 * Its semi-axes are (1/√I1, 1/√I2, 1/√I3) along principal axes.
 * ============================================================================ */

/**
 * Compute semi-axes of the inertia ellipsoid.
 * @param I  Inertia tensor
 * @param a  Output: semi-axis along 1st principal axis (1/√I1)
 * @param b  Output: semi-axis along 2nd principal axis (1/√I2)
 * @param c  Output: semi-axis along 3rd principal axis (1/√I3)
 */
void inertia_ellipsoid_semiaxes(const InertiaTensor *I, double *a, double *b, double *c);

/**
 * Generate sample points on the inertia ellipsoid surface for visualization.
 * @param I         Inertia tensor
 * @param n_theta   Number of polar angle samples
 * @param n_phi     Number of azimuthal angle samples
 * @param points    Output array of size n_theta*n_phi (caller allocates)
 */
void inertia_ellipsoid_points(const InertiaTensor *I, int n_theta, int n_phi, vec3 points[]);

/**
 * Evaluate the inertia ellipsoid quadratic form at a point.
 * f(ρ) = ρᵀ · I · ρ
 * Returns 1.0 if ρ lies on the inertia ellipsoid surface.
 */
double inertia_ellipsoid_eval(const InertiaTensor *I, vec3 rho);

/**
 * Compute the radius of gyration about an arbitrary axis through CM.
 * k = √(I_axis / M) where I_axis = n̂ᵀ · I · n̂
 * @param I    Inertia tensor about CM
 * @param M    Total mass
 * @param axis Unit vector defining rotation axis
 * @return     Radius of gyration k
 */
double radius_of_gyration(const InertiaTensor *I, double M, vec3 axis);

#endif /* RIGID_BODY_INERTIA_H */
