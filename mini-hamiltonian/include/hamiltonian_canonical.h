/**
 * @file hamiltonian_canonical.h
 * @brief Canonical transformations and generating functions.
 *
 * Reference: Goldstein Ch.9, Arnold Ch.9, Landau-Lifshitz Vol.1 Ch.7
 *
 * L2 Core Concepts: Canonical transformations as symplectic diffeomorphisms
 * L3 Mathematical Structures: Symplectic group Sp(2n, R), generating functions
 * L4 Fundamental Laws: Invariance of Hamilton's equations under canonical transformations
 *
 * University mapping:
 *   MIT 8.012 — Canonical transformations, generating functions
 *   Stanford PHYSICS 230 — Symplectic geometry, Lie transform method
 *   Caltech Ph 106 — Canonical formalism and perturbation theory
 */

#ifndef HAMILTONIAN_CANONICAL_H
#define HAMILTONIAN_CANONICAL_H

#include "hamiltonian_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * Canonical transformation construction
 * ================================================================ */

/**
 * Create the identity canonical transform: Q = q, P = p.
 *
 * Generating function: F2(q, P) = q·P = sum_i q_i P_i.
 * Then: p = dF2/dq = P, Q = dF2/dP = q  -> identity.
 */
CanonicalTransform identity_transform_create(size_t n);

/**
 * Create a scale transform: Q = mu*q, P = nu*p with mu*nu = 1.
 *
 * This preserves the symplectic form because:
 *   dQ ^ dP = (mu dq) ^ (nu dp) = mu*nu * dq ^ dp = dq ^ dp.
 *
 * Example: For a harmonic oscillator, choosing mu = sqrt(m*omega)
 * converts H = p^2/(2m) + (1/2)m*omega^2*q^2 to H = (omega/2)(P^2+Q^2).
 */
CanonicalTransform scale_transform_create(double mu, size_t n);

/**
 * Create an exchange transform: Q = p, P = -q.
 *
 * This swaps coordinates and momenta (up to sign).
 * The generating function is F1(q, Q) = q·Q.
 * Then: p = dF1/dq = Q, P = -dF1/dQ = -q.
 *
 * This transform reveals the symmetry between q and p in Hamilton's equations.
 */
CanonicalTransform exchange_transform_create(size_t n);

/**
 * Create a point transform: Q = f(q), P = (df/dq)^{-T}·p.
 *
 * This lifts a coordinate transformation Q = f(q) to a canonical
 * transformation on the full phase space.
 * If f is a diffeomorphism, the lifted transform preserves the
 * canonical Poisson brackets.
 *
 * @param f       Coordinate transformation Q = f(q)
 * @param jac_f   Jacobian matrix df/dq [n×n], row-major
 * @param inv_jac Inverse Jacobian (df/dq)^{-1} [n×n], row-major
 */
CanonicalTransform point_transform_create(size_t n,
    void (*f)(const double *q, double *Q),
    void (*jac_f)(const double *q, double *J_out),
    void (*inv_jac)(const double *q, double *Jinv_out));

/**
 * Create a canonical transform from a type-2 generating function F2(q, P).
 *
 * Implicit equations:
 *   p = dF2/dq(q, P)    — solve for P given (q, p)
 *   Q = dF2/dP(q, P)    — then compute Q
 *
 * This is the most commonly used generator type.
 *
 * @param F2        Generating function F2(q, P) -> scalar
 * @param grad_q_F2 Gradient with respect to q: dF2/dq [n]
 * @param grad_P_F2 Gradient with respect to P: dF2/dP [n]
 * @param solve_P   Nonlinear solver: given (q, p), find P s.t. grad_q_F2(q, P) = p
 */
CanonicalTransform f2_transform_create(size_t n,
    double (*F2)(const double *q, const double *P, size_t n),
    void (*grad_q_F2)(const double *q, const double *P, double *out),
    void (*grad_P_F2)(const double *q, const double *P, double *out),
    void (*solve_P)(const double *q, const double *p, double *P_out));

/* ================================================================
 * Symplectic condition verification
 * ================================================================ */

/**
 * Verify that a transformation is canonical via fundamental Poisson brackets.
 *
 * Check: {Q_i, Q_j}=0, {P_i, P_j}=0, {Q_i, P_j}=delta_ij
 * where the brackets are computed in the original (q, p) coordinates.
 *
 * Theorem: A transformation (q, p) -> (Q, P) is canonical iff
 *   {Q_i, Q_j}_{q,p} = 0, {P_i, P_j}_{q,p} = 0, {Q_i, P_j}_{q,p} = delta_ij.
 *
 * @param tol   Tolerance for bracket deviations
 * @param err   Output: maximum absolute deviation across all checks
 * @return      true if canonical within tolerance
 */
bool verify_canonical_transform(const CanonicalTransform *ct,
                                 const double *q, const double *p,
                                 double eps, double tol,
                                 double *err);

/**
 * Verify canonical condition via symplectic Jacobian.
 *
 * Let M = d(Q, P)/d(q, p) be the Jacobian of the transformation.
 * For a canonical transformation: M^T · J · M = J.
 *
 * This is the defining condition for a symplectic matrix M ∈ Sp(2n, R).
 * The set of all such M forms the symplectic group.
 *
 * @return  ||M^T J M - J||_F  (Frobenius norm of deviation)
 */
double symplectic_jacobian_check(const CanonicalTransform *ct,
                                  const double *q, const double *p,
                                  double eps);

/**
 * Compute the Jacobian matrix M = d(Q,P)/d(q,p) numerically.
 *
 * M is a 2n×2n matrix: [dQ/dq  dQ/dp]
 *                      [dP/dq  dP/dp]
 *
 * @param M  Output: Jacobian matrix, 2n×2n, row-major
 */
void compute_canonical_jacobian(const CanonicalTransform *ct,
                                 const double *q, const double *p,
                                 double eps,
                                 double *M);

/* ================================================================
 * Transform operations
 * ================================================================ */

/**
 * Compose two canonical transforms: CT3 = CT2 ∘ CT1.
 *
 * If CT1: (q,p) -> (Q',P') and CT2: (Q',P') -> (Q,P),
 * then CT3: (q,p) -> (Q,P) with CT3 = CT2 ∘ CT1.
 *
 * Theorem: The composition of canonical transforms is canonical.
 * Proof: M3 = M2·M1, and if M1^T J M1 = J and M2^T J M2 = J,
 * then M3^T J M3 = M1^T M2^T J M2 M1 = M1^T J M1 = J.
 */
CanonicalTransform compose_canonical_transforms(const CanonicalTransform *ct1,
                                                  const CanonicalTransform *ct2);

/**
 * Invert a canonical transform.
 *
 * If ct: (q,p) -> (Q,P) then ct_inv: (Q,P) -> (q,p).
 * The inverse of a canonical transform is also canonical since
 * M^{-1} T J M^{-1} = J (from M^T J M = J, multiply by M^{-T} on left
 * and M^{-1} on right).
 */
CanonicalTransform inverse_canonical_transform(const CanonicalTransform *ct);

/**
 * Free resources held by a canonical transform.
 */
void canonical_transform_free(CanonicalTransform *ct);

/* ================================================================
 * Symplectic matrix utilities
 * ================================================================ */

/**
 * Allocate the canonical symplectic matrix J_2n = [[0, I], [-I, 0]].
 *
 * J has dimension 2n×2n and satisfies:
 *   J^2 = -I, J^T = -J, J^{-1} = -J, det(J) = 1.
 */
SymplecticMatrix symplectic_matrix_alloc(size_t n);

/**
 * Free symplectic matrix.
 */
void symplectic_matrix_free(SymplecticMatrix *J);

/**
 * Test if a matrix M ∈ R^{2n×2n} is symplectic: M^T J M = J.
 *
 * @param M  Matrix data, 2n×2n, row-major
 * @param tol Tolerance for equality check
 */
bool is_symplectic_matrix(const double *M, size_t n, double tol);

/**
 * Compute the symplectic eigenvalues of a matrix.
 *
 * Symplectic eigenvalues come in reciprocal pairs (lambda, 1/lambda).
 * This property arises from the symplectic condition M^T J M = J.
 *
 * @param M         Input matrix, 2n×2n, row-major
 * @param eig_real  Output: real parts of eigenvalues [2n]
 * @param eig_imag  Output: imaginary parts of eigenvalues [2n]
 */
void symplectic_eigenvalues(const double *M, size_t n,
                             double *eig_real, double *eig_imag);

/**
 * Compute the linear stability of a fixed point via symplectic eigenvalues.
 *
 * A fixed point is:
 *   — Elliptic if all eigenvalues are on the unit circle (stable)
 *   — Hyperbolic if any eigenvalue has |lambda| != 1 (unstable)
 *
 * For Hamiltonian systems, Lyapunov's theorem: if all eigenvalues are
 * purely imaginary and distinct, the fixed point is surrounded by
 * invariant tori (KAM stability).
 */
typedef enum { STABILITY_ELLIPTIC, STABILITY_HYPERBOLIC, STABILITY_MIXED } StabilityType;

StabilityType classify_fixed_point_stability(const double *M, size_t n, double tol);

#ifdef __cplusplus
}
#endif

#endif /* HAMILTONIAN_CANONICAL_H */
