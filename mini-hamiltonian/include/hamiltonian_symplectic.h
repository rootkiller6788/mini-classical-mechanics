/**
 * @file hamiltonian_symplectic.h
 * @brief Symplectic geometry and Hamiltonian flows.
 *
 * Reference: Arnold Ch.8, Marsden & Ratiu Ch.1-5, da Silva Lectures on
 *            Symplectic Geometry
 *
 * L2 Core Concepts: Symplectic manifold, symplectic form, Darboux theorem
 * L3 Mathematical Structures: Symplectic vector space, symplectic group Sp(2n)
 * L4 Fundamental Laws: Hamiltonian flow preserves symplectic form
 * L8 Advanced Topics: Moment map, Marsden-Weinstein reduction, Lie-Poisson
 *
 * University mapping:
 *   Berkeley PHYS 242 — Symplectic geometry in classical mechanics
 *   Cambridge Part III — Symplectic methods in mathematical physics
 *   Oxford CMT — Momentum maps and geometric mechanics
 */

#ifndef HAMILTONIAN_SYMPLECTIC_H
#define HAMILTONIAN_SYMPLECTIC_H

#include "hamiltonian_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * Symplectic form and symplectic vector space
 * ================================================================ */

/**
 * Evaluate the canonical symplectic form omega = sum_i dq_i ∧ dp_i.
 *
 * For two tangent vectors xi = (delta_q, delta_p) and eta = (delta_q', delta_p'):
 *   omega(xi, eta) = sum_i (delta_q_i * delta_p'_i - delta_p_i * delta_q'_i)
 *                 = xi^T · J · eta
 *
 * where J is the canonical symplectic matrix.
 *
 * Properties:
 *   — omega is non-degenerate: omega(xi, ·) = 0 => xi = 0
 *   — omega is closed: d(omega) = 0 (in exterior calculus)
 *   — omega is exact: omega = d(p·dq) = -d(q·dp)
 *
 * Darboux theorem: Every symplectic form locally looks like the canonical
 * form above in suitable coordinates.
 */
double symplectic_form_eval(size_t n,
                             const double *dq1, const double *dp1,
                             const double *dq2, const double *dp2);

/**
 * Compute the Hamiltonian vector field X_f = J · grad_f from a function f.
 *
 * Definition: i_{X_f} omega = -df
 * In Darboux coordinates: X_f = (df/dp, -df/dq)
 *
 * The vector field X_f generates a flow that preserves f:
 *   f(phi_t(x)) = f(x)  iff {f, f} = 0 (trivially true)
 *
 * For the Hamiltonian H, X_H generates the time evolution.
 */
void symplectic_gradient(Observable f,
                          const double *q, const double *p, size_t n,
                          double eps,
                          double *X_f);

/**
 * Compute the symplectic form omega for the standard Darboux coordinates.
 *
 * omega is represented as the 2n×2n skew-symmetric matrix J.
 * In exterior algebra: omega = sum_{i=1}^n dq_i ∧ dp_i.
 * In matrix form: J = [[0, I_n], [-I_n, 0]].
 */
void symplectic_form_matrix(size_t n, double *J_out);

/* ================================================================
 * Symplectic group and linear canonical transformations
 * ================================================================ */

/**
 * The symplectic group Sp(2n, R) consists of matrices M satisfying M^T J M = J.
 *
 * Properties:
 *   — det(M) = 1 (Liouville's theorem for linear maps)
 *   — If lambda is an eigenvalue, so are 1/lambda, lambda*, 1/lambda*
 *   — dim(Sp(2n, R)) = n(2n+1)
 *
 * Sp(2, R) = SL(2, R), the special linear group in 2D.
 */

/**
 * Check if a 2n×2n matrix is symplectic: M^T J M = J.
 */
bool is_symplectic_matrix_direct(const double *M, size_t n, double tol);

/**
 * Generate a random symplectic matrix for testing.
 *
 * Uses the Euler decomposition: any symplectic matrix can be written as
 * M = O · D · O' where O, O' are orthogonal symplectic and D is diagonal.
 *
 * Reference: Dopico & Johnson (2009), SIAM J. Matrix Anal.
 */
void random_symplectic_matrix(size_t n, double *M);

/**
 * Exponential of a Hamiltonian matrix.
 *
 * If A is a Hamiltonian matrix (A^T J + J A = 0), then exp(tA) is
 * a one-parameter subgroup of symplectic matrices.
 *
 * For a quadratic Hamiltonian H = (1/2)·z^T·S·z where z = [q; p],
 * the flow is z(t) = exp(t·J·S)·z(0).
 */
void hamiltonian_matrix_exponential(const double *S, size_t n, double t,
                                     double *exp_tJS);

/* ================================================================
 * Moment map (momentum map)
 * ================================================================ */

/**
 * Compute the moment map for a symplectic group action.
 *
 * For a Lie group G acting on phase space by symplectic transformations,
 * the momentum map J: P → g* satisfies:
 *
 *   d<J, xi> = omega(Xi_xi, ·)
 *
 * where Xi_xi is the infinitesimal generator of xi ∈ g.
 *
 * Examples:
 *   G = SO(3) (rotations):  J = r × p  (angular momentum)
 *   G = R^3 (translations): J = p      (linear momentum)
 *   G = R     (time translations): J = -H (energy, up to sign)
 *
 * Noether's theorem: J is conserved along the flow of any G-invariant
 * Hamiltonian.
 *
 * Reference: Marsden & Ratiu, Introduction to Mechanics and Symmetry, Ch.11
 */
void momentum_map_so3(size_t n_particles,
                       const double *q, const double *p,
                       double *J);

/**
 * Check equivariance of a momentum map.
 *
 * A momentum map is equivariant if:
 *   J(g·x) = Ad*_{g^{-1}} J(x)   for all g ∈ G
 *
 * Infinitesimally: {J_xi, J_eta} = J_{[xi, eta]}
 *
 * This makes (J, {·,·}) a Lie algebra homomorphism from g to C^inf(P).
 */
bool check_momentum_map_homomorphism(const MomentumMap *mmap,
                                      const double *q, const double *p,
                                      double eps, double tol);

/* ================================================================
 * Marsden-Weinstein symplectic reduction
 * ================================================================ */

/**
 * Compute the reduced phase space dimension after symplectic reduction.
 *
 * Marsden-Weinstein-Meyer theorem (1974):
 *   If G acts freely and properly on P, and mu ∈ g* is a regular value
 *   of the momentum map J, then the reduced space
 *     P_mu = J^{-1}(mu) / G_mu
 *   is a symplectic manifold of dimension dim(P) - dim(G) - dim(G_mu).
 *
 * where G_mu is the isotropy subgroup of mu under the coadjoint action.
 *
 * @param dim_phase     Dimension of original phase space (2n)
 * @param dim_group     Dimension of symmetry group G
 * @param dim_isotropy  Dimension of isotropy subgroup G_mu
 * @return              Dimension of reduced symplectic manifold
 */
int reduced_phase_space_dimension(int dim_phase, int dim_group, int dim_isotropy);

/* ================================================================
 * Lie-Poisson bracket on g*
 * ================================================================ */

/**
 * Compute the (-) Lie-Poisson bracket {F, G}(mu).
 *
 * This is the reduced Poisson bracket on g* obtained by Marsden-Weinstein
 * reduction of T*G by the left G-action.
 *
 * The equations of motion for H: g* → R are:
 *   dmu/dt = ad*_{dH(mu)} mu
 *
 * For SO(3)*, this gives the Euler rigid body equations:
 *   dm_i/dt = epsilon_{ijk} m_j (I^{-1} m)_k
 */
double lie_poisson_bracket_eval(
    double (*F)(const double *mu, size_t dim),
    double (*G)(const double *mu, size_t dim),
    const double *mu, size_t dim,
    const double *struct_const, double eps);

/**
 * SO(3) Lie-Poisson equations of motion for the rigid body.
 *
 * d/dt [m1, m2, m3] = [m2*m3*(1/I3-1/I2), m3*m1*(1/I1-1/I3), m1*m2*(1/I2-1/I1)]
 *
 * These are the Euler equations for a free rigid body.
 * The solutions lie on the intersection of the sphere |m|^2 = const
 * (Casimir surface) and the ellipsoid H(m) = const (energy surface).
 */
void so3_lie_poisson_rhs(const double *m, const double *I_diag, double *dm);

#ifdef __cplusplus
}
#endif

#endif /* HAMILTONIAN_SYMPLECTIC_H */
