/**
 * @file hamiltonian_poisson.h
 * @brief Poisson bracket algebra and canonical invariants.
 *
 * Reference: Goldstein Ch.9, Arnold Ch.8, Marsden & Ratiu Ch.10
 *
 * L2 Core Concepts: Poisson bracket as Lie algebra on C^inf(P)
 * L3 Mathematical Structures: Lie algebra, Jacobi identity, Casimir invariants
 * L4 Fundamental Laws: Constants of motion via {f, H} = 0
 */

#ifndef HAMILTONIAN_POISSON_H
#define HAMILTONIAN_POISSON_H

#include "hamiltonian_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * Poisson bracket definition
 * ================================================================ */

/** Compute Poisson bracket {f, g} at (q, p).
 *
 *  {f, g} = sum_i (df/dq_i * dg/dp_i - df/dp_i * dg/dq_i)
 *
 *  Equivalent to: grad_f^T · J · grad_g.
 *
 *  Properties:
 *    1. Antisymmetry:  {f, g} = -{g, f}
 *    2. Bilinearity:   {af+bg, h} = a{f, h} + b{g, h}
 *    3. Leibniz rule:  {fg, h} = f{g, h} + {f, h}g
 *    4. Jacobi identity: {f, {g, h}} + {g, {h, f}} + {h, {f, g}} = 0
 *
 *  Classical limit of commutator: lim_{hbar->0} (1/ihbar)[A, B] = {A, B}.
 *
 *  O(n^2) due to numerical gradient evaluations.
 */
double poisson_bracket(Observable f, Observable g,
                        const double *q, const double *p, size_t n,
                        double eps);

/** Poisson bracket for separable f(q), g(p): {f, g} = sum_i (df/dq_i)(dg/dp_i).
 *  Avoids numerical gradients when separability is known. */
double poisson_bracket_qp(Observable f_q, Observable g_p,
                           const double *q, const double *p, size_t n,
                           double eps);

/** Verify fundamental brackets: {q_i,q_j}=0, {p_i,p_j}=0, {q_i,p_j}=delta_ij.
 *  Returns maximum deviation across all three checks.
 *  This verifies that (q,p) form a Darboux chart. */
double verify_fundamental_brackets(size_t n,
                                    const double *q, const double *p,
                                    double *err_qq, double *err_pp,
                                    double *err_qp,
                                    double eps);

/** Test if f is a constant of motion: {f, H} ~ 0.
 *
 *  Noether's theorem: {f, H} = 0 iff f is conserved.
 *  Poisson's theorem: if {f, H} = 0 and {g, H} = 0, then {{f, g}, H} = 0.
 */
bool is_constant_of_motion(Observable f,
                            const HamiltonianSystem *sys,
                            const double *q, const double *p,
                            double eps, double tol);

/** Verify Jacobi identity: J(f,g,h) = {f,{g,h}} + {g,{h,f}} + {h,{f,g}} = 0.
 *  This is the defining property of a Poisson algebra. */
double verify_jacobi_identity(Observable f, Observable g, Observable h,
                               const double *q, const double *p, size_t n,
                               double eps);

/** Test if C is a Casimir function: {C, f} = 0 for ALL f.
 *  Casimirs are conserved for ANY Hamiltonian.
 *  Example: for so(3)*, C(L) = |L|^2. */
bool is_casimir_function(Observable C,
                          const double *q, const double *p, size_t n,
                          double eps, double tol, int n_test);

/* ================================================================
 * Angular momentum Poisson algebra
 * ================================================================ */

/** Angular momentum L = r x p for N particles.
 *  L = sum_a r_a x p_a  (3-component vector output). */
void angular_momentum_vector(size_t n_particles,
                              const double *q, const double *p,
                              double *L);

/** Verify SO(3) algebra: {L_i, L_j} = epsilon_{ijk} L_k.
 *  Returns max error across all 3 relations. */
double verify_so3_poisson_algebra(size_t n_particles,
                                    const double *q, const double *p,
                                    double eps, double tol);

/** L^2 = Lx^2 + Ly^2 + Lz^2 — Casimir of so(3) Poisson algebra. */
double angular_momentum_casimir(size_t n_particles,
                                 const double *q, const double *p);

/* ================================================================
 * Lie-Poisson bracket
 * ================================================================ */

/** (-) Lie-Poisson bracket on g*:
 *  {F, G}(mu) = -< mu, [grad_F(mu), grad_G(mu)] >
 *            = -sum_{i,j,k} mu_i c_{ij}^k (dF/dmu_j)(dG/dmu_k)
 *
 *  Physical examples: Euler rigid body, ideal fluids, Vlasov plasma.
 *  Reference: Marsden & Ratiu, Intro to Mechanics and Symmetry, Ch.10 */
double lie_poisson_bracket(double (*F)(const double *mu, size_t dim),
                            double (*G)(const double *mu, size_t dim),
                            const double *mu, size_t dim,
                            const double *structure_constants,
                            double eps);

/** Euler rigid body: dm/dt = m x (I^{-1} m).
 *  These are the Lie-Poisson equations on so(3)*. */
void euler_rigid_body_rhs(const double *m, const double *I_diag, double *dm);

/** Initialize SO(3) structure constants: c_{ij}^k = epsilon_{ijk}.
 *  Output: flattened 3D array [27] = c[i*9 + j*3 + k]. */
void so3_structure_constants(double *c_out);

#ifdef __cplusplus
}
#endif

#endif /* HAMILTONIAN_POISSON_H */
