/**
 * @file hamiltonian_types.h
 * @brief Core Hamiltonian mechanics types and structures.
 *
 * Reference: Goldstein Classical Mechanics (3rd ed.) Ch.8-10
 *            Arnold Mathematical Methods of Classical Mechanics Ch.3-5
 *
 * L1 Definitions: All fundamental Hamiltonian mechanics data structures
 *   - PhasePoint, HamiltonianSystem, CanonicalTransform
 *   - SymplecticMatrix, PoissonAlgebra, ActionAngleSystem
 *   - EnergySurface, MomentumMap, LieAlgebraElement
 *
 * University mapping:
 *   MIT 8.012 Lecture 25 — Hamiltonian formulation
 *   Caltech Ph 106 Ch.6 — Canonical formalism
 *   Cambridge Part II Theoretical Physics — Symplectic geometry
 */

#ifndef HAMILTONIAN_TYPES_H
#define HAMILTONIAN_TYPES_H

#include <stddef.h>   /* size_t */
#include <stdbool.h>  /* bool, true, false */

#ifdef __cplusplus
extern "C" {
#endif

/* ──────────────────────────────────────────────────────────────
 * L1: Fundamental data structures
 * ────────────────────────────────────────────────────────────── */

/** Phase space point (q, p) in R^{2n}.
 *  Goldstone Theorem: Phase space dimension = 2 × degrees of freedom.
 *  The pair (q, p) forms a Darboux coordinate chart on the symplectic manifold.
 */
typedef struct {
    size_t n;          /**< Number of degrees of freedom */
    double *q;         /**< Generalized coordinates, length n */
    double *p;         /**< Conjugate momenta, length n */
    double t;          /**< Time at this point */
} PhasePoint;

/** Phase space trajectory: time-series of phase points.
 *  For a Hamiltonian flow phi_t, this records phi_t(q0, p0) at discrete times.
 */
typedef struct {
    size_t capacity;   /**< Allocated storage size */
    size_t length;     /**< Number of recorded points */
    double *times;     /**< Time array, length length */
    double **q_array;  /**< Coordinates history, length × n */
    double **p_array;  /**< Momentum history, length × n */
    double *H_values;  /**< Hamiltonian values along trajectory */
} PhaseTrajectory;

/** Hamiltonian system descriptor.
 *  H: R^{2n} → R  — the Hamiltonian function
 *  grad_H: R^{2n} → R^{2n} — gradient (dH/dq, dH/dp)
 *
 *  Hamilton's equations:
 *    dq/dt =  dH/dp     (n equations)
 *    dp/dt = -dH/dq     (n equations)
 *
 *  Equivalent to: d/dt [q; p] = J · grad_H  where J = [[0, I], [-I, 0]]
 */
typedef struct {
    size_t n;                    /**< Degrees of freedom */
    double (*H)(const double *q, const double *p, size_t n);
    void   (*grad_H)(const double *q, const double *p, size_t n,
                     double *dH_dq, double *dH_dp);
    bool   is_autonomous;        /**< True if H does not depend explicitly on t */
    char   *name;                /**< Optional system name */
} HamiltonianSystem;

/** Canonical transformation descriptor.
 *  A map (q, p) → (Q, P) that preserves the symplectic form:
 *    dQ wedge dP = dq wedge dp
 *
 *  Generating function types (Goldstein Ch.9):
 *    F1(q, Q, t):  p = dF1/dq,  P = -dF1/dQ
 *    F2(q, P, t):  p = dF2/dq,  Q =  dF2/dP
 *    F3(p, Q, t):  q = -dF3/dp, P = -dF3/dQ
 *    F4(p, P, t):  q = -dF4/dp, Q =  dF4/dP
 */
typedef enum {
    CT_GENERATOR_F1,   /**< Type-1: F1(q, Q) */
    CT_GENERATOR_F2,   /**< Type-2: F2(q, P) — most common */
    CT_GENERATOR_F3,   /**< Type-3: F3(p, Q) */
    CT_GENERATOR_F4    /**< Type-4: F4(p, P) */
} CanonicalGeneratorType;

typedef struct {
    size_t n;                    /**< Degrees of freedom */
    CanonicalGeneratorType type; /**< Generator type */
    void (*transform)(const double *q, const double *p, size_t n,
                      double *Q, double *P);
    void (*inverse)(const double *Q, const double *P, size_t n,
                    double *q, double *p);
    double (*generator)(const double *a, const double *b, size_t n);
    char *name;
} CanonicalTransform;

/** Symplectic matrix J = [[0, I_n], [-I_n, 0]] in R^{2n×2n}.
 *  Properties: J^2 = -I, J^T = -J, det(J) = 1.
 *  The symplectic form is omega(u, v) = u^T J v.
 */
typedef struct {
    size_t n;           /**< Half-dimension (DOF) */
    double *data;       /**< Flattened 2n × 2n matrix, row-major */
} SymplecticMatrix;

/* ──────────────────────────────────────────────────────────────
 * L2: Core concepts — algebraic structures
 * ────────────────────────────────────────────────────────────── */

/** Poisson bracket algebra closure.
 *  The set of smooth functions C^inf(R^{2n}) equipped with {·,·} forms
 *  a Lie algebra satisfying: antisymmetry, bilinearity, Leibniz rule,
 *  and Jacobi identity.  Classical limit of quantum commutator [A, B] → ih{A, B}.
 */
typedef struct {
    size_t n;
    double (*poisson)(void *ctx,
                      const double *q, const double *p,
                      double (*f)(const double *, const double *, size_t),
                      double (*g)(const double *, const double *, size_t));
    void *ctx;
} PoissonAlgebra;

/** Signature of a scalar observable on phase space.
 *  f: R^{2n} → R, e.g., energy, angular momentum, Runge-Lenz component.
 */
typedef double (*Observable)(const double *q, const double *p, size_t n);

/** Action-Angle variable system.
 *  For integrable systems (Arnold-Liouville theorem), the phase space
 *  foliates into invariant tori T^n parametrized by n action variables J_i.
 *
 *  Action: J_i = (1/2π) \oint_{gamma_i} p·dq
 *  Angle:  θ_i = ω_i t + θ_i0  (evolves uniformly)
 *  omega_i = dH0/dJ_i
 */
typedef struct {
    size_t n;
    void (*to_action_angle)(const double *q, const double *p, size_t n,
                            double *J, double *theta);
    void (*to_qp)(const double *J, const double *theta, size_t n,
                  double *q, double *p);
    double (*frequency)(const double *J, size_t n, double *omega);
    double (*H0)(const double *J, size_t n);
    char *name;
} ActionAngleSystem;

/** Energy surface descriptor.
 *  Sigma_E = { (q, p) in R^{2n} : H(q, p) = E }
 *  For H(q, p) = T(p) + V(q), the turning points satisfy V(q) = E.
 */
typedef struct {
    double E;                       /**< Energy level */
    size_t n;                       /**< DOF */
    bool is_compact;                /**< True if Sigma_E is compact (bounded motion) */
    double *q_turning_min;          /**< Lower turning points, length n */
    double *q_turning_max;          /**< Upper turning points, length n */
} EnergySurface;

/* ──────────────────────────────────────────────────────────────
 * L3: Mathematical structures — Lie theory
 * ────────────────────────────────────────────────────────────── */

/** Momentum map (moment map) J: P → g*.
 *  For a Lie group G acting symplectically on phase space P,
 *  J is an equivariant map to the dual Lie algebra.
 *  Noether's theorem: each continuous symmetry → conserved quantity via J.
 */
typedef struct {
    size_t dim_g;       /**< Dimension of the Lie algebra g */
    size_t n;           /**< Phase space DOF */
    void (*momentum)(const double *q, const double *p, size_t n,
                     double *J_value);
    void (*infinitesimal_action)(const double *xi,
                                  const double *q, const double *p,
                                  double *delta_q, double *delta_p);
    char *name;
} MomentumMap;

/** Lie algebra structure constants c_{ij}^k.
 *  [e_i, e_j] = Sum_k c_{ij}^k e_k, where {e_i} is a basis of g.
 *  For SO(3): c_{ij}^k = epsilon_{ijk} (Levi-Civita symbol).
 */
typedef struct {
    size_t dim;         /**< Dimension of Lie algebra */
    double *constants;  /**< Flattened 3D array: c[i*dim*dim + j*dim + k] */
    char *name;
} LieAlgebraStructureConstants;

/* ──────────────────────────────────────────────────────────────
 * L4: Fundamental laws — time evolution structures
 * ────────────────────────────────────────────────────────────── */

/** Liouville density on phase space.
 *  rho(q, p, t) evolves according to:
 *    drho/dt = drho/dt + {rho, H} = 0   (Liouville equation)
 *
 *  Classical analogue of von Neumann equation ih dρ/dt = [H, ρ].
 */
typedef struct {
    size_t n_grid;                  /**< Grid resolution per dimension */
    double q_min, q_max;            /**< Coordinate range */
    double p_min, p_max;            /**< Momentum range */
    double *density;                /**< rho values, n_grid × n_grid flattened */
    double dq, dp;                  /**< Grid spacings */
} LiouvilleDensity;

/** Hamiltonian flow on phase space.
 *  The flow phi_t: R^{2n} → R^{2n} is the one-parameter group of
 *  diffeomorphisms generated by X_H = J·grad_H.
 */
typedef struct {
    HamiltonianSystem *sys;         /**< Underlying Hamiltonian system */
    double t;                       /**< Current flow time */
    double *q0, *p0;                /**< Initial condition */
} HamiltonianFlow;

/* ──────────────────────────────────────────────────────────────
 * L5: Computational methods — integrator descriptors
 * ────────────────────────────────────────────────────────────── */

/** Symplectic integrator method descriptor. */
typedef enum {
    INTEGRATOR_SYMPLECTIC_EULER,    /**< 1st order */
    INTEGRATOR_STORMER_VERLET,      /**< 2nd order: leapfrog */
    INTEGRATOR_RUTH_4TH,           /**< 4th order: Ruth scheme */
    INTEGRATOR_YOSHIDA_6TH,        /**< 6th order: Yoshida composition */
    INTEGRATOR_RK4                 /**< 4th order Runge-Kutta (non-symplectic) */
} IntegratorType;

typedef struct {
    IntegratorType type;
    size_t n_stages;
    double *coeffs;
    bool is_symplectic;
    bool is_time_reversible;
    int order;
} IntegratorDescriptor;

/** Perturbation theory descriptor.
 *  For H = H0 + eps H1 where H0 is integrable.
 */
typedef struct {
    double epsilon;
    size_t n_resonances;
    int *resonance_vectors;
    double *detuning;
} PerturbationDescriptor;

/* ──────────────────────────────────────────────────────────────
 * Utility: memory management
 * ────────────────────────────────────────────────────────────── */

PhasePoint *phase_point_alloc(size_t n);
PhasePoint *phase_point_copy(const PhasePoint *src);
void phase_point_free(PhasePoint *p);

PhaseTrajectory *phase_trajectory_alloc(size_t capacity, size_t n);
void phase_trajectory_append(PhaseTrajectory *traj,
                              double t, const double *q, const double *p,
                              double H_val);
void phase_trajectory_free(PhaseTrajectory *traj);

LiouvilleDensity *liouville_density_alloc(size_t n_grid,
                                           double q_min, double q_max,
                                           double p_min, double p_max);
void liouville_density_free(LiouvilleDensity *rho);

#ifdef __cplusplus
}
#endif

#endif /* HAMILTONIAN_TYPES_H */
