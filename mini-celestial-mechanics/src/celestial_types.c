/**
 * celestial_types.c — Core Data Type Implementations
 *
 * Vector/Matrix algebra, type constructors, validators, and numeric utilities.
 *
 * References:
 *   - Goldstein, Poole & Safko §3.5-3.7
 *   - Murray & Dermott §2.5-2.8
 *   - Vallado §2.2-2.4
 *
 * Knowledge: L1 Definitions, L3 Mathematical Structures
 */

#include "celestial_types.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ================================================================
 * Vector3 operations
 * ================================================================ */

Vector3 vec3_zero(void) {
    Vector3 v = {0.0, 0.0, 0.0};
    return v;
}

Vector3 vec3_new(double x, double y, double z) {
    Vector3 v = {x, y, z};
    return v;
}

double vec3_dot(Vector3 a, Vector3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3 vec3_cross(Vector3 a, Vector3 b) {
    Vector3 c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}

double vec3_norm(Vector3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vector3 vec3_unit(Vector3 v) {
    double n = vec3_norm(v);
    if (n < 1e-30) {
        return vec3_zero();
    }
    return vec3_scale(1.0 / n, v);
}

Vector3 vec3_add(Vector3 a, Vector3 b) {
    Vector3 c = {a.x + b.x, a.y + b.y, a.z + b.z};
    return c;
}

Vector3 vec3_sub(Vector3 a, Vector3 b) {
    Vector3 c = {a.x - b.x, a.y - b.y, a.z - b.z};
    return c;
}

Vector3 vec3_scale(double s, Vector3 v) {
    Vector3 r = {s * v.x, s * v.y, s * v.z};
    return r;
}

Vector3 vec3_lincomb(double s, Vector3 a, double t, Vector3 b) {
    Vector3 r;
    r.x = s * a.x + t * b.x;
    r.y = s * a.y + t * b.y;
    r.z = s * a.z + t * b.z;
    return r;
}

/* ================================================================
 * Matrix33 operations (rotation matrices for 3-1-3 Euler sequence)
 * ================================================================ */

Matrix33 mat33_identity(void) {
    Matrix33 M = {{{1,0,0},{0,1,0},{0,0,1}}};
    return M;
}

/** Rotation about x-axis by `angle` radians:
 *  R_x(θ) = [[1,0,0],[0,cosθ,-sinθ],[0,sinθ,cosθ]] */
Matrix33 mat33_rotx(double angle) {
    double c = cos(angle);
    double s = sin(angle);
    Matrix33 M = {{{1,0,0},{0,c,-s},{0,s,c}}};
    return M;
}

/** Rotation about z-axis by `angle` radians:
 *  R_z(θ) = [[cosθ,-sinθ,0],[sinθ,cosθ,0],[0,0,1]] */
Matrix33 mat33_rotz(double angle) {
    double c = cos(angle);
    double s = sin(angle);
    Matrix33 M = {{{c,-s,0},{s,c,0},{0,0,1}}};
    return M;
}

Vector3 mat33_mul_vec(Matrix33 M, Vector3 v) {
    Vector3 r;
    r.x = M.m[0][0]*v.x + M.m[0][1]*v.y + M.m[0][2]*v.z;
    r.y = M.m[1][0]*v.x + M.m[1][1]*v.y + M.m[1][2]*v.z;
    r.z = M.m[2][0]*v.x + M.m[2][1]*v.y + M.m[2][2]*v.z;
    return r;
}

Matrix33 mat33_mul(Matrix33 A, Matrix33 B) {
    Matrix33 C;
    int i, j, k;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            double sum = 0.0;
            for (k = 0; k < 3; k++) {
                sum += A.m[i][k] * B.m[k][j];
            }
            C.m[i][j] = sum;
        }
    }
    return C;
}

/* ================================================================
 * Orbit type classification
 * ================================================================ */

OrbitType orbit_type_from_e(double e) {
    if (e < 0.0) return ORBIT_ELLIPTIC; /* invalid but handle */
    if (e < 1e-12) return ORBIT_CIRCULAR;
    if (e < 1.0 - 1e-12) return ORBIT_ELLIPTIC;
    if (fabs(e - 1.0) < 1e-12) return ORBIT_PARABOLIC;
    return ORBIT_HYPERBOLIC;
}

OrbitType orbit_type_from_elements(const OrbitalElements *el) {
    if (!el) return ORBIT_ELLIPTIC;
    return orbit_type_from_e(el->e);
}

/* ================================================================
 * OrbitalElements validation
 *
 * Checks: e ≥ 0, a consistent with e for elliptic (a>0),
 * inclination in [0,π], angular elements physically meaningful.
 * ================================================================ */

int orbital_elements_valid(const OrbitalElements *el) {
    if (!el) return 0;
    if (el->e < 0.0) return 0;
    if (el->e < 1.0 && el->a <= 0.0) return 0;  /* bound orbit needs a>0 */
    if (el->i < 0.0 || el->i > M_PI) return 0;
    if (isfinite(el->a) && fabs(el->a) < 1e-30) return 0;
    if (!isfinite(el->e)) return 0;
    if (!isfinite(el->i)) return 0;
    return 1;
}

/* ================================================================
 * KeplerOrbit construction
 *
 * From elements and μ, compute all derived quantities:
 * p = a(1-e²), r_p = a(1-e), r_a = a(1+e) (elliptic),
 * ε = -μ/(2a), T = 2π√(a³/μ) (elliptic), h = √(μp)
 *
 * Reference: Goldstein §3.5-3.6, Murray & Dermott §2.5
 * ================================================================ */

KeplerOrbit kepler_orbit_new(const OrbitalElements *el, double mu) {
    KeplerOrbit ko;
    memset(&ko, 0, sizeof(ko));

    if (!el) return ko;

    ko.elements = *el;
    ko.mu = mu;

    double a = el->a;
    double e = el->e;

    ko.p = a * (1.0 - e * e);
    ko.r_p = a * (1.0 - e);

    if (e < 1.0) {
        ko.r_a = a * (1.0 + e);
        ko.T = TWOPI * sqrt(a * a * a / mu);
        ko.n = TWOPI / ko.T;
    } else {
        ko.r_a = INFINITY;
        ko.T = INFINITY;
        ko.n = 0.0;
    }

    ko.epsilon = -mu / (2.0 * a);
    ko.h_mag = sqrt(mu * ko.p);

    /* Angular momentum vector direction (orbit normal):
     * h points along z of orbit plane; actual direction determined
     * by inclination and RAAN. Here we store the magnitude-oriented
     * vector in the orbit-plane frame (before rotation). */
    ko.h_vec = vec3_new(0.0, 0.0, ko.h_mag);

    return ko;
}

int kepler_orbit_valid(const KeplerOrbit *ko) {
    if (!ko) return 0;
    if (!orbital_elements_valid(&ko->elements)) return 0;
    if (ko->mu <= 0.0) return 0;
    if (!isfinite(ko->p) || ko->p <= 0.0) return 0;
    if (!isfinite(ko->h_mag) || ko->h_mag <= 0.0) return 0;
    return 1;
}

/* ================================================================
 * Other constructors
 * ================================================================ */

StateVector state_vector_new(Vector3 r, Vector3 v) {
    StateVector sv;
    sv.r = r;
    sv.v = v;
    return sv;
}

TransferOrbit transfer_orbit_default(void) {
    TransferOrbit t;
    memset(&t, 0, sizeof(t));
    strcpy(t.name, "default");
    return t;
}

/* ================================================================
 * Numeric utility functions
 * ================================================================ */

double clamp_val(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

double wrap_2pi(double angle) {
    double a = fmod(angle, TWOPI);
    if (a < 0.0) a += TWOPI;
    return a;
}

double sign_val(double x) {
    if (x > 0.0) return 1.0;
    if (x < 0.0) return -1.0;
    return 0.0;
}

int feq(double a, double b, double tol) {
    return fabs(a - b) < tol;
}
