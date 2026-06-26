/**
 * celestial_types.h — Core Data Types for Celestial Mechanics
 *
 * Defines the fundamental data structures: orbital elements, state vectors,
 * perturbation accelerations, Lagrange points, and transfer orbits.
 *
 * References:
 *   - Goldstein, Poole & Safko, Classical Mechanics (3rd Ed.), Ch.3
 *   - Murray & Dermott, Solar System Dynamics (1999), Ch.2
 *   - Vallado, Fundamentals of Astrodynamics and Applications (4th Ed.)
 *
 * Knowledge Coverage:
 *   L1 Definitions — All core struct/typedef definitions
 *   L2 Core Concepts — Orbital elements, state vectors, perturbation models
 *   L3 Math Structures — Vector3, Matrix33, rotation matrices
 */

#ifndef CELESTIAL_TYPES_H
#define CELESTIAL_TYPES_H

#include <math.h>

/* Fundamental Constants */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TWOPI (2.0 * M_PI)
#define HALFPI (M_PI / 2.0)

/* Gravitational parameters (km^3/s^2) */
#define MU_SUN      1.32712440018e11
#define MU_EARTH    398600.4418
#define MU_MOON     4902.8001
#define MU_MARS     42828.3
#define MU_JUPITER  1.26686534e8

/* Gaussian constant (AU^3/yr^2/M_sun) */
#define G_SUN_AU (4.0 * M_PI * M_PI)

/* Earth properties */
#define R_EARTH_KM    6378.1363
#define J2_EARTH      0.0010826265
#define EARTH_FLATTENING 0.00335281066

/* Astronomical constants */
#define AU_KM          149597870.7
#define C_LIGHT_KM_S   299792.458
#define P_SUN_KM       (4.56e-15)  /* SRP at 1 AU, N/km^2 */
#define SECONDS_PER_DAY  86400.0
#define SECONDS_PER_YEAR 31557600.0

/* Orbit types */
typedef enum {
    ORBIT_CIRCULAR    = 0,
    ORBIT_ELLIPTIC    = 1,
    ORBIT_PARABOLIC   = 2,
    ORBIT_HYPERBOLIC  = 3,
    ORBIT_RECTILINEAR = 4
} OrbitType;

/* Perturbation flags (bitmask for combined perturbation models) */
typedef enum {
    PERTURB_NONE       = 0,
    PERTURB_J2         = 1 << 0,
    PERTURB_DRAG       = 1 << 1,
    PERTURB_SRP        = 1 << 2,
    PERTURB_THIRD_BODY = 1 << 3,
    PERTURB_GR         = 1 << 4,
    PERTURB_ALL        = 0x1F
} PerturbationFlags;

/* 3D Cartesian vector */
typedef struct {
    double x, y, z;
} Vector3;

/* 3x3 matrix (row-major) */
typedef struct {
    double m[3][3];
} Matrix33;

/* ---- L1: Classical Keplerian Orbital Elements ---- */
/**
 * Six classical orbital elements (Keplerian elements).
 *   a      — semi-major axis (a>0 elliptic; a<0 hyperbolic)
 *   e      — eccentricity (0 ≤ e; e<1 bound; e≥1 unbound)
 *   i      — inclination [0, π] rad
 *   Omega  — right ascension of ascending node [0, 2π) rad
 *   omega  — argument of periapsis [0, 2π) rad
 *   nu     — true anomaly [0, 2π) rad
 *
 * Ref: Murray & Dermott Eq. (2.122), Goldstein §3.5
 */
typedef struct {
    double a;
    double e;
    double i;
    double Omega;
    double omega;
    double nu;
} OrbitalElements;

/* ---- L1: Cartesian State Vector ---- */
typedef struct {
    Vector3 r;
    Vector3 v;
} StateVector;

/* ---- L2: Complete Kepler Orbit (elements + derived quantities) ---- */
typedef struct {
    OrbitalElements elements;
    double mu;
    double T;        /* orbital period (INFINITY if e>=1) */
    double n;        /* mean motion 2π/T */
    double p;        /* semi-latus rectum a(1-e^2) */
    double r_p;      /* periapsis distance */
    double r_a;      /* apoapsis distance (INFINITY if e>=1) */
    double epsilon;  /* specific orbital energy -μ/(2a) */
    Vector3 h_vec;   /* specific angular momentum vector */
    double h_mag;    /* magnitude of h */
} KeplerOrbit;

/* ---- L2: Perturbation acceleration in RSW frame ---- */
typedef struct {
    double R;  /* radial */
    double S;  /* along-track (transverse) */
    double W;  /* cross-track (normal, along h) */
} PerturbationRSW;

/* ---- L2: J2 secular rates ---- */
typedef struct {
    double Omega_dot; /* nodal precession rate [rad/s] */
    double omega_dot; /* apsidal precession rate [rad/s] */
    double M_dot;     /* mean anomaly rate [rad/s] */
} J2SecularRates;

/* ---- L2: Lagrange points ---- */
typedef struct {
    Vector3 L1, L2, L3, L4, L5;
} LagrangePoints;

/* ---- L2: Transfer Orbit Specification ---- */
typedef struct {
    char name[64];
    OrbitalElements depart;
    OrbitalElements arrive;
    OrbitalElements transfer;
    double delta_v1;
    double delta_v2;
    double delta_v3;
    double delta_v_total;
    double transfer_time;
    double phase_angle;
} TransferOrbit;

/* ---- L2: Gravity Assist Result ---- */
typedef struct {
    double v_inf_in;
    double v_inf_out;
    double turn_angle;
    double delta_v_helio;
    double r_periapsis;
    double e_hyperbolic;
} GravityAssist;

/* ---- L3: Lambert Solution ---- */
typedef struct {
    Vector3 v1;           /* velocity at departure */
    Vector3 v2;           /* velocity at arrival */
    double a;
    double e;
    double delta_v_total;
    int    num_revs;
    int    converged;
} LambertSolution;

/* ---- L5: ODE system state (for numerical integration) ---- */
typedef struct {
    double t;       /* current time */
    double y[6];    /* state vector [x,y,z,vx,vy,vz] */
    int    dim;     /* dimension (6 for full, 4 for planar) */
} ODEState;

/* ---- Vector3 constructors and operations ---- */
Vector3 vec3_zero(void);
Vector3 vec3_new(double x, double y, double z);
double  vec3_dot(Vector3 a, Vector3 b);
Vector3 vec3_cross(Vector3 a, Vector3 b);
double  vec3_norm(Vector3 v);
Vector3 vec3_unit(Vector3 v);
Vector3 vec3_add(Vector3 a, Vector3 b);
Vector3 vec3_sub(Vector3 a, Vector3 b);
Vector3 vec3_scale(double s, Vector3 v);
Vector3 vec3_lincomb(double s, Vector3 a, double t, Vector3 b);

/* ---- Matrix operations ---- */
Matrix33 mat33_identity(void);
Matrix33 mat33_rotx(double angle);
Matrix33 mat33_rotz(double angle);
Vector3   mat33_mul_vec(Matrix33 M, Vector3 v);
Matrix33  mat33_mul(Matrix33 A, Matrix33 B);

/* ---- Constructors and validators ---- */
OrbitType orbit_type_from_e(double e);
OrbitType orbit_type_from_elements(const OrbitalElements *el);
int       orbital_elements_valid(const OrbitalElements *el);
KeplerOrbit kepler_orbit_new(const OrbitalElements *el, double mu);
int         kepler_orbit_valid(const KeplerOrbit *ko);
StateVector state_vector_new(Vector3 r, Vector3 v);
TransferOrbit transfer_orbit_default(void);

/* ---- Numeric helpers ---- */
double clamp_val(double x, double lo, double hi);
double wrap_2pi(double angle);
double sign_val(double x);
int    feq(double a, double b, double tol);

#endif /* CELESTIAL_TYPES_H */
