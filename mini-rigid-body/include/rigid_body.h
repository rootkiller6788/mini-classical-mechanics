/**
 * rigid_body.h — Mini Rigid Body Mechanics Library
 *
 * Umbrella header for the full rigid body mechanics module.
 * 
 * Coverage (Goldstein Ch.4-5, Landau Vol.1 Ch.6, MIT 8.012):
 *   L1: InertiaTensor, EulerAngles, RigidBodyState, PrincipalAxes, Quaternion
 *   L2: Torque-free Euler equations, Poinsot construction, gyroscopic effect
 *   L3: SO(3) Lie group, so(3) Lie algebra, SU(2) double cover
 *   L4: Euler's equations (rotational Newton's 2nd law), conservation laws
 *   L5: RK4, Verlet, Dormand-Prince integrators, Jacobi eigen-algorithm
 *   L6: Free rigid body, heavy symmetric top, tennis racket theorem, gyroscope
 *   L7: Satellite attitude dynamics, gyroscope navigation
 *   L8: Poinsot geometric construction, energy-momentum bifurcation
 *   L9: (documented) Quantum rigid rotor, relativistic Thomas precession
 *
 * Usage:
 *   #include "rigid_body.h"
 *   // Link with -lm
 *
 * All functions are thread-safe and reentrant (no global state).
 */

#ifndef RIGID_BODY_H
#define RIGID_BODY_H

#include "rigid_body_types.h"
#include "rigid_body_inertia.h"
#include "rigid_body_euler.h"
#include "rigid_body_kinematics.h"
#include "rigid_body_energy.h"
#include "rigid_body_tops.h"

/** Library version string */
#define RIGID_BODY_VERSION "1.0.0"

/** Mathematical constants used throughout */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_2PI
#define M_2PI 6.28318530717958647692
#endif

#endif /* RIGID_BODY_H */
