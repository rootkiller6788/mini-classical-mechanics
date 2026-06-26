/**
 * newtonian.h - Newtonian Mechanics master header (L1-L6 Complete)
 *
 * Single-include header for the entire mini-newtonian C library.
 *
 * Knowledge coverage:
 *   L1 Definitions: Vec3, ParticleState, NBodySystem, Trajectory
 *   L2 Concepts: reference frames, energy, momentum, constraints
 *   L3 Math: vector calculus, coordinate transforms, Frenet frame
 *   L4 Laws: Newton's 2nd, gravitation, Hooke, drag, friction, Lorentz
 *   L5 Algorithms: Euler, RK2, RK4, RK45, Verlet, Leapfrog, Yoshida
 *   L6 Canonical: SHO, projectile, pendulum, Kepler, collisions, N-body
 *
 * Course mapping:
 *   MIT 8.012 Physics I: Classical Mechanics
 *   Goldstein, Poole, Safko: Classical Mechanics (3rd Ed.)
 *   Feynman Lectures on Physics, Vol. I
 *   Kleppner & Kolenkow: An Introduction to Mechanics
 *
 * Zero external dependencies beyond libc and libm.
 */
#ifndef NEWTONIAN_H
#define NEWTONIAN_H

#include "vec3.h"
#include "kinematics.h"
#include "forces.h"
#include "integrators.h"
#include "energy.h"
#include "momentum.h"
#include "constraints.h"
#include "analysis.h"
#include "analysis.h"

/* Physical constants (CODATA 2018 / IAU 2015) */
#ifndef G_CONST
#define G_CONST     6.67430e-11   /* gravitational constant [m^3/(kg.s^2)] */
#endif
#ifndef G_EARTH
#define G_EARTH     9.80665       /* standard gravity at Earth surface [m/s^2] */
#endif
#ifndef AU
#define AU          1.495978707e11 /* astronomical unit [m] */
#endif
#ifndef M_SUN
#define M_SUN       1.98847e30    /* solar mass [kg] */
#endif
#ifndef M_EARTH
#define M_EARTH     5.9722e24     /* Earth mass [kg] */
#endif
#ifndef R_EARTH
#define R_EARTH     6.371e6       /* Earth mean radius [m] */
#endif
#ifndef DAY_SEC
#define DAY_SEC     86400.0       /* seconds per day */
#endif
#ifndef YEAR_SEC
#define YEAR_SEC    31557600.0    /* seconds per Julian year */
#endif

#endif /* NEWTONIAN_H */
