# Mini Classical Mechanics

A collection of **from-scratch, zero-dependency C implementations** of university-level classical mechanics. Each sub-module maps to MIT (and other top-tier university) courses, translating textbook formulas from Goldstein, Landau, and Marion-Thornton into executable C simulations from first principles.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|-----------|--------|-------------|
| [mini-newtonian](mini-newtonian/) | Newton's laws, particle dynamics, harmonic oscillator, central forces, non-inertial frames | MIT 8.01, Stanford PHYSICS 41 |
| [mini-lagrangian](mini-lagrangian/) | Lagrangian formulation, generalized coordinates, calculus of variations, Euler-Lagrange equations, constraints | MIT 8.09, Landau Vol.1 |
| [mini-hamiltonian](mini-hamiltonian/) | Hamiltonian mechanics, phase space, canonical transformations, Poisson brackets, integrability, action-angle variables | MIT 8.09, Goldstein Ch.8-10 |
| [mini-rigid-body](mini-rigid-body/) | Rigid body kinematics, Euler angles, inertia tensor, Euler equations, torque-free motion, heavy top | MIT 2.003, Goldstein Ch.4-5 |
| [mini-continuum](mini-continuum/) | Continuum mechanics, stress and strain tensors, constitutive relations, elasticity, fluid preliminaries, continuum beams | MIT 2.071, Landau Vol.7 |
| [mini-celestial-mechanics](mini-celestial-mechanics/) | Kepler problem, two-body orbital dynamics, perturbation theory, restricted three-body problem, Lagrange points | MIT 8.288, Murray-Dermott |
| [mini-chaos](mini-chaos/) | Nonlinear dynamics, phase portraits, bifurcations, Lyapunov exponents, Poincare sections, KAM theory | MIT 18.S197, Strogatz |
| [mini-variational-principle](mini-variational-principle/) | Principle of least action, Hamilton's principle, Noether's theorem, symmetry and conservation laws, gauge invariance | MIT 8.09, Landau Vol.1 |

## Design Philosophy

- **Zero external dependencies** — pure C (C99/C11), only `libc` and `libm`
- **Self-contained sub-modules** — each has its own `Makefile`, `include/`, `src/`, `tests/`, `examples/`, `docs/`
- **Theory-to-code mapping** — every module includes formula derivations and course-alignment notes in `docs/`
- **Physics simulation library** — ODE integrators (RK4, Verlet, symplectic), phase space visualization, orbital propagators

## Building

Each sub-module is standalone. Navigate to a directory and run:

```bash
cd mini-newtonian
make all    # build library and examples
make test   # run tests
```

Requires **GCC** and **GNU Make**.

## Project Structure

```
0. mini-classical-mechanics/
├── mini-newtonian/               # Newton's laws, particle dynamics
├── mini-lagrangian/              # Lagrangian formulation, variations
├── mini-hamiltonian/             # Hamiltonian mechanics, phase space
├── mini-rigid-body/              # Rigid body kinematics and dynamics
├── mini-continuum/               # Continuum mechanics
├── mini-celestial-mechanics/     # Orbital and celestial dynamics
├── mini-chaos/                   # Nonlinear dynamics and chaos theory
├── mini-variational-principle/   # Least action and Noether's theorem

```

## License

MIT
