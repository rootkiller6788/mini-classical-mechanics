# Hamiltonian Mechanics Cheatsheet

## Hamilton's Equations

```
dq/dt = ∂H/∂p    (n equations)
dp/dt = -∂H/∂q   (n equations)
```

For standard systems H = T + U = p'M^{-1}p/2 + U(q):
```
dq/dt = M^{-1}p    → p = M qdot
dp/dt = -∇U + gradient terms from M(q)
```

## Poisson Bracket Properties

| Property | Formula |
|----------|---------|
| Antisymmetry | {f,g} = -{g,f} |
| Bilinearity | {af+bg, h} = a{f,h} + b{g,h} |
| Leibniz rule | {fg, h} = f{g,h} + {f,h}g |
| Jacobi identity | {f,{g,h}} + {g,{h,f}} + {h,{f,g}} = 0 |
| Fundamental | {q_i, q_j}=0, {p_i, p_j}=0, {q_i, p_j}=δ_ij |

## Canonical Transformations

| Type | Generator | Equations |
|------|-----------|-----------|
| F1(q,Q) | F1 | p=∂F1/∂q, P=-∂F1/∂Q |
| F2(q,P) | F2 | p=∂F2/∂q, Q=∂F2/∂P |
| F3(p,Q) | F3 | q=-∂F3/∂p, P=-∂F3/∂Q |
| F4(p,P) | F4 | q=-∂F4/∂p, Q=∂F4/∂P |

## Liouville's Theorem

Phase space density ρ(q,p,t) satisfies:
```
dρ/dt = ∂ρ/∂t + {ρ, H} = 0
```
→ Phase space volume is conserved under Hamiltonian flow.

## Action-Angle Variables

For 1D integrable systems:
```
J = (1/2π) ∮ p dq    (action)
θ = ω t + θ₀         (angle, evolves uniformly)
ω = ∂H/∂J = 2π/T     (frequency)
```

## Common Hamiltonian Systems

| System | H(q,p) |
|--------|--------|
| Free particle | p²/(2m) |
| SHO | p²/(2m) + ½mω²q² |
| Pendulum | p²/(2mL²) - mgL cos θ |
| Central force | p_r²/(2m) + p_θ²/(2mr²) + U(r) |
| Charged particle in B | (p - qA)²/(2m) |
| Kepler | p²/(2m) - k/r |

## Symplectic Integrators

| Method | Order | Key Property |
|--------|-------|-------------|
| Symplectic Euler | 1 | Phase area preserved exactly for linear systems |
| Stormer-Verlet | 2 | Time-reversible, symplectic |
| Ruth 4th order | 4 | Higher-order symplectic |
