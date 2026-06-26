# Lagrangian Mechanics Cheatsheet

## Common Lagrangian Densities

| System | Lagrangian L | Generalized Coords | Notes |
|--------|-------------|-------------------|-------|
| Free particle | ½mv² | (x,y,z) | |
| SHO | ½mẋ² - ½kx² | x | |
| Simple pendulum | ½mL²θ̇² + mgL cos θ | θ | |
| Double pendulum | ½(m₁+m₂)L₁²θ̇₁² + ½m₂L₂²θ̇₂² + m₂L₁L₂θ̇₁θ̇₂cos(θ₁-θ₂) + (m₁+m₂)gL₁cos θ₁ + m₂gL₂cos θ₂ | (θ₁,θ₂) | |
| Spherical pendulum | ½mL²(θ̇² + sin²θ φ̇²) + mgL cos θ | (θ,φ) | |
| Central force (polar) | ½m(ṙ² + r²θ̇²) - U(r) | (r,θ) | |
| Bead on rotating hoop | ½mR²(θ̇² + ω² sin²θ) - mgR(1-cos θ) | θ | Effective U_eff included |
| Charged particle in EM | ½mv² - q(φ - v·A) | (x,y,z) | Velocity-dependent potential |
| Atwood machine | ½(m₁+m₂)ẋ² + (m₁-m₂)gx | x | x = displacement of m₁ downward |
| Coupled oscillators | ½m(ẋ₁²+ẋ₂²) - ½k(x₁²+x₂²) - ½kc(x₁-x₂)² | (x₁,x₂) | kc = coupling spring |

## Noether's Theorem Correspondence

| Symmetry | Generator Q | Conserved Quantity | Condition |
|----------|------------|-------------------|-----------|
| Time translation t→t+ε | — | Energy E = p·q̇ - L | ∂L/∂t = 0 |
| Spatial translation q→q+ε·n | Q = n | p·n (linear momentum along n) | ∂L/∂q_trans = 0 |
| Rotation q→q+ε·(n×r) | Q = n×r | L·n (angular momentum along n) | ∂L/∂θ = 0 |
| Gauge transform A→A+∇χ | — | Charge conservation | Gauge invariance of L_EM |

## Euler-Lagrange Equation

```
d/dt (∂L/∂q̇) - ∂L/∂q = 0     (no constraints)
d/dt (∂L/∂q̇) - ∂L/∂q = Σ λ_k·(∂f_k/∂q)  (with constraints f_k=0)
```

## Legendre Transform: L → H

```
H(q,p,t) = p·q̇ - L(q,q̇,t)    where p = ∂L/∂q̇
```

If L = T - U and T is homogeneous quadratic in q̇: H = T + U = E (total energy)

## References

- Goldstein Ch.1-2 — Lagrangian formulation
- Goldstein Ch.13 — Noether's theorem
- MIT 8.012 — Lectures 20-23
- Landau & Lifshitz Vol.1 — Mechanics, Ch.1
