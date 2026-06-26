# Hamiltonian.jl — 哈密顿力学主模块

module Hamiltonian

using LinearAlgebra

export PhasePoint, PhaseTrajectory, record!, phase_volume_element, poincare_section
export HamiltonianSystem, analytic_hamiltonian_system, numerical_hamiltonian_system
export hamiltons_ode, symplectic_euler_step, stormer_verlet_step, solve_hamiltonian
export poisson_bracket, fundamental_poisson_brackets, is_constant_of_motion
export verify_jacobi_identity, time_derivative_via_poisson
export CanonicalTransform, identity_transform, point_transform, scale_transform
export exchange_transform, verify_canonical
export generate_ensemble, verify_liouville
export HarmonicOscillatorActionAngle, action_from_energy_ho, energy_from_action_ho
export action_to_qp_ho, qp_to_action_angle_ho
export PendulumActionAngle, pendulum_energy, numerical_action, numerical_frequency
# H-J 方程
export hamilton_jacobi_1d, solve_p_from_H, harmonic_oscillator_HJ, is_separable
# 可积系统
export KeplerProblem, kepler_energy, kepler_angular_momentum
export runge_lenz_vector, kepler_orbit_params, verify_kepler_integrals
export are_in_involution

include("phase_space.jl")
include("hamiltons_equations.jl")
include("poisson.jl")
include("canonical_transform.jl")
include("liouville.jl")
include("action_angle.jl")
export PerturbedSystem, average_perturbation, corrected_frequencies
export DuffingOscillator, duffing_H0, duffing_H1, duffing_frequency
export DrivenOscillator, resonance_detuning, pendulum_approximation_width

include("hamilton_jacobi.jl")
include("integrable_systems.jl")
include("perturbation.jl")
include("hamiltonian_flows.jl")

end
