# Lagrangian.jl — 拉格朗日力学主模块
# 依赖: mini-newtonian (Vec3 + integrators)
# 参考: Goldstein, Landau & Lifshitz Vol.1, MIT 8.012

module Lagrangian

using LinearAlgebra

# 导出核心类型
export GeneralizedCoords, GeneralizedVelocities, GeneralizedState
export EulerLagrangeSystem, StandardLagrangian
# 导出坐标变换
export polar_to_cartesian, polar_jacobian
export spherical_to_cartesian, spherical_jacobian
export cylindrical_to_cartesian, cylindrical_jacobian
export elliptic_to_cartesian, elliptic_jacobian
export double_pendulum_cartesian, double_pendulum_jacobian
export mass_matrix, generalized_momenta
export christoffel_symbols, geodesic_acceleration
export double_pendulum_mass_matrix, spherical_pendulum_mass_matrix
export central_force_mass_matrix, em_mass_matrix
# 导出 EL 求解
export el_to_first_order, standard_el_to_ode
# 导出作用量
export action, midpoint_discrete_lagrangian, verify_least_action
export first_variation, second_variation_sign
export symplectic_midpoint_step, discrete_action
export linearized_action_hessian
export solve_variational, variational_midpoint_step
# 导出 Noether
export noether_charge, energy_from_lagrangian
export translation_symmetry_Q, rotation_symmetry_Q_z
export verify_conservation, boost_symmetry_Q, scale_symmetry_Q
export field_noether_current, lie_algebra_closure
export verify_symmetry_invariance
# 导出约束系统
export Constraint, ConstrainedLagrangianSystem
export holonomic_constraint, nonholonomic_constraint, constrained_el_to_ode
export spherical_constraint, rigid_rod_constraint, rolling_disk_constraints
export curve_constraint
# 导出 Legendre 变换
export LegendreTransform, standard_legendre
export qdot_to_momentum, momentum_to_qdot, hamiltonian
export hamiltons_equations, poisson_bracket, routhian
export verify_canonical_poisson, routhian_equations
export is_legendre_regular, fenchel_conjugate
# 导出小振动
export SmallOscillationSystem, find_equilibrium, solve_normal_modes
export numerical_hessian, numerical_gradient
export normal_mode_solution, analyze_small_oscillations
export forced_response, frequency_sweep
export damping_ratio, quality_factor, resonant_frequency, rayleigh_damping
# 导出 ODE 积分器
export rk4_step, rk4_integrate, verlet_step, verlet_integrate
export symplectic_euler_step, velocity_verlet_step, leapfrog_step
# 导出中心力问题
export CentralForceSystem, effective_potential, orbital_equation
export kepler_orbit, solve_central_force, binet_formula
export scattering_angle, differential_cross_section
export laplace_runge_lenz, verify_kepler_laws
# 导出刚体动力学
export EulerAngles, rotation_matrix_euler, body_to_space
export inertia_tensor_point, inertia_tensor_parallel_axis
export rigid_body_kinetic_energy, euler_equations_ode
export torque_free_precession, heavy_top_lagrangian, angular_velocity_body
# 导出变分法
export GateauxDerivative, functional_derivative
export euler_lagrange_variational, hamiltons_principle_action
export maupertuis_action, jacobi_metric, beltrami_identity
# 导出相对论拉格朗日量
export relativistic_particle_L, relativistic_energy
export relativistic_momentum, proper_time_lagrangian
# 导出规范理论
export gauge_transform_A, gauge_invariant_L
export verify_gauge_invariance, chern_simons_charge
# 导出微扰论
export secular_perturbation, averaging_lagrangian
export adiabatic_invariant_1d, canonical_perturbation_series

include("generalized.jl")
include("euler_lagrange.jl")
include("action.jl")
include("constraints.jl")
include("legendre.jl")
include("small_oscillations.jl")
include("noether.jl")
include("integrators.jl")
include("central_force.jl")
include("rigid_body.jl")
include("variational_calculus.jl")
include("relativistic.jl")
include("gauge_theory.jl")
include("perturbation.jl")

end
