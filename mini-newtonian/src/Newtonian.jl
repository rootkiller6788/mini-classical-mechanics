# Newtonian.jl — 牛顿力学主模块
# 零外部依赖，全部从零手写
#
# 用法:
#   include("src/Newtonian.jl")
#   using .Newtonian
#
# 参考: MIT 8.012, Goldstein Classical Mechanics

module Newtonian

# 导出核心类型
export Vec3, State, NBodyState, Trajectory, ODEProblem, ParticleODE, ODESolution
# 导出常量
export G, G_AU, G_EARTH, DAY, YEAR, AU, M_SUN, M_EARTH, R_EARTH
# 导出向量运算
export dot, cross, norm, norm2, normalize, distance, angle
# 导出运动学
export projectile_motion, projectile_range, projectile_flight_time
export circular_position, circular_velocity, centripetal_acceleration
# 导出 ODE 积分器
export solve_fixed_step, solve_adaptive
export euler_step, rk2_step, rk4_step, rk45_step
export euler_cromer_step, velocity_verlet_step, rk4_second_order_step
export leapfrog_init, leapfrog_step
export nbody_verlet_step
# 导出力学
export newton_gravity, gravity_nbody, gravity_acceleration, uniform_gravity
export hooke_force, kinetic_friction, linear_drag, quadratic_drag
export lorentz_force
# 导出能量
export kinetic_energy, gravitational_potential, elastic_potential
export total_gravitational_potential, total_kinetic_energy
export energy_drift
# 导出动量
export linear_momentum, angular_momentum, center_of_mass
export elastic_collision_1d, elastic_collision_3d, inelastic_collision_1d
export tsiolkovsky_delta_v
# 导出约束系统
export incline_acceleration, angle_of_repose, normal_force_incline
export atwood_machine, pulley_system_mechanical_advantage
export conical_pendulum_parameters, conical_pendulum_critical_omega
export loop_the_loop_min_speed, surface_normal_force, detachment_condition
export incline_pulley_system, incline_pulley_direction
export effective_gravity, apparent_weight_in_elevator
export banked_curve_angle, max_safe_speed_banked_curve

# ========== 类型系统 ==========
include("types.jl")

# ========== 运动学 ==========
include("kinematics.jl")

# ========== 力定律 ==========
include("forces.jl")

# ========== ODE 积分器 ==========
include("integrators.jl")

# ========== 能量 ==========
include("energy.jl")

# ========== 动量与碰撞 ==========
include("momentum.jl")

# ========== 约束系统（滑轮/斜面/张力） ==========
include("constraints.jl")

end # module Newtonian
