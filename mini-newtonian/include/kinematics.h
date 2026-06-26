/**
 * kinematics.h — 经典运动学 (L1 Definitions / L4 Fundamental Laws)
 *
 * 位置、速度、加速度之间的关系。不含力——仅描述运动。
 *
 * 课程对标:
 *   MIT 8.012 §1.4-2.4 Kinematics in 1D/2D/3D
 *   Goldstein §1.1-1.2 Kinematics of a particle
 *   Feynman Lectures Vol.I §8 Motion
 *
 * 核心定理:
 *   1. 牛顿第一定律 (惯性定律): ΣF=0 ⇒ v=const
 *   2. 伽利略相对性原理: 物理定律在所有惯性系中相同
 *   3. 匀速圆周运动: a = v²/r (向心加速度)
 */

#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "vec3.h"

/* ================================================================
 * L1: 一维匀加速运动方程 (SUVAT)
 * 定理来源: Galileo Galilei, Discorsi (1638) — 匀加速运动
 * ================================================================ */

double suvat_position(double r0, double v0, double a, double t);
double suvat_velocity(double v0, double a, double t);
double suvat_velocity_from_displacement(double v0, double a, double dr);

/* ================================================================
 * L1: 三维匀加速运动
 * ================================================================ */

Vec3 uniform_motion(Vec3 r0, Vec3 v, double t);          /* r(t) = r₀ + vt */
Vec3 accelerated_motion(Vec3 r0, Vec3 v0, Vec3 a, double t); /* r(t) = r₀+v₀t+½at² */
Vec3 accelerated_velocity(Vec3 v0, Vec3 a, double t);   /* v(t) = v₀ + at */

/* ================================================================
 * L2: 抛体运动 (均匀重力场中的运动)
 * 知识点: 抛物线轨迹 — 水平匀速 + 竖直匀加速 的独立叠加
 * 定理来源: Galileo's projectile theorem (1638)
 * 课程: MIT 8.012 §2.4
 * ================================================================ */

typedef struct {
    Vec3 position;
    Vec3 velocity;
} ProjectileState;

ProjectileState projectile_at_time(Vec3 r0, Vec3 v0, double t, Vec3 gravity);
double projectile_flight_time(double v0z, double z0, double g);
double projectile_range(double v0, double theta, double g);
double projectile_max_height(double v0, double theta, double g);
double projectile_range_from_height(double v0, double theta, double h0, double g);
double optimal_launch_angle(double v0, double target_range, double g);
double projectile_trajectory_height(double x, double v0, double theta, double g);

/* ================================================================
 * L3: 圆周运动
 * 知识点: 匀速圆周运动 → 向心加速度 a = ω²r = v²/r
 * ================================================================ */

Vec3 circular_position(Vec3 center, Vec3 axis_normal,
                       double radius, double omega, double t);
Vec3 circular_velocity(Vec3 axis_normal, double radius, double omega, double t);
Vec3 circular_acceleration(Vec3 center, Vec3 position,
                           double omega); /* 向心 + 切向 */
double centripetal_acceleration(double v, double r);
double centripetal_acceleration_omega(double omega, double r);
double circular_period(double omega);                    /* T = 2π/ω */
double circular_frequency(double period);                /* ω = 2π/T */

/* ================================================================
 * L3: 非均匀圆周运动
 * 知识点: 切向加速度 + 法向加速度的分解
 * a = α×r + ω×(ω×r) = a_tangential + a_centripetal
 * ================================================================ */

typedef struct {
    double centripetal_mag;   /* 向心加速度大小 v²/r */
    double tangential_mag;    /* 切向加速度大小 r*α */
    Vec3   centripetal_dir;   /* 向心方向 (指向圆心) */
    Vec3   tangential_dir;    /* 切向方向 (速度方向) */
} CurvilinearAccel;

CurvilinearAccel decompose_acceleration(Vec3 position, Vec3 velocity,
                                        Vec3 acceleration, Vec3 center);

/* ================================================================
 * L2: 伽利略变换 (惯性参考系之间的变换)
 * 知识点: Galileo 相对性原理
 * ================================================================ */

Vec3 galilean_position_transform(Vec3 r_prime, Vec3 V_frame, double t);
Vec3 galilean_velocity_transform(Vec3 v_prime, Vec3 V_frame);
Vec3 galilean_acceleration_invariance(Vec3 a_prime); /* a = a' (不变) */

/* ================================================================
 * L4: 旋转参考系 (非惯性系)
 * 知识点: 科里奥利力 + 离心力 + Euler 力
 * 定理: a_inertial = a_rot + 2ω×v_rot + ω×(ω×r) + dω/dt×r
 * 课程: MIT 8.012 §7.1-7.3, Goldstein §4.9-4.10
 * ================================================================ */

Vec3 coriolis_acceleration(Vec3 v_rot, Vec3 omega);     /* -2ω×v  */
Vec3 centrifugal_acceleration(Vec3 r_rot, Vec3 omega);  /* -ω×(ω×r) */
Vec3 euler_acceleration(Vec3 r, Vec3 omega_dot);        /* -ω̇×r */
Vec3 rotating_frame_accel(Vec3 a_inertial, Vec3 v_rot,
                          Vec3 r_rot, Vec3 omega);
Vec3 inertial_to_rotating_velocity(Vec3 v_inertial, Vec3 r, Vec3 omega);

/* ================================================================
 * L5: 轨迹曲率与 Frenet-Serret 标架
 * 知识点: 空间曲线微分几何 — 曲率/挠率/Frenet 标架
 * 课程: MIT 8.012 supplemental, Calculus on Manifolds
 * ================================================================ */

double trajectory_curvature(Vec3 velocity, Vec3 acceleration);
double trajectory_radius_of_curvature(Vec3 velocity, Vec3 acceleration);
Vec3 frenet_tangent(Vec3 velocity);        /* T = v/|v| (单位切向量) */
Vec3 frenet_normal(Vec3 velocity, Vec3 acceleration);  /* N = dT/ds 方向 */
Vec3 frenet_binormal(Vec3 velocity, Vec3 acceleration);/* B = T×N */

#endif /* KINEMATICS_H */
