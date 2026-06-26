/**
 * chaos_flows.h — 连续时间混沌系统与数值积分
 *
 * 系统:
 *   Lorenz (1963)  — 对流模型, 第一个数值发现的混沌吸引子
 *   Rössler (1976) — 最简单的连续混沌系统之一
 *   Chua (1983)    — 第一个物理实现的混沌电路
 *   Duffing (1918) — 非线性强迫振荡
 *   Van der Pol    — 自激振荡
 *   Sprott (1994)  — 极简混沌流
 *   Chen (1999)    — 对偶 Lorenz 系统
 *   Hénon-Heiles   — Hamiltonian 混沌
 *
 * 积分器: RK4 (经典 Runge-Kutta 4阶)
 *
 * 参考: Strogatz Ch.9, Sprott "Chaos and Time-Series Analysis" (2003)
 *       Lorenz "Deterministic Nonperiodic Flow" J.Atmos.Sci. (1963)
 */

#ifndef CHAOS_FLOWS_H
#define CHAOS_FLOWS_H

#include "chaos.h"

/* ──────────────────────────────────────────────────────────
 * L6: 经典混沌系统 — 右侧函数 (RHS)
 *
 * 每个函数实现: dx/dt 的计算
 * 签名: void rhs(double *dx, const double *x, double t, const double *params)
 * ────────────────────────────────────────────────────────── */

/**
 * Lorenz 系统 (Lorenz, 1963)
 *
 * 方程: ẋ = σ(y - x)
 *       ẏ = x(ρ - z) - y
 *       ż = xy - βz
 *
 * 经典参数: σ=10, ρ=28, β=8/3
 * 知识点: 奇怪吸引子, 蝴蝶效应, 对流不稳定性
 *
 * MIT 8.012: 相空间分析
 * Caltech Ph 106: KAM 理论前置
 */
void lorenz_rhs(double dx[3], const double x[3], double t, const double p[3]);

/**
 * Rössler 系统 (Rössler, 1976)
 *
 * 方程: ẋ = -y - z
 *       ẏ = x + ay
 *       ż = b + z(x - c)
 *
 * 经典参数: a=0.2, b=0.2, c=5.7
 * 知识点: 螺旋混沌, 漏斗吸引子
 */
void rossler_rhs(double dx[3], const double x[3], double t, const double p[3]);

/**
 * Chua 电路 (Chua, 1983)
 *
 * 方程: ẋ = α(y - x - f(x))
 *       ẏ = x - y + z
 *       ż = -βy
 *
 * f(x) = m₁x + ½(m₀ - m₁)(|x+1| - |x-1|)
 * 知识点: 分段线性混沌, 物理可实现的混沌电路, 双涡卷吸引子
 */
void chua_rhs(double dx[3], const double x[3], double t, const double p[5]);

/**
 * Duffing 振子 (Duffing, 1918)
 *
 * 方程: ẋ = v
 *       v̇ = -δv + βx - αx³ + γ cos(ωt)
 *       ṫ = 1  (自治化)
 *
 * 知识点: 非线性共振, 弹簧硬化/软化, 混沌鞍
 * Goldstein §11.5: 受驱阻尼非线性振子
 */
void duffing_rhs(double dx[3], const double x[3], double t, const double p[5]);

/**
 * 受迫 Van der Pol 振子
 *
 * 方程: ẋ = v
 *       v̇ = μ(1 - x²)v - x + F cos(ωt)
 *       ṫ = 1
 *
 * 知识点: 自激振荡, 极限环, 锁频, 准周期通往混沌
 * Strogatz §8.5: Van der Pol 方程
 */
void forced_vdp_rhs(double dx[3], const double x[3], double t, const double p[4]);

/**
 * Sprott B 系统: 仅5项, 二次非线性
 *
 * ẋ = yz, ẏ = x - y, ż = 1 - xy
 * 知识点: 最简混沌流, 代数简单性 ≠ 动力学简单性
 */
void sprott_b_rhs(double dx[3], const double x[3], double t, const double p[1]);

/**
 * Sprott C 系统
 *
 * ẋ = yz, ẏ = x - y, ż = 1 - x²
 */
void sprott_c_rhs(double dx[3], const double x[3], double t, const double p[1]);

/**
 * Chen 系统 (Chen & Ueta, 1999)
 *
 * ẋ = a(y - x)
 * ẏ = (c - a)x - xz + cy
 * ż = xy - bz
 *
 * 知识点: 对偶 Lorenz, 反混沌
 */
void chen_rhs(double dx[3], const double x[3], double t, const double p[3]);

/**
 * Hénon-Heiles 系统 (1964) — Hamiltonian 混沌
 *
 * H = ½(p_x² + p_y²) + ½(x² + y²) + x²y - ⅓y³
 *
 * ẋ = p_x, ẏ = p_y
 * ṗ_x = -x - 2xy, ṗ_y = -y - x² + y²
 *
 * 知识点: 可积 vs 非可积, 庞加莱截面, KAM 定理数值验证
 * Goldstein §11.3: KAM 定理
 */
void henon_heiles_rhs(double dx[4], const double x[4], double t, const double p[1]);

/**
 * 4D 超混沌 Rössler 系统
 *
 * ẋ = -y - z
 * ẏ = x + 0.25y + w
 * ż = 3 + xz
 * ẇ = -0.5z + 0.05w
 *
 * 知识点: 超混沌 (多个正 Lyapunov 指数)
 */
void hyperchaos_rossler_rhs(double dx[4], const double x[4], double t, const double p[1]);

/* ──────────────────────────────────────────────────────────
 * L5: 数值积分算法
 * ────────────────────────────────────────────────────────── */

/**
 * 流系统右侧函数签名
 *
 * @param dx    [out] 导数
 * @param x     [in]  当前状态
 * @param t     [in]  当前时间
 * @param p     [in]  参数数组
 */
typedef void (*flow_rhs_t)(double *dx, const double *x, double t, const double *p);

/**
 * RK4 单步积分
 *
 * 算法: k₁ = h f(xₙ)
 *       k₂ = h f(xₙ + k₁/2)
 *       k₃ = h f(xₙ + k₂/2)
 *       k₄ = h f(xₙ + k₃)
 *       xₙ₊₁ = xₙ + (k₁ + 2k₂ + 2k₃ + k₄)/6
 *
 * 局部截断误差: O(h⁵)
 * 全局截断误差: O(h⁴)
 *
 * 参考: Butcher (2008) "Numerical Methods for Ordinary Differential Equations"
 */
void rk4_step(double *x_next, flow_rhs_t f, const double *x,
              double t, double h, int dim, const double *p);

/**
 * 自适应步长 RK4 (基于 Richardson 外推)
 *
 * 每步用 h 和 h/2 两步计算, 估计局部误差,
 * 按 tol 调整步长
 *
 * 返回: 实际使用的步数
 */
int  rk4_adaptive(double *traj_out, flow_rhs_t f, const double *x0,
                  int dim, const double *p, double t0, double t_end,
                  double dt_init, double tol, int max_steps);

/**
 * 完整轨迹积分 (定步长)
 *
 * 返回轨迹: 包含 (t_end/dt + 1) 个点
 *
 * @return ChaosTrajectory (调用者负责 free)
 */
ChaosTrajectory* integrate_flow(flow_rhs_t f, const double *x0,
                                 int dim, const double *p,
                                 double t_end, double dt);

/**
 * Poincaré 截面提取
 *
 * 截面条件: x[section_dim] 穿过 section_value 且符号变化
 *
 * @param section_dim  截面变量索引
 * @param section_value 截面值 (通常=0)
 *
 * 参考: Strogatz §8.7, Goldstein §11.2
 */
PoincareSection* poincare_section_extract(flow_rhs_t f, const double *x0,
                                           int dim, const double *p,
                                           double t_end, double dt,
                                           int section_dim,
                                           double section_value);

/**
 * 绘制双参数分岔图的轨迹采样
 *
 * 返回: 瞬态消除后的渐近解采样
 */
double* sample_asymptotic(flow_rhs_t f, const double *x0, int dim,
                           const double *p, double t_transient,
                           double t_sample, double dt,
                           int n_samples, int var_idx);

#endif /* CHAOS_FLOWS_H */
