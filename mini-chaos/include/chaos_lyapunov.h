/**
 * chaos_lyapunov.h — Lyapunov 指数计算
 *
 * 定义: λ = lim_{t→∞} (1/t) log(||δx(t)|| / ||δx(0)||)
 *       其中 δx 是无穷小扰动
 *
 * 物理意义:
 *   λ > 0: 对初始条件敏感依赖 → 混沌
 *   λ = 0: 周期/拟周期运动
 *   λ < 0: 渐近稳定不动点
 *
 * 算法:
 *   两粒子法 — 最直观, 仅计算最大指数
 *   Benettin 法 — 完整 Lyapunov 谱 + Kaplan-Yorke 维数
 *
 * 参考:
 *   Benettin, Galgani, Strelcyn (1980) "Lyapunov Characteristic Exponents..."
 *   Wolf, Swift, Swinney, Vastano (1985) "Determining Lyapunov Exponents..."
 *   Strogatz §9.3
 */

#ifndef CHAOS_LYAPUNOV_H
#define CHAOS_LYAPUNOV_H

#include "chaos.h"
#include "chaos_flows.h"
#include "chaos_maps.h"

/* ──────────────────────────────────────────────────────────
 * L3: 数学结构 — 切空间动力学
 *
 * 扰动 δx 的演化方程:
 *   d(δx)/dt = J(x(t)) · δx
 *
 * 其中 J(x) = ∂f/∂x 是流场 f 的 Jacobian 矩阵
 * ────────────────────────────────────────────────────────── */

/**
 * 数值计算 Jacobian 矩阵 (中心差分)
 *
 * J_{ij} = ∂f_i/∂x_j ≈ [f_i(x+εe_j) - f_i(x-εe_j)] / (2ε)
 *
 * @param J  [out] dim×dim Jacobian 矩阵, 按行存储: J[i*dim+j]
 */
void jacobian_numerical(double *J, flow_rhs_t f, const double *x,
                        double t, int dim, const double *p,
                        double eps);

/**
 * 离散映射的 Jacobian (中心差分)
 *
 * @param map 映射函数 f: R^n → R^n
 * @param J   [out] n×n Jacobian
 */
void jacobian_map(double *J, void (*map)(double*,const double*),
                  const double *x, int dim, double eps);

/* ──────────────────────────────────────────────────────────
 * L5: Lyapunov 指数算法
 * ────────────────────────────────────────────────────────── */

/**
 * 两粒子法 — 最大 Lyapunov 指数 (连续系统)
 *
 * 原理: 追踪两条邻近轨道的分离率
 *
 * 算法:
 *   1. 初始化两条轨道, 间距 d₀
 *   2. 每次积分一步后计算当前间距 d
 *   3. λ ≈ (1/t) Σ log(d/d₀)
 *   4. 每次观测后将扰动轨道重正化到距离 d₀
 *
 * @param d0  初始扰动大小 (推荐 1e-8)
 * @return 最大 Lyapunov 指数估计
 *
 * 复杂度: O(n_steps * dim)
 */
double lyapunov_two_particle(flow_rhs_t f, const double *x0,
                              int dim, const double *p,
                              double t_end, double dt, double d0);

/**
 * Wolf 算法 — 改进的最大 Lyapunov 指数
 *
 * 与两粒子法类似, 但仅在每隔 evol_steps 步观测和重正化,
 * 避免了频繁重正化带来的数值误差
 *
 * 参考: Wolf et al. (1985)
 */
double lyapunov_wolf(flow_rhs_t f, const double *x0,
                     int dim, const double *p,
                     double t_end, double dt,
                     double d0, int evol_steps);

/**
 * Benettin 全 Lyapunov 谱算法
 *
 * 原理: 同时追踪 dim 个正交扰动方向, 周期性地 Gram-Schmidt 重正化
 *
 * 算法 (Strogatz §9.3):
 *   1. 在参考轨道上积分
 *   2. 同时积分 dim 个扰动向量 (构成切空间的一组基)
 *   3. 每隔 n_renorm 步, 对扰动向量施以 QR 分解 (或 Gram-Schmidt)
 *   4. λ_i = (1/t_norm) Σ log|R_ii|
 *
 * @param spect [out] 预先分配的 LyapunovSpectrum
 *
 * 复杂度: O(n_steps * dim³)
 */
int  lyapunov_spectrum_benettin(LyapunovSpectrum *spect,
                                 flow_rhs_t f, const double *x0,
                                 int dim, const double *p,
                                 double t_end, double dt,
                                 double d0, int n_transient,
                                 int n_renorm);

/* ──────────────────────────────────────────────────────────
 * L6: 离散映射的 Lyapunov 指数
 * ────────────────────────────────────────────────────────── */

/**
 * Logistic 映射的 Lyapunov 指数
 *
 * 解析导数: f'(x) = r(1-2x)
 * λ = (1/n) Σ_{k=1}^n log|r(1-2x_k)|
 *
 * @param n_iter       迭代次数 (不含瞬态)
 * @param n_transient  瞬态消除
 *
 * 验证:
 *   r=2.0  → λ=0  (超稳定周期2, 验证: λ=0)
 *   r=3.5  → λ≈0.318
 *   r=4.0  → λ=ln(2)≈0.693 (含非典型轨道的测度零集)
 */
double logistic_lyapunov(double r, int n_iter, int n_transient);

/**
 * Logistic 映射 Lyapunov 指数参数扫描
 *
 * @param r_vals  参数值数组
 * @param n_r     参数值个数
 * @param lyaps   [out] 对应的 Lyapunov 指数
 */
void logistic_lyapunov_scan(const double *r_vals, int n_r,
                             double *lyaps, int n_iter);

/**
 * 一维映射的 Lyapunov 指数 (数值导数)
 *
 * @param f      映射函数 f(x,r) → x_{n+1}
 * @param x0     初始值
 * @param param  参数
 * @param n_iter 迭代次数
 *
 * 用数值微分估计 f'(x) = (f(x+dx)-f(x-dx))/(2dx)
 */
double map1d_lyapunov(map1d_t f, double x0, double param,
                       int n_iter, int n_transient);

/**
 * Hénon 映射的 Lyapunov 谱
 *
 * Jacobian: J(x,y) = | -2ax  1 |
 *                     |  b   0 |
 *
 * QR 迭代: 每次迭代对 J·Q 做 QR 分解
 *
 * @param lam1,lam2 [out] 两个 Lyapunov 指数
 */
void henon_lyapunov_spectrum(double *lam1, double *lam2,
                              double a, double b,
                              int n_iter, int n_transient);

/**
 * 条件 Lyapunov 指数 (同步稳定性判据)
 *
 * 原理: 驱动系统 → 响应系统
 * 响应系统沿驱动信号的横向 Lyapunov 指数 < 0 ⇒ 同步
 *
 * @param f_drive     驱动系统
 * @param f_response  响应系统
 * @param x0          驱动系统初值
 * @param y0          响应系统初值
 * @param d0          初始横向扰动
 *
 * 参考: Pecora & Carroll (1990)
 */
double conditional_lyapunov(flow_rhs_t f_drive, flow_rhs_t f_response,
                             const double *x0, const double *y0,
                             int dim_x, int dim_y,
                             const double *p, double t_end,
                             double dt, double d0);

/* ──────────────────────────────────────────────────────────
 * L2: 核心概念 — Kaplan-Yorke (Lyapunov) 维数
 *
 * D_KY = k + Σ_{i=1}^k λ_i / |λ_{k+1}|
 * 其中 k 是满足 Σ_{i=1}^k λ_i ≥ 0 的最大整数
 *
 * 参考: Frederickson, Kaplan, Yorke, Yorke (1983)
 * ────────────────────────────────────────────────────────── */

/**
 * 从 Lyapunov 谱计算 Kaplan-Yorke 维数
 */
double kaplan_yorke_dimension(const double *lyaps, int n);

/**
 * Gram-Schmidt 正交化 (就地操作)
 *
 * @param vectors  dim×dim 矩阵, 列向量存储
 * @return 各列归一化因子的对数值 (用于 Lyapunov 累加)
 */
void gram_schmidt(double *vectors, double *log_norms, int dim);

#endif /* CHAOS_LYAPUNOV_H */
