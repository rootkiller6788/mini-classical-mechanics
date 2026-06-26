/**
 * chaos_bifurcation.h — 分岔理论与分析
 *
 * 核心概念:
 *   分岔 (Bifurcation): 动力学系统定性行为随参数变化的突变
 *
 * 局部余维1分岔:
 *   鞍结点    (saddle-node):     ẋ = r + x²
 *   跨临界    (transcritical):   ẋ = rx - x²
 *   叉形      (pitchfork):       ẋ = rx - x³  (超临界) / ẋ = rx + x³ (亚临界)
 *   Hopf      (Hopf):           极坐标: ṙ = r(μ - r²), θ̇ = ω
 *
 * 全局分岔:
 *   周期倍化 (period-doubling):  不动点→周期2→周期4→混沌
 *   同宿     (homoclinic):        稳定与不稳定流形相切
 *
 * 参考:
 *   Strogatz Ch.3, Ch.8
 *   Kuznetsov "Elements of Applied Bifurcation Theory" (2004) §2.1–2.4
 *   Wiggins "Introduction to Applied Nonlinear Dynamical Systems" (2003)
 */

#ifndef CHAOS_BIFURCATION_H
#define CHAOS_BIFURCATION_H

#include "chaos.h"
#include "chaos_maps.h"

/* ──────────────────────────────────────────────────────────
 * L2: 分岔检测 — 不动点与稳定性
 * ────────────────────────────────────────────────────────── */

/**
 * 1D 映射的不动点 Newton 求解
 *
 * 求 x 使得 f(x, p) = x
 * g(x) = f(x,p) - x = 0
 *
 * 迭代: x_{n+1} = x_n - g(x_n)/g'(x_n)
 * 其中 g'(x) = f'(x,p) - 1 (数值差分)
 *
 * @param f    映射函数 f(x,p) → x_{n+1}
 * @param p    参数值
 * @param x0   初始猜测
 * @param tol  收敛容忍度
 * @param max_iter 最大迭代次数
 * @param converged [out] 是否收敛
 *
 * 参考: Strogatz §3.1, Burden & Faires §2.3
 */
double find_fixed_point(map1d_t f, double p, double x0,
                        double tol, int max_iter, int *converged);

/**
 * 不动点稳定性判定 (1D)
 *
 * 乘子: μ = |f'(x*)| 
 *  |μ| < 1: 稳定
 *  |μ| > 1: 不稳定
 *  |μ| = 1: 分岔点
 *
 * @return 0=稳定, 1=不稳定, 2=分岔(临界)
 */
int fixed_point_stability(map1d_t f, double x_star, double p, double dx);

/**
 * 平衡点稳定性 (nD 流系统)
 *
 * @param eigenvalues  Jacobian 的特征值 (dim 个)
 * @param dim          维度
 * @return FixedPointType 枚举值
 */
FixedPointType classify_equilibrium(const double *eigenvalues, int dim);

/* ──────────────────────────────────────────────────────────
 * L5: 分岔图生成
 * ────────────────────────────────────────────────────────── */

/**
 * Logistic 映射的分岔图数据
 *
 * 对每个 r 值:
 *   1. 迭代 n_transient 次丢弃瞬态
 *   2. 继续迭代 n_plot 次记录稳态值
 *
 * @param r_min, r_max  参数范围
 * @param n_r           参数采样个数
 * @param n_transient   瞬态消除
 * @param n_plot        每个参数值记录的稳态个数
 * @return BifurcationDiagram (调用者 free)
 *
 * 知识点: 周期倍化级联, 周期窗口 (如周期3窗口在 r≈3.828)
 *
 * 参考: Strogatz §10.4
 */
BifurcationDiagram* logistic_bifurcation_diagram(double r_min, double r_max,
                                                  int n_r, int n_transient,
                                                  int n_plot);

/**
 * 通用 1D 映射的分岔图
 *
 * @param f            映射函数
 * @param param_range  参数值数组
 * @param n_params     参数个数
 * @param x0           初始值
 */
BifurcationDiagram* bifurcation_diagram_1d(map1d_t f,
                                            const double *param_range,
                                            int n_params,
                                            double x0,
                                            int n_transient, int n_plot);

/**
 * 2D 映射的分岔图 (投影到指定变量)
 *
 * @param var_idx 要记录的变量索引 (0 or 1)
 */
BifurcationDiagram* bifurcation_diagram_2d(map2d_t f,
                                            const double *param_range,
                                            int n_params,
                                            double x0, double y0,
                                            double p1_fixed,
                                            int var_idx,
                                            int n_transient, int n_plot);

/* ──────────────────────────────────────────────────────────
 * L2: 周期检测
 * ────────────────────────────────────────────────────────── */

/**
 * 检测离散映射的周期
 *
 * 方法: 寻找重复值, 周期 = 两次出现之间的步数
 *
 * @param tol     重复判定容忍度
 * @param period  [out] 检测到的周期 (0=未检测到, 1=不动点)
 * @return 0=成功, -1=未检测到
 *
 * 知识点: 周期轨道定义 — x_{n+p} = x_n
 */
int detect_period_discrete(map1d_t f, double x0, double param,
                            int max_iter, double tol, int *period);

/**
 * 在给定参数处查找周期 (自动去除瞬态)
 *
 * @return 周期数, 0=未找到
 */
int find_period_at_param(map1d_t f, double param, double x0,
                          int n_transient, int n_iter, double tol);

/* ──────────────────────────────────────────────────────────
 * L4: Feigenbaum 常数与普适性
 * ────────────────────────────────────────────────────────── */

/**
 * 寻找周期倍增分岔点序列
 *
 * 从 r_start 开始, 搜索 r₂, r₄, r₈, ...
 * 方法: 递增 r, 检测当前周期
 *
 * @param r_vals     [out] 分岔点数组 (由调用者分配, 长度 n_periods)
 * @param n_periods  要找的分岔点个数
 * @param tol        周期检测容忍度
 * @return 实际找到的分岔点个数
 *
 * 参考: Strogatz §10.7
 */
int find_period_doubling_points(double *r_vals, int n_periods,
                                 double r_start, double dr_init,
                                 double tol);

/**
 * 从分岔点序列估计 Feigenbaum δ
 *
 * δ_n = (r_n - r_{n-1}) / (r_{n+1} - r_n)
 * δ = lim_{n→∞} δ_n ≈ 4.669202
 *
 * @param deltas [out] 估计的 δ 序列 (长度 n-2)
 * @return 最后一个 (最佳) 估计值
 *
 * 参考: Feigenbaum (1978)
 */
double estimate_feigenbaum_delta(const double *r_vals, int n,
                                  double *deltas);

/**
 * 分类分岔类型 (1D 映射)
 *
 * 用导数条件:
 *   saddle-node:   f(x*,r*)=x*,  ∂f/∂x=1,  ∂²f/∂x²≠0,  ∂f/∂r≠0
 *   transcritical: f(x*,r*)=x*,  ∂f/∂x=1,  ∂²f/∂x²≠0,  f(0,r)=0 ∀r
 *   pitchfork:     f(x*,r*)=x*,  ∂f/∂x=1,  ∂²f/∂x²=0,  ∂³f/∂x³≠0
 *   period-doubling: f²(x*,r*)=x*, ∂(f²)/∂x = -1
 *
 * @param x_star  不动点
 * @param r_star  分岔参数
 * @param dx      数值差分步长
 * @return BifurcationType
 */
BifurcationType classify_bifurcation(map1d_t f, double x_star,
                                      double r_star, double dx);

/**
 * 一维映射的 Sharkovskii 序查找
 *
 * 给定观测到的周期, 返回按 Sharkovskii 序下一个必须存在的周期
 * 如果当前周期在 Sharkovskii 算法意义下"最大", 返回0
 *
 * 应用: 如果观测到周期3, 则必定存在所有其他周期
 */
int sharkovskii_next(int observed_period);

#endif /* CHAOS_BIFURCATION_H */
