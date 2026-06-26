/**
 * chaos_maps.h — 离散时间混沌映射
 *
 * 一维映射:
 *   Logistic — 人口动力学经典, 倍周期分岔进入混沌
 *   Tent     — 分段线性, 拓扑共轭于移位映射
 *   Sine     — 圆映射特例
 *   Gauss    — 连分数变换, 遍历理论
 *
 * 二维映射:
 *   Hénon    — 面积收缩, 奇怪吸引子
 *   Standard (Chirikov-Taylor) — Hamiltonian 混沌
 *   Lozi     — 分段线性 Hénon 变体
 *   Ikeda    — 光学环形腔
 *   Circle   — 锁频与 Arnold 舌头
 *
 * 参考: Strogatz Ch.10, May (1976), Chirikov (1979)
 *       Hénon "A Two-dimensional Mapping with a Strange Attractor" (1976)
 */

#ifndef CHAOS_MAPS_H
#define CHAOS_MAPS_H

#include "chaos.h"

/* ──────────────────────────────────────────────────────────
 * L6: 一维离散映射
 * ────────────────────────────────────────────────────────── */

/**
 * Logistic 映射: x_{n+1} = r x_n (1 - x_n)
 *
 * 定义域: x ∈ [0,1], r ∈ [0,4]
 * 动力学:
 *   0 < r ≤ 1:  不动点 x*=0 (稳定)
 *   1 < r ≤ 3:  不动点 x*=1-1/r (稳定)
 *   3 < r ≤ 3.449: 周期2
 *   3.449 < r ≤ 3.544: 周期4
 *   ...
 *   r ≈ 3.569946: Feigenbaum 点 (混沌开始)
 *   r = 4:  完全混沌, λ=log 2
 *
 * 参考: May (1976), Strogatz §10.1–10.3
 */
double logistic_map(double x, double r);

/**
 * 立方映射: x_{n+1} = r x_n - x_n³
 *
 * 知识点: 对称性破缺, 双稳态
 */
double cubic_map(double x, double r);

/**
 * Sine 映射: x_{n+1} = r sin(π x_n)
 *
 * 知识点: 与 Logistic 映射的拓扑共轭
 */
double sine_map(double x, double r);

/**
 * Tent 映射: x_{n+1} = μ min(x, 1-x)
 *
 * μ=2: 混沌, λ=log 2, 不变测度均匀
 * 知识点: 分段线性, Markov 分割, 符号动力学
 */
double tent_map(double x, double mu);

/**
 * Gauss 映射 (连分数映射)
 *
 * x_{n+1} = 1/x_n mod 1, x∈[0,1]
 *
 * 不变测度: ρ(x) = 1/((1+x)log 2)
 * 知识点: 遍历理论, 连分数展开
 */
double gauss_map(double x);

/**
 * Bernoulli 移位: x_{n+1} = 2x_n mod 1
 *
 * 符号动力学表示: s_n = ⌊2x_n⌋ ∈ {0,1}
 * 知识点: 确定论混沌 = 符号序列的随机性
 */
double bernoulli_shift(double x);

/* ──────────────────────────────────────────────────────────
 * L6: 二维离散映射
 * ────────────────────────────────────────────────────────── */

/**
 * Hénon 映射: x_{n+1} = 1 - a x_n² + y_n
 *             y_{n+1} = b x_n
 *
 * 经典参数: a=1.4, b=0.3
 * Jacobian: det J = -b (面积收缩)
 * 知识点: 奇怪吸引子, 马蹄映射, Smale 马蹄
 *
 * 参考: Hénon (1976), Strogatz §12.2
 */
void henon_map(double *xy_out, double x, double y, double a, double b);

/**
 * Lozi 映射: x_{n+1} = 1 - a|x_n| + y_n
 *            y_{n+1} = b x_n
 *
 * 知识点: 分段线性奇怪吸引子, 可严格证明混沌
 */
void lozi_map(double *xy_out, double x, double y, double a, double b);

/**
 * 标准映射 (Chirikov-Taylor)
 *
 * p_{n+1} = p_n + K sin θ_n
 * θ_{n+1} = θ_n + p_{n+1}  (mod 2π)
 *
 * K ≪ 1: 可积, K ≈ 1: KAM 环面破坏, K ≫ 1: 全局混沌
 * 知识点: Hamiltonian 混沌, KAM 定理, 共振重叠准则
 *
 * 参考: Chirikov (1979), Lichtenberg & Lieberman (1992)
 */
void standard_map(double *thetap_out, double theta, double p, double K);

/**
 * Ikeda 映射 (光学环形腔)
 *
 * E_{n+1} = A + B E_n exp(i(φ - p/(1+|E|²)))
 *
 * 知识点: 光学双稳态, 延迟反馈混沌
 */
void ikeda_map(double *xy_out, double x, double y,
               double A, double B, double phi, double p);

/**
 * 圆映射: θ_{n+1} = θ_n + Ω - (K/2π) sin(2πθ_n)  (mod 1)
 *
 * 知识点: 锁频, Arnold 舌头, 拟周期 → 混沌 (通过频率锁)
 * 魔鬼阶梯: 锁频区间在有理数处形成完整测度集
 *
 * 参考: Arnold (1965), Strogatz §8.6
 */
double circle_map(double theta, double Omega, double K);

/**
 * Arnold 猫映射: x_{n+1} = 2x_n + y_n (mod 1)
 *                y_{n+1} = x_n + y_n (mod 1)
 *
 * Lyapunov: λ₁,₂ = log((3±√5)/2) ≈ ±0.962
 * 知识点: Anosov 微分同胚, 遍历性, 混合性
 */
void arnold_cat_map(double *xy_out, double x, double y);

/**
 * Baker 映射: x_{n+1} = 2x_n mod 1
 *             y_{n+1} = y_n/2 + ⌊2x_n⌋/2
 *
 * 知识点: Bernoulli 过程的二维嵌入
 */
void bakers_map(double *xy_out, double x, double y);

/**
 * Gingerbreadman 映射: x_{n+1} = 1 - y_n + |x_n|
 *                       y_{n+1} = x_n
 *
 * 知识点: 保守混沌, 六边形周期岛
 */
void gingerbreadman_map(double *xy_out, double x, double y);

/* ──────────────────────────────────────────────────────────
 * L5: 迭代工具
 * ────────────────────────────────────────────────────────── */

/**
 * 离散映射函数签名: x_{n+1} = f(x_n, param)
 */
typedef double (*map1d_t)(double x, double param);

/**
 * 二维映射函数签名
 */
typedef void (*map2d_t)(double *out, double x, double y, double p1, double p2);

/**
 * 迭代映射并收集轨迹
 *
 * @param n_transient  丢弃的前若干个点 (瞬态消除)
 * @param n_iter       收集的迭代点数
 * @return 长度为 n_iter 的数组 (调用者 free)
 */
double* iterate_map1d(map1d_t f, double x0, double param,
                      int n_transient, int n_iter);

/**
 * 迭代二维映射并收集轨迹
 *
 * @return 长度为 2*n_iter 的数组 [x0,y0,x1,y1,...] (调用者 free)
 */
double* iterate_map2d(map2d_t f, double x0, double y0,
                      double p1, double p2,
                      int n_transient, int n_iter);

/**
 * 透镜映射 (Lens Map) — 可严格证明的混沌映射
 *
 * x_{n+1} = (a x_n) mod 1
 * 当 a>1: 正 Lyapunov, 混沌
 */
double lens_map(double x, double a);

#endif /* CHAOS_MAPS_H */
