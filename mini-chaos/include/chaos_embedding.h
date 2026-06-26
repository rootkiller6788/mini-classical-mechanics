/**
 * chaos_embedding.h — 相空间重构与时间序列分析
 *
 * Takens 嵌入定理 (1981):
 *   对于 d 维紧致流形上的光滑动力学系统, 通用观测函数 h,
 *   延迟嵌入映射 Φ(x) = (h(x), h(F^{-τ}(x)), ..., h(F^{-(m-1)τ}(x)))
 *   是嵌入 (embedding), 只要 m > 2d.
 *
 *   → 从标量时间序列可重构与原系统微分同胚的相空间
 *
 * 参考:
 *   Takens (1981) "Detecting Strange Attractors in Turbulence"
 *   Sauer, Yorke, Casdagli (1991) "Embedology"
 *   Kantz & Schreiber "Nonlinear Time Series Analysis" (2004)
 */

#ifndef CHAOS_EMBEDDING_H
#define CHAOS_EMBEDDING_H

#include "chaos.h"

/* ──────────────────────────────────────────────────────────
 * L3: 相空间重构 — 延迟嵌入
 * ────────────────────────────────────────────────────────── */

/**
 * 时间延迟嵌入
 *
 * X(t) = [s(t), s(t-τ), s(t-2τ), ..., s(t-(m-1)τ)]
 *
 * @param signal   标量时间序列
 * @param n        信号长度
 * @param tau      时间延迟 (采样步数)
 * @param m        嵌入维数
 * @return Embedding 结构 (调用者 free)
 *
 * 知识点: Takens 定理, 延迟坐标, 重构相空间
 *
 * 参考: Takens (1981), Kantz & Schreiber §3.2
 */
Embedding* time_delay_embedding(const double *signal, int n,
                                 int tau, int m);

/* ──────────────────────────────────────────────────────────
 * L5: 嵌入参数选择
 * ────────────────────────────────────────────────────────── */

/**
 * 平均互信息 (AMI) — 确定最优延迟 τ
 *
 * I(τ) = Σ P(s_n, s_{n+τ}) log[P(s_n, s_{n+τ}) / (P(s_n)P(s_{n+τ}))]
 *
 * 首个 AMI 极小值 → 最优 τ
 *
 * @param signal  信号
 * @param n       信号长度
 * @param tau     延迟
 * @param n_bins  直方图箱数
 * @return 互信息值 (单位: nat)
 *
 * 参考: Fraser & Swinney (1986), Kantz & Schreiber §3.3
 */
double average_mutual_information(const double *signal, int n,
                                   int tau, int n_bins);

/**
 * 寻找最优延迟 τ (首个 AMI 极小值)
 *
 * @param tau_max  搜索上界
 * @return 最优 τ
 */
int find_optimal_tau(const double *signal, int n, int tau_max);

/**
 * 假近邻法 (FNN) — 确定最小嵌入维数 m
 *
 * 原理: 当嵌入维数不足时, 在低维空间中看似邻近的点
 *       在高维中可能变得远离 ("假近邻")
 *
 * 准则: R_i = |x_{i+mτ} - x_{NN(i)+mτ}| / |x_i - x_{NN(i)}|
 *       如果 R_i > R_tol, 两点是假近邻
 *
 * @param signal  信号
 * @param n       信号长度
 * @param tau     延迟
 * @param m_max   最大测试维数
 * @param R_tol   阈值 (典型值 10)
 * @param fnn_ratio [out] 每个 m 的假近邻比例
 * @return 假近邻比例降至 ~0 时的 m 值
 *
 * 参考: Kennel, Brown, Abarbanel (1992)
 */
int false_nearest_neighbors(const double *signal, int n, int tau,
                             int m_max, double R_tol,
                             double *fnn_ratio);

/**
 * 关联维数 (Grassberger-Procaccia 算法)
 *
 * 关联积分: C(ε) = (2/(N(N-1))) Σ_{i<j} Θ(ε - |x_i - x_j|)
 * 关联维数: D₂ = lim_{ε→0} log C(ε) / log ε
 *
 * @param points  延迟嵌入点集
 * @param n       点数
 * @param dim     每点维度
 * @param epsilons   ε 值数组
 * @param n_eps      ε 值个数
 * @return D₂ 估计
 *
 * 参考: Grassberger & Procaccia (1983)
 */
double correlation_dimension(const double *points, int n, int dim,
                              const double *epsilons, int n_eps);

/**
 * 替代数据生成 (IAAFT 风格)
 *
 * 目的: 零假设检验 — "数据来自线性高斯过程"
 *
 * 方法: 随机重排保持幅度分布, 用 Fourier 变换调整功率谱
 *
 * @param surrogate [out] 替代数据 (长度 n)
 * @param signal    [in]  原始信号 (长度 n)
 *
 * 参考: Schreiber & Schmitz (1996)
 */
void surrogate_data_shuffle(double *surrogate, const double *signal, int n);

/**
 * 非线性预测误差 — 确定性 vs 随机性判别
 *
 * 如果数据是确定论混沌:
 *   用嵌入空间中近邻的外推可预测未来值
 *   预测误差随嵌入维数 m 的改善而降低
 *
 * 如果数据是随机噪声:
 *   近邻外推不优于均值预测
 *
 * @param T_pred  向前预测步数
 * @return 归一化预测误差
 *
 * 参考: Sugihara & May (1990)
 */
double nonlinear_prediction_error(const double *signal, int n,
                                   int tau, int m, int T_pred);

#endif /* CHAOS_EMBEDDING_H */
