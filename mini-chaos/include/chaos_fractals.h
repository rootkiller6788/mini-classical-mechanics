/**
 * chaos_fractals.h — 分形几何
 *
 * 分形定义 (Mandelbrot):
 *   集合其 Hausdorff-Besicovitch 维数严格大于拓扑维数
 *
 * 核心概念:
 *   - 自相似性 (self-similarity): 部分与整体相似
 *   - 迭代函数系统 (IFS)
 *   - 逃逸时间算法 (escape-time)
 *
 * Mandelbrot 集:
 *   M = {c ∈ ℂ : z_{n+1}=z_n²+c 不趋于无穷, z₀=0}
 *   边界具有无穷复杂度, 是连通集 (Douady & Hubbard 1982)
 *
 * Julia 集:
 *   J_c = {z₀ ∈ ℂ : z_{n+1}=z_n²+c 不趋于无穷}
 *   Fatou 集: 补集
 *
 * 参考:
 *   Mandelbrot "The Fractal Geometry of Nature" (1983)
 *   Peitgen & Richter "The Beauty of Fractals" (1986)
 *   Falconer "Fractal Geometry: Mathematical Foundations" (2013)
 */

#ifndef CHAOS_FRACTALS_H
#define CHAOS_FRACTALS_H

#include "chaos.h"

/* ──────────────────────────────────────────────────────────
 * L2: 逃逸时间算法 — 复动力学基础
 * ────────────────────────────────────────────────────────── */

/**
 * Mandelbrot 集迭代
 *
 * 动力学: z₀ = 0,  z_{n+1} = z_n² + c
 *
 * 逃逸准则: |z| > 2 ⇒ 发散
 * 基本性质: 若 |z_n| > 2, 则轨道必逃逸到无穷
 *
 * @param cx, cy    复参数 c 的实部和虚部
 * @param max_iter  最大迭代次数
 * @return 逃逸时的迭代次数 (max_iter = 未逃逸, 在M集内)
 *
 * 复杂度: O(max_iter)
 */
int mandelbrot_iter(double cx, double cy, int max_iter);

/**
 * Julia 集迭代
 *
 * 动力学: z_{n+1} = z_n² + c (固定 c, 变 z₀)
 *
 * @param zx, zy    初始值 z₀ 的实部和虚部
 * @param cx, cy    Julia 常数 c
 * @param max_iter  最大迭代次数
 * @return 逃逸迭代次数
 */
int julia_iter(double zx, double zy, double cx, double cy, int max_iter);

/**
 * 距离估计法 — 平滑着色 Mandelbrot
 *
 * 原理: 利用导数估计到M集边界的距离
 * d(c) = ½|z_n| log|z_n| / |z'_n|
 *
 * 参考: Peitgen & Richter §4.5
 *
 * @return 距离估计值 (0 = 在M集内)
 */
double mandelbrot_distance(double cx, double cy, int max_iter);

/* ──────────────────────────────────────────────────────────
 * L6: 分形图像生成
 * ────────────────────────────────────────────────────────── */

/**
 * Mandelbrot 集图像
 *
 * @param x_min, x_max  实部范围
 * @param y_min, y_max  虚部范围
 * @param width, height 图像分辨率
 * @param max_iter      最大迭代次数
 * @return FractalImage (调用者 free)
 */
FractalImage* mandelbrot_set_image(double x_min, double x_max,
                                    double y_min, double y_max,
                                    int width, int height, int max_iter);

/**
 * Julia 集图像
 *
 * @param cx, cy  Julia 常数
 */
FractalImage* julia_set_image(double x_min, double x_max,
                               double y_min, double y_max,
                               int width, int height,
                               double cx, double cy, int max_iter);

/* ──────────────────────────────────────────────────────────
 * L1: 定义 — Newton 分形
 *
 * 方程: z³ - 1 = 0, 三个根: 1, e^{2πi/3}, e^{4πi/3}
 * 迭代: z_{n+1} = z_n - (z_n³ - 1)/(3z_n²) = ⅔z_n + 1/(3z_n²)
 *
 * 吸引盆: 复平面上三个根的吸引域具有分形边界 (Wada 性质)
 * ────────────────────────────────────────────────────────── */

/**
 * Newton 分形迭代
 *
 * @return 收敛到的根的索引 (1/2/3) 或 0 (未收敛)
 */
int newton_fractal(double zx, double zy, int max_iter, double tol);

/**
 * Newton 分形图像
 */
FractalImage* newton_fractal_image(double x_min, double x_max,
                                    double y_min, double y_max,
                                    int width, int height,
                                    int max_iter, double tol);

/* ──────────────────────────────────────────────────────────
 * L6: 经典分形 — 确定性迭代
 * ────────────────────────────────────────────────────────── */

/**
 * Sierpinski 三角 (混沌游戏法)
 *
 * 算法:
 *   1. 设定三个顶点
 *   2. 从任意点开始
 *   3. 每步随机选择一个顶点, 移动到中点
 *   4. 记录位置
 *
 * @param points   [out] 长度为 2*n 的数组 [x1,y1,x2,y2,...]
 * @param n        点数
 *
 * 参考: Barnsley "Fractals Everywhere" (1993)
 */
void sierpinski_triangle(double *points, int n);

/**
 * Koch 雪花曲线
 *
 * L-system: F → F+F--F+F, θ=60°
 * 每代长度 = 3×4^n 段
 *
 * 维数: D = log(4)/log(3) ≈ 1.262
 *
 * @param seg_x, seg_y  [out] 线段端点坐标
 * @param max_segs      分配的最大线段数 (= 2·4^n + 1)
 * @param n_iter        迭代次数
 * @return 实际线段数
 */
int koch_snowflake(double *seg_x, double *seg_y, int max_segs, int n_iter);

/* ──────────────────────────────────────────────────────────
 * L5: 分形维数 — 盒计数法
 *
 * 定义: D₀ = lim_{ε→0} log N(ε) / log(1/ε)
 *
 * N(ε) = 覆盖集合所需的最小 ε-盒数
 *
 * 实现: 对一组递减的 ε 值, 统计非空盒数, 线性回归 log N ~ log(1/ε)
 * ────────────────────────────────────────────────────────── */

/**
 * 盒计数维数
 *
 * @param points      点坐标数组 [x1,y1,x2,y2,...]
 * @param n_points    点数
 * @param epsilons    尺度序列
 * @param n_eps       尺度个数
 * @return 维数估计 (最小二乘斜率)
 */
double box_counting_dimension(const double *points, int n_points,
                               const double *epsilons, int n_eps);

/**
 * 自相关函数 C(ε) 维数估计
 *
 * C(ε) = (1/N²) Σ_{i,j} Θ(ε - |x_i - x_j|) ≈ ε^D₂
 * D₂ = 关联维数 (Grassberger & Procaccia 1983)
 *
 * @return D₂ 估计值
 */
double correlation_dimension_gp(const double *points, int n_points,
                                 int dim, const double *epsilons,
                                 int n_eps);

/**
 * Koch 雪花理论维数
 */
double koch_dimension(void);

/**
 * Sierpinski 三角理论维数
 */
double sierpinski_dimension(void);

/**
 * Cantor 集理论维数: D = log(2)/log(3) ≈ 0.6309
 */
double cantor_dimension(void);

#endif /* CHAOS_FRACTALS_H */
