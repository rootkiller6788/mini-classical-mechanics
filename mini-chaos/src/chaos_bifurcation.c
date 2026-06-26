/**
 * chaos_bifurcation.c — 分岔理论与分析实现
 *
 * 涵盖:
 *   L2: 不动点稳定性, 周期检测, 分岔分类
 *   L4: Feigenbaum普适性, Sharkovskii定理
 *   L5: 分岔图生成, 周期倍增点搜索
 *
 * 参考:
 *   Strogatz Ch.3, Ch.8, Ch.10
 *   Kuznetsov "Elements of Applied Bifurcation Theory" (2004)
 *   Feigenbaum (1978) "Quantitative Universality"
 */

#include "chaos_bifurcation.h"
#include <stdio.h>
#include <math.h>

/* ============================================================
 * L2: 不动点求解与稳定性
 * ============================================================ */

/**
 * 1D映射不动点Newton求解
 *
 * 求x使得f(x,p) = x, 即g(x)=f(x,p)-x=0
 *
 * Newton迭代: x_{n+1} = x_n - g(x_n)/g'(x_n)
 * g'(x) = f'(x,p) - 1
 *
 * 数值导数: g'(x) ≈ (g(x+dx) - g(x-dx)) / (2dx)
 *
 * 定理(Newton法的局部收敛):
 *   若x*是简单根(g'(x*)≠0)且初始猜测充分接近,
 *   则迭代二次收敛: |x_{n+1}-x*| ≤ C|x_n-x*|^2
 *
 * @param f         映射函数
 * @param p         参数值
 * @param x0        初始猜测
 * @param tol       收敛容忍度
 * @param max_iter  最大迭代次数
 * @param converged [out] 是否收敛
 * @return 不动点估计
 *
 * 复杂度: O(max_iter)
 * 参考: Burden & Faires §2.3
 */
double find_fixed_point(map1d_t f, double p, double x0,
                        double tol, int max_iter, int *converged)
{
    double x = x0;
    double dx = 1e-6;
    int iter;

    *converged = 0;
    for (iter = 0; iter < max_iter; iter++) {
        double fx = f(x, p);
        double gx = fx - x;
        if (fabs(gx) < tol) {
            *converged = 1;
            return x;
        }

        /* 数值导数g'(x) */
        double gp = (f(x + dx, p) - (x + dx) - (f(x - dx, p) - (x - dx))) / (2.0 * dx);
        if (fabs(gp) < 1e-15) break;  /* 导数接近零, 无法继续 */

        x = x - gx / gp;

        /* 保持在合理范围 */
        if (x < -100.0 || x > 100.0) break;
    }

    return x;
}

/**
 * 不动点稳定性判定(1D映射)
 *
 * 乘子: mu = f'(x*)
 *
 * 分类(Strogatz §10.1):
 *   |mu| < 1: 线性稳定 (稳定不动点)
 *   |mu| > 1: 线性不稳定
 *   |mu| = 1: 分岔点
 *     mu =  1: 鞍结点/跨临界/叉形分岔
 *     mu = -1: 周期倍增分岔
 *
 * 定理: 对于双曲不动点(|mu|≠1), 非线性映射在其附近
 *       拓扑共轭于其线性化(Hartman-Grobman定理的映射版本)
 *
 * @param f     映射函数
 * @param x_star 不动点
 * @param p     参数
 * @param dx    数值差分步长
 * @return 0=稳定, 1=不稳定, 2=临界(分岔点)
 *
 * 复杂度: O(1) + 2次f求值
 */
int fixed_point_stability(map1d_t f, double x_star, double p, double dx)
{
    double fp = f(x_star + dx, p);
    double fm = f(x_star - dx, p);
    double mu = (fp - fm) / (2.0 * dx);
    double abs_mu = fabs(mu);

    if (abs_mu < 0.9999) return 0;       /* 稳定 */
    if (abs_mu > 1.0001) return 1;       /* 不稳定 */
    return 2;                             /* 临界(分岔) */
}

/**
 * 平衡点稳定性分类(n维流系统)
 *
 * 特征值分类:
 *   所有Re(λ)<0, 所有Im(λ)=0:  稳定结点
 *   所有Re(λ)<0, 部分Im(λ)≠0:  稳定焦点
 *   所有Re(λ)>0, 所有Im(λ)=0:  不稳定结点
 *   所有Re(λ)>0, 部分Im(λ)≠0:  不稳定焦点
 *   Re(λ)符号不一:           鞍点
 *   Re(λ)=0, Im(λ)≠0:       中心
 *   某个Re(λ)=0:              退化
 *
 * 参考: Strogatz §5.2 "Classification of Linear Systems"
 */
FixedPointType classify_equilibrium(const double *eigenvalues, int dim)
{
    if (dim < 1 || !eigenvalues) return DEGENERATE;

    int has_pos = 0, has_neg = 0, has_zero = 0, has_imag = 0;
    int i;

    /* 注: 此简化版假设输入为实数特征值。
     * 对于复特征值, 需另行处理(当前函数用于从已分解的实特征值判断) */
    for (i = 0; i < dim; i++) {
        double re = eigenvalues[i];
        /* 简化: 假设输入为实数。复特征值应成对出现, 
           此处按实部正负判别 */
        if (fabs(re) < 1e-10) has_zero++;
        else if (re > 0.0) has_pos++;
        else has_neg++;
    }

    /* 补充: 复数对检测 — 在实际使用中, 调用者应提供实特征值 */
    if (has_zero > 0) return DEGENERATE;

    if (has_neg == dim) {
        return (has_imag > 0) ? STABLE_FOCUS : STABLE_NODE;
    }
    if (has_pos == dim) {
        return (has_imag > 0) ? UNSTABLE_FOCUS : UNSTABLE_NODE;
    }
    if (has_pos > 0 && has_neg > 0) return SADDLE_POINT;
    if (has_pos == 0 && has_neg == 0 && has_imag > 0) return CENTER;

    return DEGENERATE;
}

/* ============================================================
 * L5: 分岔图生成
 * ============================================================ */

/**
 * Logistic映射分岔图
 *
 * 算法:
 *   for each r in [r_min, r_max]:
 *     x = 0.5
 *     for i=1..n_transient: x = f(x, r)     # 丢弃瞬态
 *     for i=1..n_plot:     x = f(x, r); 记录(r, x)
 *
 * 知识点: 倍周期级联的可视化, 周期窗口(如r≈3.828的周期3窗口),
 *         Feigenbaum普适性的几何表现
 *
 * 复杂度: O(n_r * (n_transient + n_plot))
 */
BifurcationDiagram* logistic_bifurcation_diagram(double r_min, double r_max,
                                                  int n_r, int n_transient,
                                                  int n_plot)
{
    BifurcationDiagram *bd = (BifurcationDiagram*)chaos_calloc(1, sizeof(BifurcationDiagram));
    bd->n_params          = n_r;
    bd->n_steady_per_param = n_plot;
    bd->param_vals  = (double*)chaos_calloc(n_r, sizeof(double));
    bd->steady_vals = (double*)chaos_calloc((size_t)n_r * n_plot, sizeof(double));

    int i, j;
    double dr = (r_max - r_min) / (n_r - 1);

    for (i = 0; i < n_r; i++) {
        double r = r_min + i * dr;
        bd->param_vals[i] = r;

        double x = 0.5;

        /* 瞬态消除 */
        for (j = 0; j < n_transient; j++)
            x = logistic_map(x, r);

        /* 稳态采样 */
        for (j = 0; j < n_plot; j++) {
            x = logistic_map(x, r);
            bd->steady_vals[i * n_plot + j] = x;
        }
    }

    return bd;
}

/**
 * 通用1D映射分岔图
 *
 * @param f            映射函数f(x,p)
 * @param param_range  参数值数组
 * @param n_params     参数值个数
 * @param x0           各参数重用的初始值
 */
BifurcationDiagram* bifurcation_diagram_1d(map1d_t f,
                                            const double *param_range,
                                            int n_params,
                                            double x0,
                                            int n_transient, int n_plot)
{
    BifurcationDiagram *bd = (BifurcationDiagram*)chaos_calloc(1, sizeof(BifurcationDiagram));
    bd->n_params          = n_params;
    bd->n_steady_per_param = n_plot;
    bd->param_vals  = (double*)chaos_calloc(n_params, sizeof(double));
    bd->steady_vals = (double*)chaos_calloc((size_t)n_params * n_plot, sizeof(double));

    int i, j;
    for (i = 0; i < n_params; i++) {
        double p = param_range[i];
        bd->param_vals[i] = p;

        double x = x0;
        for (j = 0; j < n_transient; j++)
            x = f(x, p);

        for (j = 0; j < n_plot; j++) {
            x = f(x, p);
            bd->steady_vals[i * n_plot + j] = x;
        }
    }

    return bd;
}

/**
 * 2D映射的分岔图(投影到指定变量)
 *
 * @param f        2D映射函数
 * @param var_idx  记录的变量索引(0=x, 1=y)
 * @param p1_fixed 固定参数(如Henon的b=0.3, 另一个参数扫描)
 */
BifurcationDiagram* bifurcation_diagram_2d(map2d_t f,
                                            const double *param_range,
                                            int n_params,
                                            double x0, double y0,
                                            double p1_fixed,
                                            int var_idx,
                                            int n_transient, int n_plot)
{
    BifurcationDiagram *bd = (BifurcationDiagram*)chaos_calloc(1, sizeof(BifurcationDiagram));
    bd->n_params          = n_params;
    bd->n_steady_per_param = n_plot;
    bd->param_vals  = (double*)chaos_calloc(n_params, sizeof(double));
    bd->steady_vals = (double*)chaos_calloc((size_t)n_params * n_plot, sizeof(double));

    int i, j;
    double xy[2];

    for (i = 0; i < n_params; i++) {
        double p = param_range[i];
        bd->param_vals[i] = p;

        double x = x0, y = y0;

        for (j = 0; j < n_transient; j++) {
            f(xy, x, y, p, p1_fixed);
            x = xy[0]; y = xy[1];
        }

        for (j = 0; j < n_plot; j++) {
            f(xy, x, y, p, p1_fixed);
            x = xy[0]; y = xy[1];
            bd->steady_vals[i * n_plot + j] = (var_idx == 0) ? x : y;
        }
    }

    return bd;
}

/* ============================================================
 * L2: 周期检测
 * ============================================================ */

/**
 * 检测离散映射的周期
 *
 * 方法: 寻找重复值
 * 周期p满足: x_{n+p} = x_n (在tol范围内)
 *
 * @param period [out] 检测到的周期(0=未检测到, 1=不动点, 2=周期2, ...)
 * @return 0=成功, -1=未检测到
 *
 * 知识点: 周期轨道定义, 数值检测的tol敏感性问题
 * 复杂度: O(max_iter^2) (朴素查找)
 */
int detect_period_discrete(map1d_t f, double x0, double param,
                            int max_iter, double tol, int *period)
{
    int max_visit = max_iter;
    double *visited = (double*)chaos_calloc(max_visit, sizeof(double));
    double x = x0;
    int n_visited = 0;
    int i;

    *period = 0;

    for (i = 0; i < max_iter; i++) {
        x = f(x, param);

        /* 与已访问值比较 */
        int j;
        for (j = 0; j < n_visited; j++) {
            if (fabs(x - visited[j]) < tol) {
                *period = n_visited - j;
                free(visited);
                return 0;
            }
        }

        if (n_visited < max_visit) {
            visited[n_visited++] = x;
        }
    }

    free(visited);
    return -1;
}

/**
 * 在给定参数处查找周期(自动去瞬态)
 *
 * @return 周期数, 0=未找到
 */
int find_period_at_param(map1d_t f, double param, double x0,
                          int n_transient, int n_iter, double tol)
{
    int i;
    double x = x0;

    for (i = 0; i < n_transient; i++)
        x = f(x, param);

    int period;
    if (detect_period_discrete(f, x, param, n_iter, tol, &period) == 0)
        return period;
    return 0;
}

/* ============================================================
 * L4: Feigenbaum普适性与Sharkovskii定理
 * ============================================================ */

/**
 * 寻找周期倍增分岔点序列
 *
 * 从r_start开始搜索, 递增r直至检测到周期2,4,8,...
 *
 * 算法:
 *   r = r_start, dr = dr_init
 *   for each target period p=2,4,8,...:
 *     while period(r) != p: r += dr
 *     record r, reduce dr
 *
 * @param r_vals      [out] 分岔点数组(调用者分配, 长度n_periods)
 * @param n_periods   要找的分岔点个数
 * @param r_start     起始参数值(应略小于第一个分岔点)
 * @param dr_init     初始步长
 * @param tol         周期检测容忍度
 * @return 实际找到的分岔点个数
 *
 * 知识点: 周期倍增级联, 通向混沌的倍周期路径
 * 复杂度: O(n_periods * steps)
 */
int find_period_doubling_points(double *r_vals, int n_periods,
                                 double r_start, double dr_init,
                                 double tol)
{
    double r = r_start;
    double dr = dr_init;
    int found = 0;
    int p;

    for (p = 0; p < n_periods; p++) {
        int target_period = 1 << (p + 1);  /* 2, 4, 8, ... */

        int max_attempts = 10000;
        int attempt = 0;
        while (attempt < max_attempts) {
            int period = find_period_at_param(logistic_map, r, 0.5, 500, 2000, tol);
            if (period == target_period) {
                r_vals[found++] = r;
                dr /= 5.0;
                break;
            }
            r += dr;
            if (r > 4.0) break;  /* Logistic映射参数上限 */
            attempt++;
        } /* while */

        if (attempt >= max_attempts) break;
    }

    return found;
}

/**
 * 从分岔点序列估计Feigenbaum delta
 *
 * delta_n = (r_n - r_{n-1}) / (r_{n+1} - r_n)
 *
 * 理论极限: lim_{n->inf} delta_n = FEIGENBAUM_DELTA ≈ 4.669202
 *
 * 示例(Logistic映射):
 *   n=1(r=3): delta_1 = (3.44949-3) / (3.54409-3.44949) ≈ 4.751
 *   n=2(r=3.44949): delta_2 ≈ 4.656
 *   n=3: delta_3 ≈ 4.668
 *   n=4: delta_4 ≈ 4.669
 *
 * 知识点: 普适性 — 对所有单峰映射, delta相同
 *         重正化群解释: 函数空间中的不动点
 *
 * 参考: Feigenbaum (1978), Strogatz §10.7
 *       Coullet & Tresser (1978) 独立发现
 *
 * @param deltas [out] delta估计值序列
 * @return 最后一个(最佳)估计值
 */
double estimate_feigenbaum_delta(const double *r_vals, int n,
                                  double *deltas)
{
    if (n < 3) return 0.0;

    int i;
    for (i = 1; i < n - 1; i++) {
        double num = r_vals[i] - r_vals[i - 1];
        double den = r_vals[i + 1] - r_vals[i];
        if (fabs(den) > 1e-15)
            deltas[i - 1] = num / den;
        else
            deltas[i - 1] = 0.0;
    }

    return deltas[n - 3];
}

/**
 * 分类分岔类型(1D映射)
 *
 * 条件(Kuznetsov §2.1):
 *   saddle-node:   f(x*,r*)=x*, f_x=1, f_xx≠0, f_r≠0
 *   transcritical: f(x*,r*)=x*, f_x=1, f_xx≠0, f(0,r)=0
 *   pitchfork:     f(x*,r*)=x*, f_x=1, f_xx=0, f_xxx≠0 (对称中)
 *   period-doubling: f^2(x*,r*)=x*, (f^2)_x=-1
 *
 * @param x_star   不动点
 * @param r_star   分岔参数
 * @param dx       数值差分步长
 * @return BifurcationType
 *
 * 复杂度: 需要多次f求值以估计高阶导数
 */
BifurcationType classify_bifurcation(map1d_t f, double x_star,
                                      double r_star, double dx)
{
    /* 一阶导数 f_x */
    double fp = f(x_star + dx, r_star);
    double fm = f(x_star - dx, r_star);
    double fx = (fp - fm) / (2.0 * dx);

    /* 在分岔点应有 |f_x| ≈ 1 */
    if (fabs(fabs(fx) - 1.0) > 0.01) return BIFURCATION_NONE;

    int is_period_doubling = (fx < -0.9);

    /* 二阶导数 f_xx */
    double f0 = f(x_star, r_star);
    double f_xx = (fp - 2.0 * f0 + fm) / (dx * dx);

    if (fabs(f_xx) > 1e-5 && !is_period_doubling) {
        /* 检查跨临界条件: f(0,r)=0 对所有r */
        double f_zero = f(0.0, r_star);
        if (fabs(f_zero) < 1e-10)
            return BIFURCATION_TRANSCRITICAL;
        else
            return BIFURCATION_SADDLE_NODE;
    }

    if (fabs(f_xx) < 1e-5 && !is_period_doubling) {
        return BIFURCATION_PITCHFORK;
    }

    if (is_period_doubling)
        return BIFURCATION_PERIOD_DOUBLING;

    /* Hopf分岔需要复特征值, 一维映射无Hopf */
    return BIFURCATION_NONE;
}

/**
 * Sharkovskii序 — 下一级周期
 *
 * 定理(Sharkovskii 1964): 若f有周期k, 则对所有k ▷ m, f必有周期m
 *
 * Sharkovskii序:
 *   3 ▷ 5 ▷ 7 ▷ 9 ▷ ... (所有大于1的奇数)
 *   ▷ 2·3 ▷ 2·5 ▷ 2·7 ▷ ... (2*奇数)
 *   ▷ 2²·3 ▷ 2²·5 ▷ ... (4*奇数)
 *   ▷ ... ▷ 2³ ▷ 2² ▷ 2 ▷ 1
 *
 * 推论(Li-Yorke 1975): 周期3 ⇒ 所有周期 ⇒ 混沌
 *
 * 此函数返回给定周期的Sharkovskii下一级
 *
 * @param observed_period 观测到的周期
 * @return 下一个"必然存在"的周期, 0=当前已是最大
 */
int sharkovskii_next(int observed_period)
{
    if (observed_period <= 1) return 0;

    /* 如果周期是奇数且>1, 下一个是2*周期 */
    if (observed_period % 2 == 1 && observed_period > 1) {
        if (observed_period == 3) return 0;  /* 3是Sharkovskii序的极大元! */
        return observed_period - 2;  /* 递减奇数序列: 7→5→3 */
    }

    /* 如果周期是2*(奇数), 下一个*/
    if (observed_period % 2 == 0) {
        int factor = 1;
        int rem = observed_period;
        while (rem % 2 == 0) {
            factor *= 2;
            rem /= 2;
        }
        if (rem == 1) {
            /* 2^n: 下一个是2^{n-1} */
            return observed_period / 2;
        } else {
            /* 2^n * (奇数): 下一个是2^n * (下一个奇数) */
            if (rem == 3) {
                /* 2^n*3 → 2^{n+1}*odd */
                return factor * 2;
            }
            return observed_period - 2 * factor;  /* 下一个更小的奇数 */
        }
    }

    return 0;
}
