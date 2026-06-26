/**
 * chaos_lyapunov.c — Lyapunov指数算法实现
 *
 * 定义: lambda = lim_{t->inf} (1/t) * log(||delta x(t)|| / ||delta x(0)||)
 *
 * 算法:
 *   两粒子法 — 最大Lyapunov指数
 *   Wolf算法  — 改进的最大Lyapunov
 *   Benettin法 — 完整Lyapunov谱 + Kaplan-Yorke维数
 *
 * 参考:
 *   Benettin, Galgani, Strelcyn (1980) Meccanica 15:9, 15:21
 *   Wolf, Swift, Swinney, Vastano (1985) Physica D 16:285
 *   Strogatz §9.3
 */

#include "chaos_lyapunov.h"
#include "chaos_maps.h"
#include <stdio.h>
#include <math.h>

/* ──────────────────────────────────────────────────────────
 * L3: 数学结构 — Jacobian矩阵与切空间
 * ────────────────────────────────────────────────────────── */

/**
 * 数值Jacobian矩阵(中心差分)
 *
 * J_{ij} = df_i/dx_j ≈ [f_i(x+eps*e_j) - f_i(x-eps*e_j)] / (2*eps)
 *
 * 用于: Lyapunov指数计算, 平衡点稳定性分析, 分岔检测
 *
 * 复杂度: O(dim^2 * dim) = O(n^3) — 每个元素需要两次RHS求值
 *
 * 参考: Press et al. "Numerical Recipes" §5.7
 */
void jacobian_numerical(double *J, flow_rhs_t f, const double *x,
                        double t, int dim, const double *p,
                        double eps)
{
    int i, j;
    double *xp = (double*)chaos_calloc(dim, sizeof(double));
    double *xm = (double*)chaos_calloc(dim, sizeof(double));
    double *fp = (double*)chaos_calloc(dim, sizeof(double));
    double *fm = (double*)chaos_calloc(dim, sizeof(double));

    for (j = 0; j < dim; j++) {
        for (i = 0; i < dim; i++) {
            xp[i] = x[i];
            xm[i] = x[i];
        }
        xp[j] += eps;
        xm[j] -= eps;
        f(fp, xp, t, p);
        f(fm, xm, t, p);
        for (i = 0; i < dim; i++)
            J[i * dim + j] = (fp[i] - fm[i]) / (2.0 * eps);
    }

    free(xp); free(xm); free(fp); free(fm);
}

/**
 * 离散映射的Jacobian
 *
 * J_{ij} = df_i/dx_j
 *
 * @param map  映射函数: out = map(in), out和in均有dim维
 */
void jacobian_map(double *J, void (*map)(double*, const double*),
                  const double *x, int dim, double eps)
{
    int i, j;
    double *xp = (double*)chaos_calloc(dim, sizeof(double));
    double *xm = (double*)chaos_calloc(dim, sizeof(double));
    double *fp = (double*)chaos_calloc(dim, sizeof(double));
    double *fm = (double*)chaos_calloc(dim, sizeof(double));

    for (j = 0; j < dim; j++) {
        for (i = 0; i < dim; i++) {
            xp[i] = x[i];
            xm[i] = x[i];
        }
        xp[j] += eps;
        xm[j] -= eps;
        map(fp, xp);
        map(fm, xm);
        for (i = 0; i < dim; i++)
            J[i * dim + j] = (fp[i] - fm[i]) / (2.0 * eps);
    }

    free(xp); free(xm); free(fp); free(fm);
}

/* ──────────────────────────────────────────────────────────
 * L5: 两粒子法 — 最大Lyapunov指数
 * ────────────────────────────────────────────────────────── */

/**
 * 两粒子法
 *
 * 原理: 比较邻近轨道的指数分离率
 *
 * 算法步骤 (Strogatz §9.3):
 *   1. xA(0) = x0, xB(0) = x0 + d0*u (u为任意单位向量)
 *   2. 积分一步: xA, xB
 *   3. 计算当前距离: d = |xB - xA|
 *   4. 累计: sum += log(d/d0)
 *   5. 重正化: xB = xA + (d0/d)*(xB - xA)
 *   6. 重复
 *
 *   lambda_max = sum / (count * dt)
 *
 * 优点: 最直观, 易于理解和验证
 * 缺点: 只给出最大指数, 扰动方向可能偏转
 *
 * @param d0   初始扰动大小(推荐1e-8)
 * @param t_end 积分总时间
 * @param dt   积分步长
 * @return 最大Lyapunov指数
 *
 * 复杂度: O(n_steps * dim)
 */
double lyapunov_two_particle(flow_rhs_t f, const double *x0,
                              int dim, const double *p,
                              double t_end, double dt, double d0)
{
    int i, n_steps = (int)(t_end / dt);
    if (n_steps < 1) n_steps = 1;

    double *xA  = (double*)chaos_calloc(dim, sizeof(double));
    double *xB  = (double*)chaos_calloc(dim, sizeof(double));
    double *xAn = (double*)chaos_calloc(dim, sizeof(double));
    double *xBn = (double*)chaos_calloc(dim, sizeof(double));

    for (i = 0; i < dim; i++) xA[i] = x0[i];

    /* 生成随机初始扰动(用简单LCG) */
    unsigned int seed = 12345;
    for (i = 0; i < dim; i++) {
        seed = (1103515245 * seed + 12345) & 0x7fffffff;
        xB[i] = xA[i] + d0 * ((double)seed / 0x7fffffff - 0.5);
    }
    /* 归一化到d0 */
    double norm = 0.0;
    for (i = 0; i < dim; i++) {
        double diff = xB[i] - xA[i];
        norm += diff * diff;
    }
    norm = sqrt(norm);
    if (norm > 1e-30) {
        for (i = 0; i < dim; i++)
            xB[i] = xA[i] + (xB[i] - xA[i]) * d0 / norm;
    }

    double log_sum = 0.0;
    int count = 0;
    double t = 0.0;
    int step;

    for (step = 0; step < n_steps; step++) {
        rk4_step(xAn, f, xA, t, dt, dim, p);
        rk4_step(xBn, f, xB, t, dt, dim, p);

        /* 计算分离距离 */
        double d = 0.0;
        for (i = 0; i < dim; i++) {
            double diff = xBn[i] - xAn[i];
            d += diff * diff;
        }
        d = sqrt(d);

        if (d > 1e-15) {
            log_sum += log(d / d0);
            count++;

            /* 重正化: xB = xAn + (d0/d)*(xB - xA) */
            for (i = 0; i < dim; i++)
                xBn[i] = xAn[i] + (xBn[i] - xAn[i]) * d0 / d;
        }

        for (i = 0; i < dim; i++) {
            xA[i] = xAn[i];
            xB[i] = xBn[i];
        }
        t += dt;
    }

    free(xA); free(xB); free(xAn); free(xBn);

    return (count > 0) ? log_sum / (count * dt) : 0.0;
}

/**
 * Wolf算法 — 改进的最大Lyapunov指数
 *
 * 与两粒子法类似, 但仅在每隔evol_steps步时
 * 观测并重正化, 避免频繁重正化带来的数值漂移
 *
 * 参考: Wolf et al. (1985) "Determining Lyapunov Exponents from a Time Series"
 */
double lyapunov_wolf(flow_rhs_t f, const double *x0,
                     int dim, const double *p,
                     double t_end, double dt,
                     double d0, int evol_steps)
{
    int i, n_steps = (int)(t_end / dt);
    if (n_steps < 1) n_steps = 1;
    if (evol_steps < 1) evol_steps = 1;

    double *x     = (double*)chaos_calloc(dim, sizeof(double));
    double *x_pert = (double*)chaos_calloc(dim, sizeof(double));

    for (i = 0; i < dim; i++) x[i] = x0[i];

    /* 生成随机扰动 */
    unsigned int seed = 54321;
    double norm = 0.0;
    for (i = 0; i < dim; i++) {
        seed = (1103515245 * seed + 12345) & 0x7fffffff;
        x_pert[i] = (double)seed / 0x7fffffff - 0.5;
        norm += x_pert[i] * x_pert[i];
    }
    norm = sqrt(norm);
    if (norm > 1e-30) {
        for (i = 0; i < dim; i++)
            x_pert[i] = x[i] + x_pert[i] * d0 / norm;
    }

    double log_sum = 0.0;
    int count = 0;
    double t = 0.0;
    double *x_next = (double*)chaos_calloc(dim, sizeof(double));
    double *xp_next = (double*)chaos_calloc(dim, sizeof(double));
    int step;

    for (step = 0; step < n_steps; step++) {
        rk4_step(x_next,  f, x,      t, dt, dim, p);
        rk4_step(xp_next, f, x_pert, t, dt, dim, p);

        for (i = 0; i < dim; i++) {
            x[i]      = x_next[i];
            x_pert[i] = xp_next[i];
        }

        if ((step + 1) % evol_steps == 0) {
            double d = 0.0;
            for (i = 0; i < dim; i++) {
                double diff = x_pert[i] - x[i];
                d += diff * diff;
            }
            d = sqrt(d);

            if (d > 1e-15) {
                log_sum += log(d / d0);
                count++;
                for (i = 0; i < dim; i++)
                    x_pert[i] = x[i] + (x_pert[i] - x[i]) * d0 / d;
            }
        }
        t += dt;
    }

    free(x); free(x_pert); free(x_next); free(xp_next);

    if (count == 0) return 0.0;
    return log_sum / (count * evol_steps * dt);
}

/* ──────────────────────────────────────────────────────────
 * L3: Gram-Schmidt正交化
 * ────────────────────────────────────────────────────────── */

/**
 * Gram-Schmidt正交归一化(列向量)
 *
 * 算法:
 *   for j = 0..dim-1:
 *     for k = 0..j-1:  v_j = v_j - (v_j·v_k)*v_k
 *     r_jj = |v_j|
 *     v_j = v_j / r_jj
 *
 * @param vectors   dim×dim矩阵, 列优先存储(vectors[j*dim+i] = 第j列第i行)
 * @param log_norms [out] 每列的log(|v_j|), 用于Lyapunov累加
 *
 * 复杂度: O(dim^3)
 */
void gram_schmidt(double *vectors, double *log_norms, int dim)
{
    int i, j, k;

    for (j = 0; j < dim; j++) {
        /* 对前j-1列正交化 */
        for (k = 0; k < j; k++) {
            double dot = 0.0;
            for (i = 0; i < dim; i++)
                dot += vectors[j * dim + i] * vectors[k * dim + i];
            for (i = 0; i < dim; i++)
                vectors[j * dim + i] -= dot * vectors[k * dim + i];
        }

        /* 归一化 */
        double norm = 0.0;
        for (i = 0; i < dim; i++)
            norm += vectors[j * dim + i] * vectors[j * dim + i];
        norm = sqrt(norm);

        if (norm > 1e-30) {
            log_norms[j] = log(norm);
            for (i = 0; i < dim; i++)
                vectors[j * dim + i] /= norm;
        } else {
            log_norms[j] = -1e10;  /* 数值零 */
        }
    }
}

/* ──────────────────────────────────────────────────────────
 * L5: Benettin全Lyapunov谱算法
 * ────────────────────────────────────────────────────────── */

/**
 * Benettin全谱算法
 *
 * 原理: 同时追踪dim个正交扰动方向(切空间的一组基)
 *       周期性Gram-Schmidt重正化, 累积对角线元素的log值
 *
 * 步骤:
 *   1. 在参考轨道上积分
 *   2. 用数值Jacobian积分每个扰动基向量
 *   3. 每隔n_renorm步: Gram-Schmidt正交化W
 *   4. lambda_i = sum(log|W的第i个对角线|) / (count * dt * n_renorm)
 *
 * 收敛: O(1/sqrt(t)), 需长时间积分
 *
 * @param spect        [out] 预分配结构
 * @param n_transient  丢弃前n_transient步(瞬态)
 * @param n_renorm     重正化间隔
 * @return 0=成功, -1=失败
 *
 * 复杂度: O(n_steps * dim^3)
 */
int lyapunov_spectrum_benettin(LyapunovSpectrum *spect,
                                flow_rhs_t f, const double *x0,
                                int dim, const double *p,
                                double t_end, double dt,
                                double d0, int n_transient,
                                int n_renorm)
{
    if (!spect || dim < 1 || dim > CHAOS_MAX_DIM) return -1;

    int i, j, n_steps = (int)(t_end / dt);
    if (n_steps < 1) n_steps = 1;
    if (n_renorm < 1) n_renorm = 1;

    /* 参考轨道 */
    double *x = (double*)chaos_calloc(dim, sizeof(double));
    for (i = 0; i < dim; i++) x[i] = x0[i];

    /* 扰动基 — 初始化为单位矩阵(列向量存储) */
    double *W = (double*)chaos_calloc(dim * dim, sizeof(double));
    for (i = 0; i < dim; i++) W[i * dim + i] = 1.0;

    /* Lyapunov累加器 */
    double *sums   = (double*)chaos_calloc(dim, sizeof(double));
    double *lnorms = (double*)chaos_calloc(dim, sizeof(double));
    int count = 0;

    double *x_new = (double*)chaos_calloc(dim, sizeof(double));
    double *x_p   = (double*)chaos_calloc(dim, sizeof(double));
    double *x_m   = (double*)chaos_calloc(dim, sizeof(double));
    double *fp    = (double*)chaos_calloc(dim, sizeof(double));
    double *fm    = (double*)chaos_calloc(dim, sizeof(double));
    double t = 0.0;
    int step;

    for (step = 0; step < n_steps; step++) {
        /* 积分参考轨道 */
        rk4_step(x_new, f, x, t, dt, dim, p);

        /* 计算Jacobian (固定eps) */
        double eps = d0;
        double *Jmat = (double*)chaos_calloc(dim * dim, sizeof(double));
        int k_col;
        for (k_col = 0; k_col < dim; k_col++) {
            for (i = 0; i < dim; i++) {
                x_p[i] = x[i]; x_m[i] = x[i];
            }
            x_p[k_col] += eps;
            x_m[k_col] -= eps;
            f(fp, x_p, t, p);
            f(fm, x_m, t, p);
            for (i = 0; i < dim; i++)
                Jmat[i * dim + k_col] = (fp[i] - fm[i]) / (2.0 * eps);
        }

        /* 用J推进每个扰动基: W_j = J * W_j */
        double *W_new = (double*)chaos_calloc(dim * dim, sizeof(double));
        for (j = 0; j < dim; j++) {
            int k_row;
            for (k_row = 0; k_row < dim; k_row++) {
                double sum = 0.0;
                for (k_col = 0; k_col < dim; k_col++)
                    sum += Jmat[k_row * dim + k_col] * W[j * dim + k_col];
                W_new[j * dim + k_row] = sum;
            }
        }
        for (i = 0; i < dim * dim; i++) W[i] = W_new[i];
        free(W_new);
        free(Jmat);

        /* 周期性重正化 */
        if ((step + 1) % n_renorm == 0) {
            gram_schmidt(W, lnorms, dim);
            if (step >= n_transient) {
                for (i = 0; i < dim; i++)
                    sums[i] += lnorms[i];
                count++;
            }
        }

        for (i = 0; i < dim; i++) x[i] = x_new[i];
        t += dt;
    }

    /* 计算Lyapunov指数 */
    if (count > 0) {
        for (i = 0; i < dim; i++)
            sums[i] /= (count * n_renorm * dt);
    }

    /* 降序排序 */
    for (i = 0; i < dim - 1; i++) {
        for (j = i + 1; j < dim; j++) {
            if (sums[j] > sums[i]) {
                double tmp = sums[i]; sums[i] = sums[j]; sums[j] = tmp;
            }
        }
    }

    /* 填充输出 */
    spect->n_exponents = dim;
    spect->exponents   = (double*)chaos_calloc(dim, sizeof(double));
    for (i = 0; i < dim; i++)
        spect->exponents[i] = sums[i];

    spect->sum_positive  = 0.0;
    for (i = 0; i < dim && sums[i] > 0.0; i++)
        spect->sum_positive += sums[i];

    spect->kaplan_yorke = kaplan_yorke_dimension(sums, dim);

    free(x); free(W); free(sums); free(lnorms);
    free(x_new); free(x_p); free(x_m); free(fp); free(fm);
    return 0;
}

/* ──────────────────────────────────────────────────────────
 * L2: Kaplan-Yorke维数
 * ────────────────────────────────────────────────────────── */

/**
 * Kaplan-Yorke (Lyapunov)维数
 *
 * D_KY = k + (lambda_1 + ... + lambda_k) / |lambda_{k+1}|
 *
 * 其中k是满足 sum_{i=1}^k lambda_i >= 0 的最大整数
 *
 * 性质:
 *   - D_KY >= 拓扑维数
 *   - 对典型奇怪吸引子, D_KY ≈ 信息维数 ≈ 关联维数
 *
 * Lorenz: D_KY ≈ 2.06
 * Rossler: D_KY ≈ 2.01
 *
 * 参考: Frederickson, Kaplan, Yorke, Yorke (1983)
 *       Kaplan & Yorke (1979)
 */
double kaplan_yorke_dimension(const double *lyaps, int n)
{
    if (n < 1 || !lyaps) return 0.0;

    double sum = 0.0;
    int k;

    for (k = 0; k < n; k++) {
        sum += lyaps[k];
        if (sum < 0.0) break;
    }

    if (k == 0) return 0.0;                 /* 所有指数为负 -> 不动点 */
    if (k == n) return (double)n;           /* 所有非负 -> 拟周期环面 */

    return k + (sum - lyaps[k]) / fabs(lyaps[k]);
}

/* ──────────────────────────────────────────────────────────
 * L6: 离散映射的Lyapunov指数
 * ────────────────────────────────────────────────────────── */

/**
 * Logistic映射Lyapunov指数
 *
 * 解析: f'(x) = r*(1-2x)
 * lambda = (1/n) * sum_{k=1}^n log|r*(1-2*x_k)|
 *
 * 理论值:
 *   r=2.0:  lambda=0 (超稳定不动点)
 *   r=3.5:  lambda≈0.318
 *   r=3.569946: lambda=0 (Feigenbaum点)
 *   r=4.0:  lambda=log(2)≈0.693 (除去测度零的非典型轨道)
 *
 * 复杂度: O(n_iter + n_transient)
 */
double logistic_lyapunov(double r, int n_iter, int n_transient)
{
    double x = 0.37;
    int i;

    for (i = 0; i < n_transient; i++)
        x = logistic_map(x, r);

    double lam_sum = 0.0;
    for (i = 0; i < n_iter; i++) {
        x = logistic_map(x, r);
        double deriv = fabs(r * (1.0 - 2.0 * x));
        lam_sum += (deriv < 1e-300) ? -30.0 : log(deriv);
    }

    return lam_sum / n_iter;
}

/**
 * Logistic映射Lyapunov扫描
 *
 * @param r_vals  参数值数组
 * @param n_r     参数值个数
 * @param lyaps   [out] Lyapunov指数(长度n_r)
 */
void logistic_lyapunov_scan(const double *r_vals, int n_r,
                             double *lyaps, int n_iter)
{
    int i;
    int n_trans = (n_iter > 1000) ? 500 : n_iter / 10;
    for (i = 0; i < n_r; i++)
        lyaps[i] = logistic_lyapunov(r_vals[i], n_iter, n_trans);
}

/**
 * 通用1D映射Lyapunov指数(数值导数版)
 *
 * 用数值差分估计f'(x), 适用于无解析导数的映射
 *
 * @param dx  差分步长
 */
double map1d_lyapunov(map1d_t f, double x0, double param,
                       int n_iter, int n_transient)
{
    double x = x0;
    int i;
    double dx = 1e-8;

    for (i = 0; i < n_transient; i++)
        x = f(x, param);

    double lam_sum = 0.0;
    for (i = 0; i < n_iter; i++) {
        double fp = f(x + dx, param);
        double fm = f(x - dx, param);
        double deriv = (fp - fm) / (2.0 * dx);
        double aderiv = fabs(deriv);
        lam_sum += (aderiv < 1e-300) ? -30.0 : log(aderiv);
        x = f(x, param);
    }

    return lam_sum / n_iter;
}

/**
 * Henon映射Lyapunov谱(QR迭代法)
 *
 * Jacobian: J = | -2a*x  1 |
 *              |   b    0 |
 *
 * 每次迭代: W = J * Q, then QR分解W = Q*R
 * lambda_i = (1/n) * sum log|R_{ii}|
 *
 * @param lam1, lam2 [out] 两个Lyapunov指数
 *
 * 理论值(a=1.4, b=0.3): lambda1≈0.419, lambda2≈-1.623
 * D_KY ≈ 1 + 0.419/1.623 ≈ 1.258
 */
void henon_lyapunov_spectrum(double *lam1, double *lam2,
                              double a, double b,
                              int n_iter, int n_transient)
{
    double x = 0.5, y = 0.5;
    int i;

    for (i = 0; i < n_transient; i++) {
        double xn = 1.0 - a * x * x + y;
        double yn = b * x;
        x = xn; y = yn;
    }

    /* 初始化正交基 */
    double Q[4] = {1.0, 0.0, 0.0, 1.0};  /* 列向量: Q[c*2+r] */
    double sum1 = 0.0, sum2 = 0.0;

    for (i = 0; i < n_iter; i++) {
        /* J = | -2a*x  1 | */
        /*     |   b    0 | */
        double Jxx = -2.0 * a * x;
        double Jxy = 1.0;
        double Jyx = b;
        double Jyy = 0.0;

        /* W = J * Q (W列向量) */
        double Wxx = Jxx * Q[0] + Jxy * Q[2];
        double Wyx = Jyx * Q[0] + Jyy * Q[2];
        double Wxy = Jxx * Q[1] + Jxy * Q[3];
        double Wyy = Jyx * Q[1] + Jyy * Q[3];

        /* Gram-Schmidt: Q1 = W1/|W1|, Q2 = (W2-(W2·Q1)Q1)/|...| */
        double r11 = sqrt(Wxx * Wxx + Wyx * Wyx);
        if (r11 > 1e-30) {
            Q[0] = Wxx / r11;
            Q[2] = Wyx / r11;
        }
        sum1 += (r11 > 1e-30) ? log(r11) : -30.0;

        double dot = Wxy * Q[0] + Wyy * Q[2];
        double r22_x = Wxy - dot * Q[0];
        double r22_y = Wyy - dot * Q[2];
        double r22 = sqrt(r22_x * r22_x + r22_y * r22_y);
        if (r22 > 1e-30) {
            Q[1] = r22_x / r22;
            Q[3] = r22_y / r22;
        }
        sum2 += (r22 > 1e-30) ? log(r22) : -30.0;

        /* 迭代映射 */
        double xn = 1.0 - a * x * x + y;
        double yn = b * x;
        x = xn; y = yn;
    }

    *lam1 = sum1 / n_iter;
    *lam2 = sum2 / n_iter;
}

/**
 * 条件Lyapunov指数(同步稳定性判据)
 *
 * 横向Lyapunov指数: lambda_perp < 0 => 同步稳定
 *
 * 方法:
 *   - 驱动系统x用其自身方程积分
 *   - 响应系统y的RHS依赖驱动信号x_drive
 *   - 在两个响应系统副本间追踪扰动增长率
 *
 * @param f_drive     驱动系统的RHS
 * @param f_response  响应系统的RHS (含耦合项, 用x_drive)
 * @param x0          驱动系统初值(dim_x维)
 * @param y0          响应系统初值(dim_y维)
 * @param d0          初始横向扰动
 *
 * 参考: Pecora & Carroll (1990) PRL 64:821
 */
double conditional_lyapunov(flow_rhs_t f_drive, flow_rhs_t f_response,
                             const double *x0, const double *y0,
                             int dim_x, int dim_y,
                             const double *p, double t_end,
                             double dt, double d0)
{
    int i, n_steps = (int)(t_end / dt);
    if (n_steps < 1) n_steps = 1;

    double *x_drive = (double*)chaos_calloc(dim_x, sizeof(double));
    double *x_next  = (double*)chaos_calloc(dim_x, sizeof(double));
    double *yA = (double*)chaos_calloc(dim_y, sizeof(double));
    double *yB = (double*)chaos_calloc(dim_y, sizeof(double));
    double *yAn = (double*)chaos_calloc(dim_y, sizeof(double));
    double *yBn = (double*)chaos_calloc(dim_y, sizeof(double));

    for (i = 0; i < dim_x; i++) x_drive[i] = x0[i];
    for (i = 0; i < dim_y; i++) yA[i] = y0[i];

    /* 扰动yB = yA + 随机方向 */
    unsigned int seed = 99999;
    double norm = 0.0;
    for (i = 0; i < dim_y; i++) {
        seed = (1103515245 * seed + 12345) & 0x7fffffff;
        double r = (double)seed / 0x7fffffff - 0.5;
        yB[i] = yA[i] + d0 * r;
        norm += r * r;
    }
    if (norm > 1e-30) {
        norm = sqrt(norm);
        for (i = 0; i < dim_y; i++)
            yB[i] = yA[i] + (yB[i] - yA[i]) * d0 / norm;
    }

    double log_sum = 0.0;
    int count = 0;
    double t = 0.0;
    int step;

    for (step = 0; step < n_steps; step++) {
        /* 积分驱动系统 */
        rk4_step(x_next, f_drive, x_drive, t, dt, dim_x, p);
        for (i = 0; i < dim_x; i++) x_drive[i] = x_next[i];

        /* 将驱动信号合并到响应系统的参数中 */
        double p_combined[CHAOS_MAX_DIM + 4];
        int k;
        for (k = 0; k < 4; k++) p_combined[k] = p[k];
        for (k = 0; k < dim_x && k < CHAOS_MAX_DIM; k++)
            p_combined[4 + k] = x_drive[k];

        /* 积分响应系统 */
        rk4_step(yAn, f_response, yA, t, dt, dim_y, p_combined);
        rk4_step(yBn, f_response, yB, t, dt, dim_y, p_combined);

        double d = 0.0;
        for (i = 0; i < dim_y; i++) {
            double diff = yBn[i] - yAn[i];
            d += diff * diff;
        }
        d = sqrt(d);

        if (d > 1e-15) {
            log_sum += log(d / d0);
            count++;
            for (i = 0; i < dim_y; i++)
                yBn[i] = yAn[i] + (yBn[i] - yAn[i]) * d0 / d;
        }

        for (i = 0; i < dim_y; i++) yA[i] = yAn[i];
        t += dt;
    }

    free(x_drive); free(x_next);
    free(yA); free(yB); free(yAn); free(yBn);

    return (count > 0) ? log_sum / (count * dt) : 0.0;
}
