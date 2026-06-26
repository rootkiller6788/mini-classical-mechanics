/**
 * chaos_fractals.c — 分形几何实现
 *
 * 涵盖:
 *   L1: Mandelbrot集/Julia集定义, Newton分形
 *   L6: 经典分形(Sierpinski三角, Koch雪花)
 *   L5: 盒计数维数, 关联维数
 *
 * 参考:
 *   Mandelbrot "The Fractal Geometry of Nature" (1983)
 *   Peitgen & Richter "The Beauty of Fractals" (1986)
 *   Falconer "Fractal Geometry" (2013)
 */

#include "chaos_fractals.h"
#include <stdio.h>
#include <math.h>

/* ============================================================
 * L2: 逃逸时间算法 — Mandelbrot集和Julia集
 * ============================================================ */

/**
 * Mandelbrot集迭代
 *
 * 动力学: z_0 = 0, z_{n+1} = z_n^2 + c
 *
 * 定理: 若|z_n| > 2 对某个n成立, 则轨道必逃逸到无穷
 * 证明: 若|z| > 2, 则|z^2 + c| >= |z|^2 - |c| > 2|z| - |z| = |z|
 *       (当|z|>2且|c|<|z|时), 所以|z|单调增 → 发散
 *
 * 基本性质:
 *   - M集是连通的(Douady & Hubbard, 1982)
 *   - M集的边界是分形, 维数=2(Shishikura 1998)
 *   - M集在实轴上[-2, 0.25], 虚轴上[-1.2, 1.2]
 *
 * 参数空间(c平面)上的分形
 *
 * @param cx, cy    复参数c
 * @param max_iter  最大迭代次数
 * @return 逃逸迭代数(max_iter=在M集内)
 *
 * 复杂度: O(max_iter) 对非M集点, 通常<50
 */
int mandelbrot_iter(double cx, double cy, int max_iter)
{
    double zx = 0.0, zy = 0.0;
    int n;

    for (n = 0; n < max_iter; n++) {
        double zx2 = zx * zx;
        double zy2 = zy * zy;
        if (zx2 + zy2 > 4.0) break;  /* 逃逸条件: |z|>2 => |z|^2>4 */
        zy = 2.0 * zx * zy + cy;
        zx = zx2 - zy2 + cx;
    }
    return n;
}

/**
 * Julia集迭代
 *
 * 动力学: z_{n+1} = z_n^2 + c (固定c, 变z_0)
 *
 * Julia集J_c: 排斥不动点集的闭包
 * Fatou集F_c: C \ J_c (稳定集)
 *
 * 二分法: 对给定c, Julia集要么连通(当c在M集内)要么是Cantor集
 *
 * 知识: Julia和Fatou集基于Montel正规族理论
 * 参考: Milnor "Dynamics in One Complex Variable" (2006)
 *
 * 复杂度: O(max_iter)
 */
int julia_iter(double zx, double zy, double cx, double cy, int max_iter)
{
    int n;
    for (n = 0; n < max_iter; n++) {
        double zx2 = zx * zx;
        double zy2 = zy * zy;
        if (zx2 + zy2 > 4.0) break;
        zy = 2.0 * zx * zy + cy;
        zx = zx2 - zy2 + cx;
    }
    return n;
}

/**
 * 距离估计法 — 平滑着色Mandelbrot
 *
 * 公式: d(c) = 0.5 * |z_n| * log|z_n| / |z'_n|
 * 其中z'_n = dz_n/dc 是导数
 *
 * 迭代: z'_{n+1} = 2*z_n*z'_n + 1 (对c的导数)
 *
 * 知识点: 势函数与等势线, 距离估计的外边界定理
 * 参考: Peitgen & Richter §4.5
 */
double mandelbrot_distance(double cx, double cy, int max_iter)
{
    double zx = 0.0, zy = 0.0;
    double dzx = 1.0, dzy = 0.0;  /* z' 初始: d/dc(z_0) = d/dc(0) = 0, 但第一步后=1 */
    int n;
    int zprime_inited = 0;

    for (n = 0; n < max_iter; n++) {
        double zx2 = zx * zx, zy2 = zy * zy;
        double mod2 = zx2 + zy2;
        if (mod2 > 4.0) {
            double mod = sqrt(mod2);
            double dmod = sqrt(dzx * dzx + dzy * dzy);
            if (dmod > 1e-300)
                return 0.5 * mod * log(mod) / dmod;
            return 0.0;
        }

        /* 更新导数: dz/dc */
        if (!zprime_inited) {
            dzx = 1.0; dzy = 0.0;
            zprime_inited = 1;
        }
        double dzx_new = 2.0 * (zx * dzx - zy * dzy) + 1.0;
        double dzy_new = 2.0 * (zx * dzy + zy * dzx);
        dzx = dzx_new;
        dzy = dzy_new;

        /* 更新z */
        zy = 2.0 * zx * zy + cy;
        zx = zx2 - zy2 + cx;
    }
    return 0.0;  /* 在M集内 */
}

/* ============================================================
 * L6: 分形图像生成
 * ============================================================ */

/**
 * Mandelbrot集图像
 *
 * 知识点: 图像网格采样, 复平面映射
 * 复杂度: O(width * height * max_iter)
 */
FractalImage* mandelbrot_set_image(double x_min, double x_max,
                                    double y_min, double y_max,
                                    int width, int height, int max_iter)
{
    FractalImage *img = (FractalImage*)chaos_calloc(1, sizeof(FractalImage));
    img->width  = width;
    img->height = height;
    img->x_min  = x_min;
    img->x_max  = x_max;
    img->y_min  = y_min;
    img->y_max  = y_max;
    img->max_iter = max_iter;
    img->data   = (int*)chaos_calloc((size_t)width * height, sizeof(int));

    double dx = (x_max - x_min) / (width - 1);
    double dy = (y_max - y_min) / (height - 1);
    int px, py;

    for (py = 0; py < height; py++) {
        double cy = y_min + py * dy;
        for (px = 0; px < width; px++) {
            double cx = x_min + px * dx;
            img->data[py * width + px] = mandelbrot_iter(cx, cy, max_iter);
        }
    }
    return img;
}

/**
 * Julia集图像
 */
FractalImage* julia_set_image(double x_min, double x_max,
                               double y_min, double y_max,
                               int width, int height,
                               double cx, double cy, int max_iter)
{
    FractalImage *img = (FractalImage*)chaos_calloc(1, sizeof(FractalImage));
    img->width  = width;
    img->height = height;
    img->x_min  = x_min;
    img->x_max  = x_max;
    img->y_min  = y_min;
    img->y_max  = y_max;
    img->max_iter = max_iter;
    img->data   = (int*)chaos_calloc((size_t)width * height, sizeof(int));

    double dx = (x_max - x_min) / (width - 1);
    double dy = (y_max - y_min) / (height - 1);
    int px, py;

    for (py = 0; py < height; py++) {
        double zy = y_min + py * dy;
        for (px = 0; px < width; px++) {
            double zx = x_min + px * dx;
            img->data[py * width + px] = julia_iter(zx, zy, cx, cy, max_iter);
        }
    }
    return img;
}

/* ============================================================
 * L1: Newton分形
 * ============================================================ */

/**
 * Newton分形: z^3 - 1 = 0
 *
 * 迭代: z_{n+1} = z_n - (z_n^3 - 1)/(3*z_n^2) = (2*z_n + 1/z_n^2)/3
 *
 * 三个根: 1, exp(2*pi*i/3), exp(4*pi*i/3)
 * 即: (1,0), (-1/2, sqrt(3)/2), (-1/2, -sqrt(3)/2)
 *
 * 吸引盆具有Wada性质: 每个边界点同时属于所有三个盆的边界
 *
 * 知识点: 复分析中的Newton法收敛性, 分形盆边界
 * 参考: Cayley问题(1879), Hubbard, Schleicher, Sutherland (2001)
 *
 * @return 收敛到的根索引(1,2,3)或0(未收敛)
 */
int newton_fractal(double zx, double zy, int max_iter, double tol)
{
    double roots[6] = {
        1.0, 0.0,
        -0.5,  0.8660254037844386,
        -0.5, -0.8660254037844386
    };
    int n, k;

    for (n = 0; n < max_iter; n++) {
        double zx2 = zx * zx, zy2 = zy * zy;
        double denom = 3.0 * (zx2 + zy2);
        if (fabs(denom) < 1e-300) break;

        /* z^3 - 1 */
        double num_x = zx * (zx2 - 3.0 * zy2) - 1.0;
        double num_y = zy * (3.0 * zx2 - zy2);

        /* z - (z^3-1)/(3z^2) = (2z^3+1)/(3z^2) */
        zx = zx - num_x / denom;
        zy = zy - num_y / denom;

        /* 检查收敛到哪个根 */
        for (k = 0; k < 3; k++) {
            double dx = zx - roots[2 * k];
            double dy = zy - roots[2 * k + 1];
            if (sqrt(dx * dx + dy * dy) < tol)
                return k + 1;
        }
    }
    return 0;
}

/**
 * Newton分形图像
 */
FractalImage* newton_fractal_image(double x_min, double x_max,
                                    double y_min, double y_max,
                                    int width, int height,
                                    int max_iter, double tol)
{
    FractalImage *img = (FractalImage*)chaos_calloc(1, sizeof(FractalImage));
    img->width  = width;
    img->height = height;
    img->x_min  = x_min;
    img->x_max  = x_max;
    img->y_min  = y_min;
    img->y_max  = y_max;
    img->max_iter = max_iter;
    img->data   = (int*)chaos_calloc((size_t)width * height, sizeof(int));

    double dx = (x_max - x_min) / (width - 1);
    double dy = (y_max - y_min) / (height - 1);
    int px, py;

    for (py = 0; py < height; py++) {
        double zy = y_min + py * dy;
        for (px = 0; px < width; px++) {
            double zx = x_min + px * dx;
            img->data[py * width + px] = newton_fractal(zx, zy, max_iter, tol);
        }
    }
    return img;
}

/* ============================================================
 * L6: 经典分形 — 确定性迭代
 * ============================================================ */

/**
 * Sierpinski三角(混沌游戏法)
 *
 * 算法:
 *   顶点: (0,0), (1,0), (0.5, sqrt(3)/2)
 *   1. 从任意点开始
 *   2. 随机选一个顶点
 *   3. 移动到当前点与选中顶点的中点
 *   4. 记录位置
 *
 * 结论: 点的极限集合 = Sierpinski三角
 *
 * 原理: 三个收缩映射(比例1/2)构成的IFS的不变集
 *
 * @param points [out] 长度为2*n的数组[x1,y1,x2,y2,...]
 * @param n     点数
 *
 * 参考: Barnsley "Fractals Everywhere" (1993)
 * Cambridge Part III: IFS理论
 */
void sierpinski_triangle(double *points, int n)
{
    double vx[3] = {0.0, 1.0, 0.5};
    double vy[3] = {0.0, 0.0, 0.8660254037844386};
    double px = 0.5, py = 0.25;
    int i;
    unsigned int seed = 1234;

    /* 简单LCG随机 */
    #define CHAOS_RAND() (seed = (1103515245 * seed + 12345) & 0x7fffffff, seed)

    for (i = 0; i < n; i++) {
        int v = CHAOS_RAND() % 3;
        px = (px + vx[v]) * 0.5;
        py = (py + vy[v]) * 0.5;
        points[2 * i]     = px;
        points[2 * i + 1] = py;
    }
}

/**
 * Koch雪花曲线
 *
 * L-system: 公理=F++F++F, F→F-F++F-F, θ=60°
 *
 * 每代替换: 每条线段 → 4条线段(原长的1/3)
 * n=0: 等边三角形
 * n=1: 六角星
 * n->inf: Koch雪花, 周长趋于无穷, 面积收敛
 *
 * 维数: D = log(4)/log(3) ≈ 1.26186
 *
 * @param seg_x    [out] x坐标数组
 * @param seg_y    [out] y坐标数组
 * @param max_segs 预分配大小
 * @param n_iter   迭代次数
 * @return 实际线段数(也是点数)
 */
int koch_snowflake(double *seg_x, double *seg_y, int max_segs, int n_iter)
{
    if (n_iter < 0 || max_segs < 4) return 0;

    /* 初始等边三角形 */
    double tri_x[4] = {0.0, 1.0, 0.5, 0.0};
    double tri_y[4] = {0.0, 0.0, 0.8660254037844386, 0.0};
    double *cur_x = (double*)chaos_calloc(max_segs, sizeof(double));
    double *cur_y = (double*)chaos_calloc(max_segs, sizeof(double));
    double *next_x = (double*)chaos_calloc(max_segs, sizeof(double));
    double *next_y = (double*)chaos_calloc(max_segs, sizeof(double));

    int n_cur = 4;
    int i, iter;

    for (i = 0; i < 4; i++) {
        cur_x[i] = tri_x[i];
        cur_y[i] = tri_y[i];
    }

    for (iter = 0; iter < n_iter; iter++) {
        int n_next = 0;
        for (i = 0; i < n_cur - 1; i++) {
            if (n_next + 5 > max_segs) break;
            double x1 = cur_x[i],   y1 = cur_y[i];
            double x2 = cur_x[i+1], y2 = cur_y[i+1];
            double dx = (x2 - x1) / 3.0;
            double dy = (y2 - y1) / 3.0;

            /* A = p1 + (p2-p1)/3 */
            double ax = x1 + dx, ay = y1 + dy;
            /* B = p1 + (p2-p1)/2 + sqrt(3)/6 * rotated90(p2-p1) */
            double cos60 = 0.5, sin60 = 0.8660254037844386;
            double rx = dx * cos60 - dy * sin60;
            double ry = dx * sin60 + dy * cos60;
            double bx = x1 + dx + rx, by = y1 + dy + ry;
            /* C = p1 + 2*(p2-p1)/3 */
            double cx = x1 + 2.0 * dx, cy = y1 + 2.0 * dy;

            next_x[n_next] = x1; next_y[n_next] = y1; n_next++;
            next_x[n_next] = ax; next_y[n_next] = ay; n_next++;
            next_x[n_next] = bx; next_y[n_next] = by; n_next++;
            next_x[n_next] = cx; next_y[n_next] = cy; n_next++;
        }
        next_x[n_next] = cur_x[n_cur-1];
        next_y[n_next] = cur_y[n_cur-1];
        n_next++;

        /* 复制回cur */
        for (i = 0; i < n_next && i < max_segs; i++) {
            cur_x[i] = next_x[i];
            cur_y[i] = next_y[i];
        }
        n_cur = (n_next < max_segs) ? n_next : max_segs;
    }

    for (i = 0; i < n_cur; i++) {
        seg_x[i] = cur_x[i];
        seg_y[i] = cur_y[i];
    }

    free(cur_x); free(cur_y); free(next_x); free(next_y);
    return n_cur;
}

/* ============================================================
 * L5: 分形维数 — 盒计数法
 * ============================================================ */

/**
 * 盒计数维数
 *
 * 算法:
 *   1. 对每个eps, 将空间划分成eps×eps的盒子
 *   2. 统计包含至少一个点的盒子数N(eps)
 *   3. 线性回归: log(N) ~ D*log(1/eps)
 *
 * 盒计数维数 = 容量维数 = Kolmogorov维数
 *
 * @param points   点数组[x1,y1,x2,y2,...]
 * @param n_points 点数
 * @param epsilons 尺度序列
 * @param n_eps    尺度个数
 * @return 维数估计(最小二乘斜率)
 *
 * 复杂度: O(n_eps * n_points)
 */
double box_counting_dimension(const double *points, int n_points,
                               const double *epsilons, int n_eps)
{
    if (n_points < 2 || n_eps < 2) return 0.0;

    double *log_N = (double*)chaos_calloc(n_eps, sizeof(double));
    double *log_inv = (double*)chaos_calloc(n_eps, sizeof(double));
    int ie;

    for (ie = 0; ie < n_eps; ie++) {
        double eps = epsilons[ie];
        if (eps <= 0.0) { log_N[ie] = 0.0; log_inv[ie] = 0.0; continue; }

        /* 用哈希集合统计非空盒子 */
        /* 简化: 使用int坐标映射到(ix,iy)对, 用O(n^2)去重 */
        int n_boxes = 0;
        int *box_ids = (int*)chaos_calloc(n_points, sizeof(int));
        int i, j;
        for (i = 0; i < n_points; i++) {
            int ix = (int)floor(points[2 * i] / eps);
            int iy = (int)floor(points[2 * i + 1] / eps);
            int box_id = ix * 1000000 + iy;  /* 简化: 假设范围不太大 */

            /* 检查是否已存在 */
            int found = 0;
            for (j = 0; j < n_boxes; j++) {
                if (box_ids[j] == box_id) { found = 1; break; }
            }
            if (!found) box_ids[n_boxes++] = box_id;
        }

        log_N[ie] = log((double)n_boxes);
        log_inv[ie] = log(1.0 / eps);
        free(box_ids);
    }

    /* 最小二乘线性回归: log_N = D * log_inv + b */
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
    for (ie = 0; ie < n_eps; ie++) {
        sum_x  += log_inv[ie];
        sum_y  += log_N[ie];
        sum_xy += log_inv[ie] * log_N[ie];
        sum_xx += log_inv[ie] * log_inv[ie];
    }

    double slope = 0.0;
    double denom = n_eps * sum_xx - sum_x * sum_x;
    if (fabs(denom) > 1e-15)
        slope = (n_eps * sum_xy - sum_x * sum_y) / denom;

    free(log_N); free(log_inv);
    return slope;
}

/**
 * 关联维数(Grassberger-Procaccia)
 *
 * 关联积分: C(eps) = (2/(N*(N-1))) * sum_{i<j} Theta(eps - |x_i - x_j|)
 * 对于小eps: C(eps) ∝ eps^{D_2}
 *
 * D_2 = lim_{eps->0} log(C(eps)) / log(eps)
 *
 * @param points  点集(按行: points[i*dim+k])
 * @param n       点数
 * @param dim     每点维度
 * @param epsilons eps值数组
 * @param n_eps   eps值个数
 * @return D_2估计
 *
 * 复杂度: O(n_eps * n^2 * dim) — 双循环
 * 参考: Grassberger & Procaccia (1983) PRL 50:346
 */
double correlation_dimension_gp(const double *points, int n_points,
                                 int dim, const double *epsilons,
                                 int n_eps)
{
    if (n_points < 2 || n_eps < 2 || dim < 1) return 0.0;

    double *log_C = (double*)chaos_calloc(n_eps, sizeof(double));
    double *log_eps = (double*)chaos_calloc(n_eps, sizeof(double));
    int ie, i, j;

    for (ie = 0; ie < n_eps; ie++) {
        double eps = epsilons[ie];
        if (eps <= 0.0) { log_C[ie] = 0.0; log_eps[ie] = 0.0; continue; }

        long long count = 0;
        for (i = 0; i < n_points - 1; i++) {
            for (j = i + 1; j < n_points; j++) {
                double dist2 = 0.0;
                int k;
                for (k = 0; k < dim; k++) {
                    double diff = points[i * dim + k] - points[j * dim + k];
                    dist2 += diff * diff;
                }
                if (dist2 < eps * eps) count++;
            }
        }

        double C_val = (2.0 * count) / ((double)n_points * (n_points - 1));
        if (C_val < 1e-300) C_val = 1e-300;
        log_C[ie]   = log(C_val);
        log_eps[ie] = log(eps);
    }

    /* 线性回归 */
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
    for (ie = 0; ie < n_eps; ie++) {
        sum_x  += log_eps[ie];
        sum_y  += log_C[ie];
        sum_xy += log_eps[ie] * log_C[ie];
        sum_xx += log_eps[ie] * log_eps[ie];
    }

    double slope = 0.0;
    double denom = n_eps * sum_xx - sum_x * sum_x;
    if (fabs(denom) > 1e-15)
        slope = (n_eps * sum_xy - sum_x * sum_y) / denom;

    free(log_C); free(log_eps);
    return slope;
}

/* ──────────────────────────────────────────────────────────
 * 理论维数公式
 * ────────────────────────────────────────────────────────── */

/** Koch雪花维数: D=log(4)/log(3) ≈ 1.26186 */
double koch_dimension(void) { return log(4.0) / log(3.0); }

/** Sierpinski三角维数: D=log(3)/log(2) ≈ 1.58496 */
double sierpinski_dimension(void) { return log(3.0) / log(2.0); }

/** Cantor集维数: D=log(2)/log(3) ≈ 0.63093 */
double cantor_dimension(void) { return log(2.0) / log(3.0); }
