/**
 * chaos_analysis.c — 混沌高级分析: 递归分析, 同步, 控制, 统计
 *
 * 涵盖:
 *   L2: 递归图(Eckmann et al. 1987)
 *   L7: 混沌控制(OGY方法), 混沌同步
 *   L8: 复杂网络分析, Kuramoto模型
 *
 * 参考:
 *   Marwan et al. (2007) "Recurrence Plots..."
 *   Ott, Grebogi, Yorke (1990) PRL 64:1196
 *   Pecora & Carroll (1990) PRL 64:821
 */

#include "chaos.h"
#include "chaos_flows.h"
#include "chaos_embedding.h"
#include "chaos_maps.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

/* ============================================================
 * L2: 递归图 (Recurrence Plot)
 * ============================================================ */

/**
 * 递归图
 *
 * 定义: R_{ij} = Theta(epsilon - |x_i - x_j|), i,j = 1..N
 *
 * 其中 Theta 是Heaviside阶跃函数(>0时为1, 否则0)
 *
 * 模式(Eckmann et al. 1987):
 *   - 均匀分布点: 均匀性/平稳性
 *   - 对角线: 确定性/周期性(R_{i+k, j+k} = 1)
 *   - 水平/垂直线: 层流状态/间歇性
 *   - 白色带: 非平稳性/突变
 *   - 棋盘格: 交替振荡
 *
 * 递归率: RR = (1/N^2) * sum_{i,j} R_{ij}
 * 确定性:  DET = sum_{l >= l_min} l * p(l) / sum_{i,j} R_{ij}
 *          其中p(l)是对角线长度l的频率
 *
 * @param points  相空间点集
 * @param n       点数
 * @param dim     每点维度
 * @param epsilon 阈值
 * @return RecurrencePlot (caller free)
 *
 * 复杂度: O(n^2 * dim)
 */
RecurrencePlot* recurrence_plot_create(const double *points,
                                        int n, int dim, double epsilon)
{
    RecurrencePlot *rp = (RecurrencePlot*)chaos_calloc(1, sizeof(RecurrencePlot));
    rp->n        = n;
    rp->epsilon  = epsilon;
    rp->matrix   = (int*)chaos_calloc((size_t)n * n, sizeof(int));

    int i, j, total_ones = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double dist2 = 0.0;
            int k;
            for (k = 0; k < dim; k++) {
                double diff = points[i * dim + k] - points[j * dim + k];
                dist2 += diff * diff;
            }
            if (dist2 < epsilon * epsilon) {
                rp->matrix[i * n + j] = 1;
                total_ones++;
            }
        }
    }

    rp->recurrence_rate = (double)total_ones / ((double)n * n);
    return rp;
}

/**
 * 从轨迹生成递归图
 */
RecurrencePlot* recurrence_plot_from_trajectory(const ChaosTrajectory *traj,
                                                  double epsilon)
{
    if (!traj || !traj->data) return NULL;
    return recurrence_plot_create(traj->data, traj->n_points,
                                   traj->dim, epsilon);
}

/**
 * RQA — 确定性(DET): 对角线结构的比例
 *
 * DET = sum_{l >= l_min} l*P(l) / sum_{l >= l_min} sum_{i,j: l对角线} R_{ij}
 *
 * 其中P(l)是长度为l的对角线的频率
 *
 * @param l_min  最小对角线长度(通常=2)
 * @return DET值
 */
double recurrence_determinism(const RecurrencePlot *rp, int l_min)
{
    if (!rp || l_min < 2) return 0.0;

    int n = rp->n;
    int *diag_hist = (int*)chaos_calloc(n + 1, sizeof(int));
    int *current_diag = (int*)chaos_calloc(n, sizeof(int));
    int i, shift;

    /* 扫描所有偏移的对角线 */
    for (shift = -(n - 1); shift < n; shift++) {
        int len = 0;
        int i_start = (shift < 0) ? -shift : 0;
        int j_start = (shift > 0) ? shift : 0;
        int max_k = n - (i_start > j_start ? i_start : j_start);
        int k;

        for (k = 0; k < max_k; k++) {
            int ii = i_start + k;
            int jj = j_start + k;
            if (rp->matrix[ii * n + jj]) {
                len++;
            } else {
                if (len >= l_min) diag_hist[len]++;
                len = 0;
            }
        }
        if (len >= l_min) diag_hist[len]++;
    }

    double diag_sum = 0.0;
    for (i = l_min; i <= n; i++) {
        diag_sum += (double)i * diag_hist[i];
    }

    int total_rp = 0;
    for (i = 0; i < n * n; i++) {
        if (rp->matrix[i]) total_rp++;
    }

    free(diag_hist); free(current_diag);

    if (total_rp == 0) return 0.0;
    return diag_sum / total_rp;
}

/**
 * 层流性(LAM): 垂直结构的比例
 *
 * LAM = sum_{v >= v_min} v*P(v) / sum_{i,j} R_{ij}
 *
 * 其中P(v)是长度为v的垂直线的频率
 *
 * 高LAM → 层流相(间歇性)
 *
 * @param v_min 最小垂直线长度(通常=2)
 */
double recurrence_laminarity(const RecurrencePlot *rp, int v_min)
{
    if (!rp || v_min < 2) return 0.0;

    int n = rp->n;
    int *vert_hist = (int*)chaos_calloc(n + 1, sizeof(int));
    int col, row;

    for (col = 0; col < n; col++) {
        int len = 0;
        for (row = 0; row < n; row++) {
            if (rp->matrix[row * n + col]) {
                len++;
            } else {
                if (len >= v_min) vert_hist[len]++;
                len = 0;
            }
        }
        if (len >= v_min) vert_hist[len]++;
    }

    double vert_sum = 0.0;
    int i;
    for (i = v_min; i <= n; i++) {
        vert_sum += (double)i * vert_hist[i];
    }

    int total_rp = 0;
    for (i = 0; i < n * n; i++) {
        if (rp->matrix[i]) total_rp++;
    }

    free(vert_hist);

    if (total_rp == 0) return 0.0;
    return vert_sum / total_rp;
}

/**
 * 平均对角线长度L
 *
 * 解释: L与最大Lyapunov指数相关: L ~ exp(lambda_max)
 */
double recurrence_avg_diagonal(const RecurrencePlot *rp, int l_min)
{
    if (!rp || l_min < 2) return 0.0;

    int n = rp->n;
    int total_len = 0, total_count = 0;
    int shift;

    for (shift = -(n - 1); shift < n; shift++) {
        int len = 0;
        int i_start = (shift < 0) ? -shift : 0;
        int j_start = (shift > 0) ? shift : 0;
        int max_k = n - (i_start > j_start ? i_start : j_start);
        int k;

        for (k = 0; k < max_k; k++) {
            if (rp->matrix[(i_start + k) * n + (j_start + k)]) {
                len++;
            } else {
                if (len >= l_min) { total_len += len; total_count++; }
                len = 0;
            }
        }
        if (len >= l_min) { total_len += len; total_count++; }
    }

    if (total_count == 0) return 0.0;
    return (double)total_len / total_count;
}

/* ============================================================
 * L7: 混沌控制 — OGY方法
 * ============================================================ */

/**
 * OGY (Ott-Grebogi-Yorke) 混沌控制
 *
 * 原理: 在奇怪吸引子上存在无穷多个嵌入的不稳定周期轨道(UPO)
 *       通过对可调参数施加微小扰动, 将轨道稳定到期望的UPO
 *
 * 关键洞察: 不动点附近的线性化动力学
 *   x_{n+1} - x* = A*(x_n - x*) + B*(p - p*)
 *
 * 控制律: 选择p_n使x_{n+1}落在稳定流形上
 *   p_n = p* - (lambda_u/(lambda_u - 1)) * f_u · (x_n - x*)
 *
 * 前提条件:
 *   - 吸引子上存在UPO
 *   - 参数p可微调
 *   - 轨道自然回到不动点附近(ergodicity保证)
 *
 * @param map_func  映射函数 f(x, p)
 * @param x_target  目标UPO上的不动点
 * @param p_nominal 标称参数
 * @param K         反馈增益
 * @param n_iter    控制迭代次数
 * @param results   [out] 受控轨道(长度n_iter)
 *
 * 参考: Ott, Grebogi, Yorke (1990) PRL 64:1196
 *       Shinbrot et al. (1993) Nature 363:411
 */
void ogy_control_1d(map1d_t map_func, double x_target, double p_nominal,
                     double K, int n_iter, double *results)
{
    double x = x_target + 0.01;  /* 微小偏移开始 */
    int i;

    for (i = 0; i < n_iter; i++) {
        double delta_x = x - x_target;
        double p = p_nominal - K * delta_x;
        /* 限制参数范围 */
        if (p < p_nominal * 0.5) p = p_nominal * 0.5;
        if (p > p_nominal * 2.0) p = p_nominal * 2.0;

        x = map_func(x, p);
        results[i] = x;
    }
}

/* ============================================================
 * L8: Kuramoto模型 — 耦合同步
 * ============================================================ */

/**
 * Kuramoto模型右侧函数
 *
 * dtheta_i/dt = omega_i + (K/N) * sum_{j} sin(theta_j - theta_i)
 *
 * 序参量: r * exp(i*psi) = (1/N) * sum_{j} exp(i*theta_j)
 * 平均场形式: dtheta_i/dt = omega_i + K*r*sin(psi - theta_i)
 *
 * 临界耦合: K_c = 2 / (pi * g(0))
 * 其中g(omega)是固有频率的分布密度
 *
 * @param theta    相角数组(长度N)
 * @param omega    固有频率数组(长度N)
 * @param K        耦合强度
 * @param dtheta   [out] 相角导数(长度N)
 * @param N        振子个数
 *
 * 参考: Kuramoto (1975, 1984), Strogatz (2000)
 *
 * 知识点: 同步相变, 平均场理论, 临界现象
 * MIT/Cambridge/Oxford 非线性动力学课程
 */
void kuramoto_rhs(const double *theta, const double *omega,
                  double K, double *dtheta, int N)
{
    int i, j;
    for (i = 0; i < N; i++) {
        double coupling = 0.0;
        for (j = 0; j < N; j++) {
            coupling += sin(theta[j] - theta[i]);
        }
        dtheta[i] = omega[i] + (K / N) * coupling;
    }
}

/**
 * Kuramoto序参量
 *
 * r = |(1/N) * sum exp(i*theta_j)|
 * psi = arg((1/N) * sum exp(i*theta_j))
 *
 * r=0: 完全去同步
 * r=1: 完全同步
 *
 * @param r   [out] 序参量幅度
 * @param psi [out] 平均相角
 */
void kuramoto_order_parameter(const double *theta, int N,
                                double *r, double *psi)
{
    double sum_x = 0.0, sum_y = 0.0;
    int i;

    for (i = 0; i < N; i++) {
        sum_x += cos(theta[i]);
        sum_y += sin(theta[i]);
    }

    sum_x /= N;
    sum_y /= N;

    *r   = sqrt(sum_x * sum_x + sum_y * sum_y);
    *psi = atan2(sum_y, sum_x);
}

/**
 * Kuramoto临界耦合强度
 *
 * K_c = 2 / (pi * g(omega_bar))
 *
 * 其中g(omega)是固有频率分布密度在均值处的值
 *
 * 对于Lorentz分布: g(omega) = gamma/(pi*(gamma^2 + (omega-omega0)^2))
 * 则 g(omega0) = 1/(pi*gamma)
 * K_c = 2*gamma
 *
 * 对于高斯分布: g(omega) = 1/(sqrt(2*pi)*sigma)
 * 则近似 K_c ≈ 2*sqrt(2/pi)*sigma
 *
 * @param g0  g(omega_bar)的值
 * @return K_c
 */
double kuramoto_critical_K(double g0)
{
    if (g0 <= 0.0) return 0.0;
    return 2.0 / (M_PI * g0);
}

/* ============================================================
 * L7: 混沌同步 — Pecora-Carroll方法
 * ============================================================ */

/**
 * Lorenz系统间的驱动-响应同步误差
 *
 * 驱动:  Lorenz x (完整3D)
 * 响应:  Lorenz y, 但y[0]被x[0]替换(单向耦合)
 *
 * 响应方程:
 *   dy[0]/dt = sigma*(y[1]-y[0]) + C*(x[0]-y[0])
 *   dy[1]/dt = y[0]*(rho-y[2]) - y[1]
 *   dy[2]/dt = y[0]*y[1] - beta*y[2]
 *
 * 同步判据: 条件Lyapunov指数 < 0
 *
 * 参考: Pecora & Carroll (1990) PRL 64:821
 *
 * @param coupling  耦合强度C
 * @param t_end     积分时间
 * @param dt        步长
 * @param errors    [out] 同步误差时间序列 ||x - y||
 * @param max_errs   最大存储误差点数
 * @return 同步建立后的平均误差
 */
double lorenz_sync_error(double sigma, double rho, double beta,
                          double coupling, double t_end, double dt,
                          double *errors, int max_errs)
{
    double p[3] = {sigma, rho, beta};
    double x[3] = {1.0, 1.0, 1.0};
    double y[3] = {1.1, 0.9, 0.9};
    int n_steps = (int)(t_end / dt);
    int i, k, err_idx = 0;
    double avg_err = 0.0;
    int skip = n_steps / max_errs;
    if (skip < 1) skip = 1;

    double *xn = (double*)chaos_calloc(3, sizeof(double));
    double *yn = (double*)chaos_calloc(3, sizeof(double));

    for (i = 0; i < n_steps; i++) {
        /* 驱动系统 */
        double dx[3], dy[3];
        lorenz_rhs(dx, x, 0.0, p);
        lorenz_rhs(dy, y, 0.0, p);
        dy[0] += coupling * (x[0] - y[0]);

        for (k = 0; k < 3; k++) {
            xn[k] = x[k] + dt * dx[k];
            yn[k] = y[k] + dt * dy[k];
        }

        for (k = 0; k < 3; k++) { x[k] = xn[k]; y[k] = yn[k]; }

        double err = 0.0;
        for (k = 0; k < 3; k++) {
            double d = x[k] - y[k];
            err += d * d;
        }
        err = sqrt(err);
        avg_err += err;

        if (errors && err_idx < max_errs && i % skip == 0)
            errors[err_idx++] = err;
    }

    free(xn); free(yn);
    return avg_err / n_steps;
}

/* ============================================================
 * L8: 复杂网络 — 小世界与无标度
 * ============================================================ */

/**
 * Watts-Strogatz小世界网络生成
 *
 * 算法:
 *   1. 构建环状近邻网络(n个节点, 每侧k个近邻)
 *   2. 以概率p重连每条边(排除自环和重复)
 *
 * 属性:
 *   p=0: 规则网络(高聚类, 长平均路径)
 *   p=1: 随机网络(低聚类, 短平均路径)
 *   小世界区域(p≈0.01-0.1): 高聚类 + 短平均路径
 *
 * @param adj     [out] n×n邻接矩阵(调用者分配)
 * @param n       节点数
 * @param k       每侧近邻数
 * @param p       重连概率
 *
 * 参考: Watts & Strogatz (1998) Nature 393:440
 */
void watts_strogatz_network(int *adj, int n, int k, double p)
{
    int i, j;
    unsigned int seed = 31337;
    #define WS_RAND() (seed = (1103515245 * seed + 12345) & 0x7fffffff, (double)seed / 0x7fffffff)

    /* 初始环 */
    memset(adj, 0, (size_t)n * n * sizeof(int));
    for (i = 0; i < n; i++) {
        for (j = 1; j <= k; j++) {
            int neighbor = (i + j) % n;
            adj[i * n + neighbor] = 1;
            adj[neighbor * n + i] = 1;
        }
    }

    /* 重连 */
    for (i = 0; i < n; i++) {
        for (j = 1; j <= k; j++) {
            int neighbor = (i + j) % n;
            if (WS_RAND() < p) {
                /* 选择新的随机邻居 */
                int new_nb;
                do {
                    new_nb = (int)(WS_RAND() * (n - 1));
                    if (new_nb >= i) new_nb++;
                } while (new_nb == i || adj[i * n + new_nb] == 1);

                adj[i * n + neighbor] = 0;
                adj[neighbor * n + i] = 0;
                adj[i * n + new_nb] = 1;
                adj[new_nb * n + i] = 1;
            }
        }
    }
}

/**
 * 网络聚类系数
 *
 * C_i = 2*E_i / (k_i*(k_i-1))
 * 其中E_i是节点i邻居之间的边数, k_i是i的度
 *
 * C = <C_i> (平均聚类系数)
 *
 * 对于规则环网络(k=2): C = 3/4
 * 对于随机网络: C ≈ <k>/N
 * 对于小世界: C >> <k>/N
 */
double network_clustering(const int *adj, int n)
{
    double sum_C = 0.0;
    int i;

    for (i = 0; i < n; i++) {
        /* 统计i的邻居 */
        int *neighbors = (int*)chaos_calloc(n, sizeof(int));
        int k_i = 0;
        int j;
        for (j = 0; j < n; j++) {
            if (adj[i * n + j]) neighbors[k_i++] = j;
        }

        if (k_i < 2) continue;

        /* 统计邻居之间的边数 */
        int triangles = 0;
        int a, b;
        for (a = 0; a < k_i - 1; a++) {
            for (b = a + 1; b < k_i; b++) {
                if (adj[neighbors[a] * n + neighbors[b]])
                    triangles++;
            }
        }

        sum_C += (2.0 * triangles) / (k_i * (k_i - 1));
        free(neighbors);
    }

    return sum_C / n;
}

/**
 * 网络平均路径长度(Floyd-Warshall)
 *
 * L = (1/(N*(N-1))) * sum_{i!=j} d(i,j)
 *
 * @return 平均路径长度, -1=不连通
 */
double network_average_path(const int *adj, int n)
{
    int i, j, k;
    int *dist = (int*)chaos_calloc((size_t)n * n, sizeof(int));

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) dist[i * n + j] = 0;
            else if (adj[i * n + j]) dist[i * n + j] = 1;
            else dist[i * n + j] = 999999;
        }
    }

    /* Floyd-Warshall */
    for (k = 0; k < n; k++) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                int ik = dist[i * n + k];
                int kj = dist[k * n + j];
                if (ik + kj < dist[i * n + j])
                    dist[i * n + j] = ik + kj;
            }
        }
    }

    double sum_d = 0.0;
    int count = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j && dist[i * n + j] < 999999) {
                sum_d += dist[i * n + j];
                count++;
            }
        }
    }

    free(dist);
    return (count > 0) ? sum_d / count : -1.0;
}

/* ============================================================
 * L2: 瞬态Lyapunov指数 — 有限时间稳定性的刻画
 * ============================================================ */

/**
 * 有限时间Lyapunov指数(FTLE)
 *
 * 定义: lambda(tau) = (1/tau) * log(||delta x(tau)|| / ||delta x(0)||)
 *
 * 用于: Lagrangian相干结构(LCS)检测, 拉格朗日湍流分析
 *
 * @param flow_f   流场RHS
 * @param x0       初始位置
 * @param tau      有限时间间隔
 * @param dt       积分步长
 * @param delta0   初始扰动大小
 * @return FTLE
 *
 * 参考: Haller (2015) "Lagrangian Coherent Structures"
 *       Shadden, Lekien, Marsden (2005)
 */
double finite_time_lyapunov(flow_rhs_t flow_f, const double *x0,
                             int dim, const double *p,
                             double tau, double dt, double delta0)
{
    double *xA = (double*)chaos_calloc(dim, sizeof(double));
    double *xB = (double*)chaos_calloc(dim, sizeof(double));
    double *xAn = (double*)chaos_calloc(dim, sizeof(double));
    double *xBn = (double*)chaos_calloc(dim, sizeof(double));
    int i;

    for (i = 0; i < dim; i++) xA[i] = x0[i];
    for (i = 0; i < dim; i++) xB[i] = x0[i];
    xB[0] += delta0;  /* 沿x1方向扰动 */

    int n_steps = (int)(tau / dt);
    if (n_steps < 1) { free(xA); free(xB); free(xAn); free(xBn); return 0.0; }

    double t = 0.0;
    int step;
    for (step = 0; step < n_steps; step++) {
        rk4_step(xAn, flow_f, xA, t, dt, dim, p);
        rk4_step(xBn, flow_f, xB, t, dt, dim, p);
        for (i = 0; i < dim; i++) { xA[i] = xAn[i]; xB[i] = xBn[i]; }
        t += dt;
    }

    double dist = 0.0;
    for (i = 0; i < dim; i++) {
        double d = xB[i] - xA[i];
        dist += d * d;
    }
    dist = sqrt(dist);

    free(xA); free(xB); free(xAn); free(xBn);
    if (dist < 1e-300 || delta0 < 1e-300) return 0.0;
    return log(dist / delta0) / tau;
}

/**
 * 0-1混沌检测 (Gottwald & Melbourne 2004)
 *
 * 原理: 对时间序列phi, 计算辅助变量
 *   p_c(n) = sum_{j=1}^n phi_j * cos(j*c)
 *   q_c(n) = sum_{j=1}^n phi_j * sin(j*c)
 *
 * 均方位移: M_c(n) = <(p_c(j+n)-p_c(j))^2 + (q_c(j+n)-q_c(j))^2>
 *
 * 渐近增长率K_c: M_c(n) ~ K_c*n  (规则运动) 或 ~ V*n^2 (混沌)
 *
 * K = median_{c} K_c
 * K ≈ 0: 规则运动
 * K ≈ 1: 混沌运动
 *
 * @param signal  时间序列
 * @param n       长度
 * @param c       频率参数(常取1.7或黄金比例相关值)
 * @return K值 (0~1之间)
 *
 * 参考: Gottwald & Melbourne (2004) Proc.Roy.Soc.A 460:603
 */
double test_0_1_chaos(const double *signal, int n, double c)
{
    if (n < 100) return 0.0;

    /* 计算p(n)和q(n) */
    double *p_vals = (double*)chaos_calloc(n, sizeof(double));
    double *q_vals = (double*)chaos_calloc(n, sizeof(double));
    double p_sum = 0.0, q_sum = 0.0;
    int i;

    for (i = 0; i < n; i++) {
        p_sum += signal[i] * cos(c * (i + 1));
        q_sum += signal[i] * sin(c * (i + 1));
        p_vals[i] = p_sum;
        q_vals[i] = q_sum;
    }

    /* 均方位移 */
    int n_cut = n / 10;
    int max_shift = n / 2;
    double *M = (double*)chaos_calloc(max_shift, sizeof(double));
    int shift, j;

    for (shift = 1; shift <= max_shift; shift++) {
        double sum = 0.0;
        for (j = n_cut; j < n - shift; j++) {
            double dp = p_vals[j + shift] - p_vals[j];
            double dq = q_vals[j + shift] - q_vals[j];
            sum += dp * dp + dq * dq;
        }
        M[shift - 1] = sum / (n - shift - n_cut);
    }

    /* 修正均方位移: D(n) = M(n) - <phi>^2 * (1-cos(n*c))/(1-cos(c)) */
    double mean_phi = 0.0;
    for (i = 0; i < n; i++) mean_phi += signal[i];
    mean_phi /= n;

    /* 计算渐近增长率K (相关法) */
    double *xi = (double*)chaos_calloc(max_shift, sizeof(double));
    for (i = 0; i < max_shift; i++) xi[i] = (double)(i + 1);

    double mean_M = 0.0, mean_xi = 0.0;
    int n_use = max_shift / 2;
    for (i = max_shift / 4; i < max_shift * 3 / 4; i++) {
        mean_M  += M[i];
        mean_xi += xi[i];
    }
    mean_M  /= n_use;
    mean_xi /= n_use;

    double cov = 0.0, var_xi = 0.0;
    for (i = max_shift / 4; i < max_shift * 3 / 4; i++) {
        cov    += (M[i] - mean_M) * (xi[i] - mean_xi);
        var_xi += (xi[i] - mean_xi) * (xi[i] - mean_xi);
    }
    double K = (var_xi > 1e-15) ? cov / var_xi : 0.0;

    free(p_vals); free(q_vals); free(M); free(xi);
    return K;
}
