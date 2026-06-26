/**
 * chaos_embedding.c — 相空间重构与时间序列分析
 *
 * 涵盖:
 *   L3: Takens嵌入定理 — 延迟坐标重构
 *   L5: AMI/FNN — 嵌入参数选择
 *   L2: 关联维数 — Grassberger-Procaccia算法
 *   L7: 替代数据检验, 非线性预测
 *
 * 参考:
 *   Takens (1981), Kantz & Schreiber (2004)
 *   Grassberger & Procaccia (1983)
 */

#include "chaos_embedding.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

/* ============================================================
 * L3: Takens嵌入定理 — 相空间重构
 * ============================================================ */

/**
 * 时间延迟嵌入
 *
 * X(t_i) = [s(t_i), s(t_i-tau), s(t_i-2*tau), ..., s(t_i-(m-1)*tau)]
 *
 * 其中:
 *   s 是标量观测时间序列
 *   tau 是时间延迟(采样步数)
 *   m 是嵌入维数
 *
 * Takens定理(1981):
 *   对于d维紧致流形M上的C^2动力学系统和C^2的观测函数h,
 *   对通有的系统(f,h)和通有的tau,
 *   当m > 2d时, 延迟映射Phi: M -> R^m 是嵌入(微分同胚子其像)
 *
 * 推论: 重构的相空间与原系统拓扑等价
 *
 * 扩展(Sauer, Yorke, Casdagli 1991): m > 2*D_0(盒计数维数) 即足
 *
 * @param signal  标量时间序列
 * @param n       信号长度
 * @param tau     延迟
 * @param m       嵌入维数
 * @return Embedding (调用者free)
 *
 * 复杂度: O((n - (m-1)*tau) * m)
 */
Embedding* time_delay_embedding(const double *signal, int n,
                                 int tau, int m)
{
    int N = n - (m - 1) * tau;
    if (N <= 0 || m < 2 || tau < 1) return NULL;

    Embedding *emb = (Embedding*)chaos_calloc(1, sizeof(Embedding));
    emb->n_points  = N;
    emb->embed_dim = m;
    emb->delay     = tau;
    emb->points    = (double*)chaos_calloc((size_t)N * m, sizeof(double));

    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < m; j++) {
            emb->points[i * m + j] = signal[i + j * tau];
        }
    }
    return emb;
}

/* ============================================================
 * L5: 嵌入参数选择 — AMI + FNN
 * ============================================================ */

/**
 * 平均互信息(AMI)
 *
 * 定义: I(tau) = sum_{s_n, s_{n+tau}} P(s_n, s_{n+tau})
 *                * log[P(s_n, s_{n+tau}) / (P(s_n) * P(s_{n+tau}))]
 *
 * 解释: 互信息度量了两个时间点之间共享的信息量
 *       首个AMI极小值对应最优延迟tau
 *
 * 原理: tau太小 → 冗余信息, 相空间沿对角线压缩
 *       tau太大 → 信息丢失, 相空间过度展开
 *
 * @param signal  时间序列
 * @param n       信号长度
 * @param tau     延迟
 * @param n_bins  直方图箱数
 * @return 互信息(nats)
 *
 * 参考: Fraser & Swinney (1986) PRA 33:1134
 *       Kantz & Schreiber §3.3
 *
 * 复杂度: O(n + n_bins^2)
 */
double average_mutual_information(const double *signal, int n,
                                   int tau, int n_bins)
{
    if (tau >= n || n_bins < 2) return 0.0;

    int N = n - tau;
    if (N <= 0) return 0.0;

    /* 找到数据范围 */
    double s_min = signal[0], s_max = signal[0];
    int i;
    for (i = 0; i < n; i++) {
        if (signal[i] < s_min) s_min = signal[i];
        if (signal[i] > s_max) s_max = signal[i];
    }
    double range = s_max - s_min;
    if (range < 1e-15) return 0.0;
    double bin_width = range / n_bins;

    /* 构建二维直方图 */
    int *hist2d = (int*)chaos_calloc((size_t)n_bins * n_bins, sizeof(int));
    for (i = 0; i < N; i++) {
        int bx = (int)((signal[i]     - s_min) / bin_width);
        int by = (int)((signal[i + tau] - s_min) / bin_width);
        if (bx >= n_bins) bx = n_bins - 1;
        if (by >= n_bins) by = n_bins - 1;
        if (bx < 0) bx = 0;
        if (by < 0) by = 0;
        hist2d[by * n_bins + bx]++;
    }

    /* 计算概率和互信息 */
    double ami = 0.0;
    double inv_N = 1.0 / N;
    int bx, by;

    /* 边缘概率 */
    double *px = (double*)chaos_calloc(n_bins, sizeof(double));
    double *py = (double*)chaos_calloc(n_bins, sizeof(double));
    for (by = 0; by < n_bins; by++) {
        for (bx = 0; bx < n_bins; bx++) {
            double p = hist2d[by * n_bins + bx] * inv_N;
            px[bx] += p;
            py[by] += p;
        }
    }

    for (by = 0; by < n_bins; by++) {
        for (bx = 0; bx < n_bins; bx++) {
            double p = hist2d[by * n_bins + bx] * inv_N;
            if (p > 0.0 && px[bx] > 0.0 && py[by] > 0.0) {
                ami += p * log(p / (px[bx] * py[by]));
            }
        }
    }

    free(hist2d); free(px); free(py);
    return ami;
}

/**
 * 寻找最优延迟tau (首个AMI极小值)
 *
 * 扫描tau=1..tau_max, 返回AMI第一个局部极小值处的tau
 *
 * @param tau_max  搜索上界
 * @return 最优tau
 */
int find_optimal_tau(const double *signal, int n, int tau_max)
{
    if (tau_max < 3 || n < tau_max + 10) return 1;

    double *amis = (double*)chaos_calloc(tau_max, sizeof(double));
    int t, best = 1;

    for (t = 1; t <= tau_max; t++)
        amis[t - 1] = average_mutual_information(signal, n, t, 20);

    double min_ami = amis[0];
    for (t = 2; t < tau_max; t++) {
        /* 局部极小: amis[t-2] > amis[t-1] 且 amis[t-1] < amis[t] */
        if (amis[t - 1] < amis[t - 2] && amis[t - 1] < amis[t]) {
            best = t;
            break;
        }
        if (amis[t - 1] < min_ami) {
            min_ami = amis[t - 1];
            best = t;
        }
    }

    free(amis);
    return best;
}

/**
 * 假近邻法(FNN) — 确定最小嵌入维数m
 *
 * 原理: 在过低维嵌入中, 投影会造成"假近邻"
 *       — 在低维空间中看似相邻, 在高一维中远离
 *
 * 判据(Kennel, Brown, Abarbanel 1992):
 *   对每个点i, 找最近邻j:
 *     R_i = |s_{i+m*tau} - s_{j+m*tau}| / |X_i - X_j|
 *   如果R_i > R_tol, (i,j)是假近邻
 *
 * 最优m: 假近邻比例降至接近0的最小维数
 *
 * @param fnn_ratio [out] 每个m的假近邻比例(长度m_max)
 * @return 推荐的最优m(首个fnn_ratio<0.01的m值)
 *
 * 参考: Kennel, Brown, Abarbanel (1992) PRA 45:3403
 *
 * 复杂度: O(m_max * N^2 * m)
 */
int false_nearest_neighbors(const double *signal, int n, int tau,
                             int m_max, double R_tol,
                             double *fnn_ratio)
{
    int m;
    for (m = 1; m <= m_max; m++) {
        Embedding *emb = time_delay_embedding(signal, n, tau, m);
        if (!emb) { fnn_ratio[m - 1] = 1.0; continue; }

        int N = emb->n_points;
        int fnn_count = 0;
        int i, j;

        for (i = 0; i < N; i++) {
            /* 找最近邻 */
            double min_dist = 1e100;
            int nn_idx = -1;
            for (j = 0; j < N; j++) {
                if (i == j) continue;
                double dist2 = 0.0;
                int k;
                for (k = 0; k < m; k++) {
                    double diff = emb->points[i * m + k] - emb->points[j * m + k];
                    dist2 += diff * diff;
                }
                if (dist2 < min_dist) {
                    min_dist = dist2;
                    nn_idx = j;
                }
            }

            if (nn_idx >= 0 && min_dist > 1e-15) {
                /* 检查m+1维中的距离扩大 */
                int idx_m1 = i + m * tau;    /* 增加一维后的坐标 */
                int nid_m1 = nn_idx + m * tau;
                if (idx_m1 < n && nid_m1 < n) {
                    double diff_m1 = signal[idx_m1] - signal[nid_m1];
                    double R = fabs(diff_m1) / sqrt(min_dist);
                    if (R > R_tol) fnn_count++;
                }
            }
        }

        fnn_ratio[m - 1] = (double)fnn_count / N;
        chaos_embedding_free(emb);
        free(emb);
    }

    /* 返回首个fnn_ratio < 0.01的m */
    for (m = 0; m < m_max; m++) {
        if (fnn_ratio[m] < 0.01)
            return m + 1;
    }
    return m_max;  /* 回退到最大测试值 */
}

/**
 * 关联维数 (Grassberger-Procaccia)
 *
 * D_2 = lim_{eps->0} log(C(eps)) / log(eps)
 *
 * @param points  延迟嵌入点集
 * @return D_2估计
 *
 * 复杂度: O(n_eps * n^2 * dim)
 */
double correlation_dimension(const double *points, int n, int dim,
                              const double *epsilons, int n_eps)
{
    if (n < 2 || n_eps < 2 || dim < 1) return 0.0;

    double *log_C = (double*)chaos_calloc(n_eps, sizeof(double));
    double *log_eps = (double*)chaos_calloc(n_eps, sizeof(double));
    int ie, i, j;

    for (ie = 0; ie < n_eps; ie++) {
        double eps = epsilons[ie];
        if (eps <= 0.0) { log_C[ie] = 0.0; log_eps[ie] = 0.0; continue; }

        long long count = 0;
        for (i = 0; i < n - 1; i++) {
            for (j = i + 1; j < n; j++) {
                double dist2 = 0.0;
                int k;
                for (k = 0; k < dim; k++) {
                    double diff = points[i * dim + k] - points[j * dim + k];
                    dist2 += diff * diff;
                }
                if (dist2 < eps * eps) count++;
            }
        }

        double C_val = (2.0 * count) / ((double)n * (n - 1));
        if (C_val < 1e-300) C_val = 1e-300;
        log_C[ie] = log(C_val);
        log_eps[ie] = log(eps);
    }

    /* 线性回归: log(C) ~ D_2 * log(eps) */
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

/* ============================================================
 * L7: 替代数据检验 — 确定性vs随机性
 * ============================================================ */

/**
 * 替代数据 — 随机重排(保持幅度分布)
 *
 * Fisher-Yates洗牌: O(n)
 *
 * 应用: 零假设 H0 = "数据来自独立同分布随机变量"
 * 如果原始数据的某个统计量(如预测误差)与替代数据显著不同,
 * 则拒绝H0 → 存在确定性结构
 *
 * 参考: Theiler et al. (1992) Physica D 58:77
 */
void surrogate_data_shuffle(double *surrogate, const double *signal, int n)
{
    int i;
    for (i = 0; i < n; i++)
        surrogate[i] = signal[i];

    /* Fisher-Yates shuffle */
    unsigned int seed = 7777;
    for (i = n - 1; i > 0; i--) {
        seed = (1103515245 * seed + 12345) & 0x7fffffff;
        int j = (int)(((double)seed / 0x7fffffff) * (i + 1));
        if (j > i) j = i;
        double tmp = surrogate[i];
        surrogate[i] = surrogate[j];
        surrogate[j] = tmp;
    }
}

/**
 * 非线性预测误差
 *
 * 方法(Sugihara & May 1990):
 *   1. 构建延迟嵌入相空间
 *   2. 对每个参考点, 找k个最近邻
 *   3. 用近邻的T_pred步演化外推预测值
 *   4. 预测误差 = <|x_pred - x_actual|> / 信号标准差
 *
 * 解释: 对确定论混沌, 预测误差随嵌入维数m增加而降低
 *       对随机噪声, 预测误差不受m影响
 *
 * @param T_pred  向前预测步数
 * @return 归一化预测误差
 */
double nonlinear_prediction_error(const double *signal, int n,
                                   int tau, int m, int T_pred)
{
    if (n < 100 || m < 2 || T_pred < 1) return 1.0;

    /* 构建嵌入 */
    Embedding *emb = time_delay_embedding(signal, n, tau, m);
    if (!emb) return 1.0;
    int N = emb->n_points;

    /* 计算信号std */
    double mean_s = 0.0;
    int i;
    for (i = 0; i < n; i++) mean_s += signal[i];
    mean_s /= n;
    double var_s = 0.0;
    for (i = 0; i < n; i++) {
        double d = signal[i] - mean_s;
        var_s += d * d;
    }
    var_s /= n;
    double std_s = sqrt(var_s);
    if (std_s < 1e-15) { chaos_embedding_free(emb); free(emb); return 0.0; }

    /* 对每个参考点做预测 */
    double err_sum = 0.0;
    int count = 0;
    int ref;

    for (ref = 0; ref < N - T_pred; ref++) {
        /* 找最近邻 */
        double min_dist = 1e100;
        int nn = -1;
        int j;
        for (j = 0; j < N; j++) {
            if (j == ref) continue;
            double dist2 = 0.0;
            int k;
            for (k = 0; k < m; k++) {
                double diff = emb->points[ref * m + k] - emb->points[j * m + k];
                dist2 += diff * diff;
            }
            if (dist2 < min_dist) { min_dist = dist2; nn = j; }
        }

        if (nn >= 0) {
            /* 用近邻的演化预测ref的未来 */
            int idx_ref = ref + (m - 1) * tau + T_pred;
            int idx_nn  = nn  + (m - 1) * tau + T_pred;
            if (idx_ref < n && idx_nn < n) {
                double predicted = signal[idx_nn];
                double actual    = signal[idx_ref];
                err_sum += fabs(predicted - actual) / std_s;
                count++;
            }
        }
    }

    chaos_embedding_free(emb);
    free(emb);

    return (count > 0) ? err_sum / count : 1.0;
}
