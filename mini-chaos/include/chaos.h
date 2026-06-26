/**
 * chaos.h — 混沌动力学主头文件 (Chaos Theory & Nonlinear Dynamics)
 *
 * 参考教材:
 *   Strogatz "Nonlinear Dynamics and Chaos" (2015)
 *   Goldstein, Poole, Safko "Classical Mechanics" Ch.11 (2002)
 *   May "Simple Mathematical Models" Nature (1976)
 *   Benettin et al. "Lyapunov Characteristic Exponents" Meccanica (1980)
 *
 * 知识覆盖: L1–L9 (详见 docs/knowledge-graph.md)
 * 零外部依赖 (仅标准 C99 库)
 */

#ifndef CHAOS_H
#define CHAOS_H

#include <stdlib.h>
#include <math.h>

/* ============================================================
 * L1: 核心类型定义 (Core Type Definitions)
 * ============================================================ */

/** 相空间维度上限 */
#define CHAOS_MAX_DIM 16

/**
 * 混沌动力学系统描述符
 *
 * 离散映射:  x_{n+1} = f(x_n, params)
 * 连续流:    dx/dt = f(x, t, params)
 *
 * 参考: Strogatz §2.0, §5.0
 */
typedef struct {
    int    dim;
    const char *name;
    int    is_nonautonomous;
} ChaosSystem;

/**
 * 动力学轨迹
 *
 * 按行存储: data[i*dim + k] = 第 i 个点的第 k 分量
 * 参考: Strogatz §2.2
 */
typedef struct {
    int          n_points;
    int          dim;
    double      *data;
    double       t_end;
    double       dt;
} ChaosTrajectory;

/**
 * Poincaré 截面
 *
 * 参考: Strogatz §8.7, Goldstein §11.2
 */
typedef struct {
    int          n_sections;
    int          dim;
    double      *points;
    int          section_var;
    double       section_value;
} PoincareSection;

/**
 * Lyapunov 指数谱
 *
 * 参考: Benettin et al. (1980), Wolf et al. (1985), Strogatz §9.3
 */
typedef struct {
    int          n_exponents;
    double      *exponents;      /**< 从大到小排序 */
    double       sum_positive;   /**< KS 熵上界 */
    double       kaplan_yorke;   /**< Lyapunov 维数 */
} LyapunovSpectrum;

/**
 * 分岔图数据
 *
 * 参考: Strogatz Ch.8, Kuznetsov (2004)
 */
typedef struct {
    int          n_params;
    double      *param_vals;
    int          n_steady_per_param;
    double      *steady_vals;    /**< 按行: [i*n_steady + k] */
} BifurcationDiagram;

/**
 * 分岔类型
 *
 * 参考: Strogatz §3.1–3.4, Kuznetsov Ch.2
 */
typedef enum {
    BIFURCATION_NONE = 0,
    BIFURCATION_SADDLE_NODE,     /**< 鞍结点分岔 */
    BIFURCATION_TRANSCRITICAL,   /**< 跨临界分岔 */
    BIFURCATION_PITCHFORK,       /**< 叉形分岔 */
    BIFURCATION_PERIOD_DOUBLING, /**< 倍周期分岔 */
    BIFURCATION_HOPF,            /**< Hopf 分岔 */
    BIFURCATION_NEIMARK_SACKER   /**< Neimark-Sacker 分岔 */
} BifurcationType;

/**
 * 分岔点信息
 */
typedef struct {
    double           param_value;
    BifurcationType  type;
    int              period_before;
    int              period_after;   /**< 0 = 混沌 */
} BifurcationPoint;

/** 不动点稳定性类型 */
typedef enum {
    STABLE_NODE = 0,
    UNSTABLE_NODE,
    SADDLE_POINT,
    STABLE_FOCUS,
    UNSTABLE_FOCUS,
    CENTER,
    DEGENERATE
} FixedPointType;

/**
 * 分形图像
 *
 * 参考: Peitgen & Richter (1986), Falconer (2013)
 */
typedef struct {
    int          width, height;
    int         *data;           /**< data[y*width + x] */
    double       x_min, x_max;
    double       y_min, y_max;
    int          max_iter;
} FractalImage;

/**
 * 延迟嵌入重构
 *
 * 参考: Takens (1981), Kantz & Schreiber (2004)
 */
typedef struct {
    int          n_points;
    int          embed_dim;      /**< 嵌入维数 m */
    int          delay;          /**< 时间延迟 τ */
    double      *points;         /**< 按行存储 */
} Embedding;

/**
 * 递归图
 *
 * 参考: Eckmann et al. (1987), Marwan et al. (2007)
 */
typedef struct {
    int          n;
    int         *matrix;         /**< n×n 二值矩阵 */
    double       epsilon;
    double       recurrence_rate;
} RecurrencePlot;

/* ============================================================
 * L4: 基本常数
 * ============================================================ */

/** Feigenbaum δ: 倍周期分岔间距比极限 */
#define FEIGENBAUM_DELTA  4.669201609102990

/** Feigenbaum α: 分岔宽度比极限 */
#define FEIGENBAUM_ALPHA  2.502907875095892

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_2PI
#define M_2PI 6.28318530717958647692
#endif

/**
 * Sharkovskii 序 — Li-Yorke "Period Three Implies Chaos" (1975)
 *
 * 序: 3 ▷ 5 ▷ 7 ▷ ... ▷ 2·3 ▷ 2·5 ▷ ... ▷ 2²·3 ▷ ... ▷ 2³ ▷ 2² ▷ 2 ▷ 1
 * 定理: 连续映射有周期 k 且 k ▷ m ⇒ 必有周期 m
 */
#define CHAOS_LYAPUNOV_THRESHOLD 0.005

/* ============================================================
 * 内存管理
 * ============================================================ */

void *chaos_calloc(size_t nmemb, size_t size);
void  chaos_free(void *ptr);
void  chaos_trajectory_free(ChaosTrajectory *traj);
void  chaos_lyapunov_free(LyapunovSpectrum *spec);
void  chaos_bifurcation_free(BifurcationDiagram *bd);
void  chaos_fractal_free(FractalImage *img);
void  chaos_embedding_free(Embedding *emb);
void  chaos_recurrence_free(RecurrencePlot *rp);
void  chaos_poincare_free(PoincareSection *ps);

#endif /* CHAOS_H */
