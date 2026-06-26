/**
 * chaos_types.c — 核心类型与内存管理实现
 *
 * 按 SKILL.md 规范: 每个函数实现独立知识点
 * 零外部依赖 (仅 stdlib)
 */

#include "chaos.h"
#include <stdio.h>
#include <string.h>

/* ──────────────────────────────────────────────────────────
 * L1: 内存分配与释放
 * ────────────────────────────────────────────────────────── */

/**
 * 安全 calloc — 检测分配失败
 *
 * 知识点: 防御性编程, 失败时打印错误并退出
 * 复杂度: O(1)
 */
void *chaos_calloc(size_t nmemb, size_t size)
{
    size_t total = nmemb * size;
    if (total == 0) return NULL;
    void *ptr = calloc(nmemb, size);
    if (!ptr) {
        fprintf(stderr, "chaos_calloc: failed to allocate %zu bytes\n", total);
        abort();
    }
    return ptr;
}

/**
 * 安全 free — 封装以备未来扩展 (如内存池)
 *
 * 知识点: RAII 概念在 C 中的模拟
 * 复杂度: O(1)
 */
void chaos_free(void *ptr)
{
    if (ptr) free(ptr);
}

/**
 * 释放轨迹内部数据
 *
 * 知识点: 结构体深释放模式
 * 复杂度: O(1)
 */
void chaos_trajectory_free(ChaosTrajectory *traj)
{
    if (!traj) return;
    if (traj->data) { free(traj->data); traj->data = NULL; }
    traj->n_points = 0;
    traj->dim = 0;
}

/**
 * 释放 Lyapunov 谱
 */
void chaos_lyapunov_free(LyapunovSpectrum *spec)
{
    if (!spec) return;
    if (spec->exponents) { free(spec->exponents); spec->exponents = NULL; }
    spec->n_exponents = 0;
}

/**
 * 释放分岔图数据
 */
void chaos_bifurcation_free(BifurcationDiagram *bd)
{
    if (!bd) return;
    if (bd->param_vals)  { free(bd->param_vals);  bd->param_vals  = NULL; }
    if (bd->steady_vals) { free(bd->steady_vals); bd->steady_vals = NULL; }
    bd->n_params = 0;
    bd->n_steady_per_param = 0;
}

/**
 * 释放分形图像
 */
void chaos_fractal_free(FractalImage *img)
{
    if (!img) return;
    if (img->data) { free(img->data); img->data = NULL; }
    img->width  = 0;
    img->height = 0;
}

/**
 * 释放嵌入结构
 */
void chaos_embedding_free(Embedding *emb)
{
    if (!emb) return;
    if (emb->points) { free(emb->points); emb->points = NULL; }
    emb->n_points  = 0;
    emb->embed_dim = 0;
    emb->delay     = 0;
}

/**
 * 释放递归图
 */
void chaos_recurrence_free(RecurrencePlot *rp)
{
    if (!rp) return;
    if (rp->matrix) { free(rp->matrix); rp->matrix = NULL; }
    rp->n = 0;
    rp->epsilon = 0.0;
    rp->recurrence_rate = 0.0;
}

/**
 * 释放 Poincaré 截面
 */
void chaos_poincare_free(PoincareSection *ps)
{
    if (!ps) return;
    if (ps->points) { free(ps->points); ps->points = NULL; }
    ps->n_sections = 0;
    ps->dim = 0;
}
