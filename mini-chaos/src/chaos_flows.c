/**
 * chaos_flows.c — 连续时间混沌系统 + RK4 积分器
 *
 * 涵盖:
 *   L4: 基本定律 — Lorenz方程, Duffing方程
 *   L6: 经典系统 — Lorenz/Rossler/Chua/Duffing/VdP/Sprott/Chen/H-H
 *   L5: 计算方法 — RK4, 自适应步长, Poincare截面
 *
 * 参考: Lorenz (1963), Rossler (1976), Strogatz Ch.9
 */

#include "chaos_flows.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

/* ============================================================
 * L6: 经典混沌系统 — 右侧函数
 * ============================================================ */

/**
 * Lorenz系统: dx/dt=sigma(y-x), dy/dt=x(rho-z)-y, dz/dt=xy-beta*z
 *
 * 参数: p[0]=sigma=10, p[1]=rho=28, p[2]=beta=8/3
 *
 * 不动点 (rho>1):
 *   C+ = (sqrt(beta*(rho-1)), sqrt(beta*(rho-1)), rho-1)
 *   C- = (-sqrt(beta*(rho-1)), -sqrt(beta*(rho-1)), rho-1)
 *   O  = (0, 0, 0) — 鞍点
 *
 * 知识点: 确定性非周期流, 蝴蝶效应, Rayleigh-Benard对流
 * MIT 8.07 / Caltech Ph 106 / Cambridge Part III
 */
void lorenz_rhs(double dx[3], const double x[3], double t, const double p[3])
{
    (void)t;
    double sigma = p[0], rho = p[1], beta = p[2];
    dx[0] = sigma * (x[1] - x[0]);
    dx[1] = x[0] * (rho - x[2]) - x[1];
    dx[2] = x[0] * x[1] - beta * x[2];
}

/**
 * Rossler系统: dx/dt=-y-z, dy/dt=x+ay, dz/dt=b+z(x-c)
 *
 * 参数: p[0]=a=0.2, p[1]=b=0.2, p[2]=c=5.7
 *
 * 知识点: 螺旋混沌, 漏斗吸引子, 鞍焦点同宿轨道
 * Stanford PHYSICS 230 / ETH 402-0800
 */
void rossler_rhs(double dx[3], const double x[3], double t, const double p[3])
{
    (void)t;
    double a = p[0], b = p[1], c = p[2];
    dx[0] = -x[1] - x[2];
    dx[1] =  x[0] + a * x[1];
    dx[2] =  b + x[2] * (x[0] - c);
}

/**
 * Chua电路 — 第一个物理实现的混沌电路
 *
 * 方程:
 *   dx/dt = alpha*(y - x - f(x))
 *   dy/dt = x - y + z
 *   dz/dt = -beta*y
 *
 * f(x) = m1*x + 0.5*(m0-m1)*(|x+1|-|x-1|)  (Chua二极管)
 *
 * 参数: p[0]=alpha=15.6, p[1]=beta=28, p[2]=m0=-1.143, p[3]=m1=-0.714
 *
 * 知识点: 分段线性混沌, 双涡卷吸引子, 物理实现
 * Berkeley PHYS 242: 非线性电路
 */
void chua_rhs(double dx[3], const double x[3], double t, const double p[5])
{
    (void)t;
    double alpha = p[0], beta = p[1];
    double m0 = p[2], m1 = p[3];
    double fx = m1 * x[0]
              + 0.5 * (m0 - m1) * (fabs(x[0] + 1.0) - fabs(x[0] - 1.0));
    dx[0] = alpha * (x[1] - x[0] - fx);
    dx[1] = x[0] - x[1] + x[2];
    dx[2] = -beta * x[1];
}

/**
 * Duffing振子: 非线性弹性 + 受迫阻尼
 *
 * 自治化形式:
 *   dx/dt = v
 *   dv/dt = -delta*v + beta*x - alpha*x^3 + gamma*cos(omega*tau)
 *   dtau/dt = 1
 *
 * 参数: p[0]=delta=0.2, p[1]=beta=1, p[2]=alpha=1, p[3]=gamma=0.3, p[4]=omega=1.0
 *
 * 物理背景: 磁弹性梁, Josephson结, 船舶横摇
 * 知识点: 双稳态, 混沌鞍, Melnikov方法
 * Goldstein §11.5 / MIT 8.012
 */
void duffing_rhs(double dx[3], const double x[3], double t, const double p[5])
{
    (void)t;
    double delta = p[0], beta0 = p[1], alpha = p[2];
    double gamma = p[3], omega = p[4];
    dx[0] = x[1];
    dx[1] = -delta * x[1] + beta0 * x[0]
            - alpha * x[0] * x[0] * x[0]
            + gamma * cos(omega * x[2]);
    dx[2] = 1.0;
}

/**
 * 受迫Van der Pol振子
 *
 * 方程: dx/dt=v, dv/dt=mu*(1-x^2)*v - x + F*cos(omega*tau), dtau/dt=1
 *
 * 参数: p[0]=mu=3, p[1]=F=1.2, p[2]=omega=1.0
 *
 * 知识点: 自激振荡, 极限环, 频率夹带, 准周期->混沌
 * Strogatz §8.5 / Oxford CMT
 */
void forced_vdp_rhs(double dx[3], const double x[3], double t, const double p[4])
{
    (void)t;
    double mu = p[0], F = p[1], omega = p[2];
    dx[0] = x[1];
    dx[1] = mu * (1.0 - x[0] * x[0]) * x[1] - x[0]
            + F * cos(omega * x[2]);
    dx[2] = 1.0;
}

/**
 * Sprott B — 最简混沌流 (5项二次)
 *
 * dx/dt = y*z,  dy/dt = x-y,  dz/dt = 1-x*y
 *
 * 知识点: 计算机搜索最简混沌系统
 * Sprott (1994) "Some Simple Chaotic Flows"
 */
void sprott_b_rhs(double dx[3], const double x[3], double t, const double p[1])
{
    (void)t; (void)p;
    dx[0] = x[1] * x[2];
    dx[1] = x[0] - x[1];
    dx[2] = 1.0 - x[0] * x[1];
}

/**
 * Sprott C — 最简混沌流 (4项二次)
 *
 * dx/dt = y*z,  dy/dt = x-y,  dz/dt = 1-x^2
 */
void sprott_c_rhs(double dx[3], const double x[3], double t, const double p[1])
{
    (void)t; (void)p;
    dx[0] = x[1] * x[2];
    dx[1] = x[0] - x[1];
    dx[2] = 1.0 - x[0] * x[0];
}

/**
 * Chen系统 — 对偶Lorenz
 *
 * dx/dt = a*(y-x)
 * dy/dt = (c-a)*x - x*z + c*y
 * dz/dt = x*y - b*z
 *
 * 参数: p[0]=a=35, p[1]=b=3, p[2]=c=28
 *
 * 知识点: 反混沌, 拓扑不等价于Lorenz
 * Chen & Ueta (1999)
 */
void chen_rhs(double dx[3], const double x[3], double t, const double p[3])
{
    (void)t;
    double a = p[0], b = p[1], c = p[2];
    dx[0] = a * (x[1] - x[0]);
    dx[1] = (c - a) * x[0] - x[0] * x[2] + c * x[1];
    dx[2] = x[0] * x[1] - b * x[2];
}

/**
 * Henon-Heiles系统 — Hamiltonian混沌
 *
 * H = 0.5*(px^2+py^2) + 0.5*(x^2+y^2) + x^2*y - y^3/3
 *
 * 方程:
 *   dx/dt = px,  dy/dt = py
 *   dpx/dt = -x - 2*x*y
 *   dpy/dt = -y - x^2 + y^2
 *
 * 知识点: 可积vs非可积, 庞加莱截面, KAM定理数值验证
 * E < 1/6: 近可积; E > 1/6: 混沌海
 *
 * Goldstein §11.3 / Caltech Ph 205
 */
void henon_heiles_rhs(double dx[4], const double x[4], double t, const double p[1])
{
    (void)t; (void)p;
    dx[0] = x[2];
    dx[1] = x[3];
    dx[2] = -x[0] - 2.0 * x[0] * x[1];
    dx[3] = -x[1] - x[0] * x[0] + x[1] * x[1];
}

/**
 * 4D超混沌Rossler — 两个正Lyapunov指数
 *
 * dx/dt = -y-z,  dy/dt = x+0.25*y+w
 * dz/dt = 3+x*z,  dw/dt = -0.5*z+0.05*w
 *
 * 知识点: 超混沌(多个正Lyapunov), 高维吸引子
 */
void hyperchaos_rossler_rhs(double dx[4], const double x[4], double t, const double p[1])
{
    (void)t; (void)p;
    dx[0] = -x[1] - x[2];
    dx[1] =  x[0] + 0.25 * x[1] + x[3];
    dx[2] =  3.0 + x[0] * x[2];
    dx[3] = -0.5 * x[2] + 0.05 * x[3];
}

/* ============================================================
 * L5: 数值积分 — RK4
 * ============================================================ */

/**
 * RK4单步积分 (4阶Runge-Kutta, Butcher 1901)
 *
 * k1 = h * f(t, x)
 * k2 = h * f(t+h/2, x+k1/2)
 * k3 = h * f(t+h/2, x+k2/2)
 * k4 = h * f(t+h, x+k3)
 * x_next = x + (k1 + 2*k2 + 2*k3 + k4)/6
 *
 * 局部截断误差: O(h^5)
 * 全局误差: O(h^4)
 *
 * 参考: Butcher (2008), Hairer et al. (1993)
 * 复杂度: O(dim) 每步, 4次RHS求值
 */
void rk4_step(double *x_next, flow_rhs_t f, const double *x,
              double t, double h, int dim, const double *p)
{
    int i;
    double *k1 = (double*)chaos_calloc(dim, sizeof(double));
    double *k2 = (double*)chaos_calloc(dim, sizeof(double));
    double *k3 = (double*)chaos_calloc(dim, sizeof(double));
    double *k4 = (double*)chaos_calloc(dim, sizeof(double));
    double *xtmp = (double*)chaos_calloc(dim, sizeof(double));

    /* k1 = h * f(t, x) */
    f(k1, x, t, p);
    for (i = 0; i < dim; i++) { k1[i] *= h; xtmp[i] = x[i] + 0.5 * k1[i]; }

    /* k2 = h * f(t+h/2, x+k1/2) */
    f(k2, xtmp, t + 0.5*h, p);
    for (i = 0; i < dim; i++) { k2[i] *= h; xtmp[i] = x[i] + 0.5 * k2[i]; }

    /* k3 = h * f(t+h/2, x+k2/2) */
    f(k3, xtmp, t + 0.5*h, p);
    for (i = 0; i < dim; i++) { k3[i] *= h; xtmp[i] = x[i] + k3[i]; }

    /* k4 = h * f(t+h, x+k3) */
    f(k4, xtmp, t + h, p);
    for (i = 0; i < dim; i++) k4[i] *= h;

    /* x_next = x + (k1 + 2*k2 + 2*k3 + k4)/6 */
    for (i = 0; i < dim; i++)
        x_next[i] = x[i] + (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;

    free(k1); free(k2); free(k3); free(k4); free(xtmp);
}

/**
 * 完整轨迹积分 (固定步长)
 *
 * @param f     右侧函数
 * @param x0    初值 (长度 dim)
 * @param dim   维度
 * @param p     参数
 * @param t_end 终止时间
 * @param dt    步长
 * @return 轨迹 (调用者负责free)
 *
 * 复杂度: O(n_steps * dim), n_steps = t_end/dt
 */
ChaosTrajectory* integrate_flow(flow_rhs_t f, const double *x0,
                                 int dim, const double *p,
                                 double t_end, double dt)
{
    int n_steps = (int)(t_end / dt);
    if (n_steps < 1) n_steps = 1;

    ChaosTrajectory *traj = (ChaosTrajectory*)chaos_calloc(1, sizeof(ChaosTrajectory));
    traj->n_points = n_steps + 1;
    traj->dim      = dim;
    traj->t_end    = t_end;
    traj->dt       = dt;
    traj->data     = (double*)chaos_calloc((size_t)(n_steps + 1) * dim, sizeof(double));

    int i, k;
    for (k = 0; k < dim; k++) traj->data[k] = x0[k];

    double *x_cur  = (double*)chaos_calloc(dim, sizeof(double));
    double *x_next = (double*)chaos_calloc(dim, sizeof(double));
    for (k = 0; k < dim; k++) x_cur[k] = x0[k];

    double t = 0.0;
    for (i = 0; i < n_steps; i++) {
        rk4_step(x_next, f, x_cur, t, dt, dim, p);
        for (k = 0; k < dim; k++) {
            x_cur[k] = x_next[k];
            traj->data[(i + 1) * dim + k] = x_cur[k];
        }
        t += dt;
    }

    free(x_cur); free(x_next);
    return traj;
}

/**
 * Poincare截面提取
 *
 * 截面条件: x[section_dim]穿过section_value且符号变化
 *
 * 应用:
 *   - Lorenz的z=rho-1截面 -> Lorenz映射
 *   - Henon-Heiles的y=0, py>0截面 -> 庞加莱截面图
 *   - Duffing的周期采样 -> 闪频映射
 *
 * 知识点: 降维分析, 周期轨道检测, KAM环面可视化
 * Strogatz §8.7 / Goldstein §11.2
 */
PoincareSection* poincare_section_extract(flow_rhs_t f, const double *x0,
                                           int dim, const double *p,
                                           double t_end, double dt,
                                           int section_dim,
                                           double section_value)
{
    ChaosTrajectory *traj = integrate_flow(f, x0, dim, p, t_end, dt);
    int i, k;

    /* 统计截面穿过次数 */
    int n_sec = 0;
    for (i = 1; i < traj->n_points; i++) {
        double v_prev = traj->data[(i - 1) * dim + section_dim] - section_value;
        double v_curr = traj->data[i * dim + section_dim] - section_value;
        if (v_prev * v_curr < 0.0) n_sec++;
    }

    PoincareSection *ps = (PoincareSection*)chaos_calloc(1, sizeof(PoincareSection));
    ps->n_sections    = n_sec;
    ps->dim           = dim;
    ps->section_var   = section_dim;
    ps->section_value = section_value;
    ps->points        = (double*)chaos_calloc((size_t)n_sec * dim, sizeof(double));

    int sec_idx = 0;
    for (i = 1; i < traj->n_points && sec_idx < n_sec; i++) {
        double v_prev = traj->data[(i - 1) * dim + section_dim] - section_value;
        double v_curr = traj->data[i * dim + section_dim] - section_value;
        if (v_prev * v_curr < 0.0) {
            for (k = 0; k < dim; k++)
                ps->points[sec_idx * dim + k] = traj->data[i * dim + k];
            sec_idx++;
        }
    }

    chaos_trajectory_free(traj);
    free(traj);
    return ps;
}

/**
 * 自适应步长RK4 (Richardson外推)
 *
 * 局部误差: err = ||x_{h/2,h/2} - x_h|| / (2^4 - 1)
 * 步长更新: h_new = h * 0.9 * (tol/err)^{1/5}
 *
 * 知识点: 自适应步长控制, 误差估计, 刚性问题回避
 *
 * @param traj_out 如非NULL, 存储轨迹 (需分配(max_steps+1)*dim doublest)
 * @return 实际步数
 */
int rk4_adaptive(double *traj_out, flow_rhs_t f, const double *x0,
                 int dim, const double *p, double t0, double t_end,
                 double dt_init, double tol, int max_steps)
{
    double t = t0;
    double h = dt_init;
    int step = 0;
    double *x   = (double*)chaos_calloc(dim, sizeof(double));
    double *xh  = (double*)chaos_calloc(dim, sizeof(double));
    double *xh1 = (double*)chaos_calloc(dim, sizeof(double));
    double *xh2 = (double*)chaos_calloc(dim, sizeof(double));
    int i;

    for (i = 0; i < dim; i++) {
        x[i] = x0[i];
        if (traj_out) traj_out[step * dim + i] = x[i];
    }

    while (t < t_end && step < max_steps) {
        if (t + h > t_end) h = t_end - t;
        if (h < 1e-12) break;

        /* 一步h, 两步h/2 */
        rk4_step(xh,  f, x, t, h, dim, p);
        rk4_step(xh1, f, x, t, h/2.0, dim, p);
        rk4_step(xh2, f, xh1, t + h/2.0, h/2.0, dim, p);

        /* 误差估计 */
        double err = 0.0;
        for (i = 0; i < dim; i++) {
            double diff = xh2[i] - xh[i];
            err += diff * diff;
        }
        err = sqrt(err) / 15.0;  /* 2^4 - 1 */

        if (err < tol || h < dt_init * 0.01) {
            for (i = 0; i < dim; i++) x[i] = xh2[i];
            t += h;
            step++;
            if (traj_out) {
                for (i = 0; i < dim; i++)
                    traj_out[step * dim + i] = x[i];
            }
        }

        /* 步长控制 */
        if (err > 0)
            h *= 0.9 * pow(tol / err, 0.2);  /* 1/(p+1) = 1/5 */
        else
            h *= 1.5;
        if (h > dt_init * 10.0) h = dt_init * 10.0;
        if (h < dt_init * 1e-4) h = dt_init * 1e-4;
    }

    free(x); free(xh); free(xh1); free(xh2);
    return step;
}

/**
 * 渐近解采样 — 用于Poincare截面/分岔图数据
 *
 * 先积分t_transient时间丢弃瞬态, 然后在t_sample时间内
 * 均匀采样n_samples个点(取var_idx分量)
 *
 * 知识点: 瞬态消除, 渐近行为采样
 * 复杂度: O((n_trans + n_samp) * dim)
 */
double* sample_asymptotic(flow_rhs_t f, const double *x0, int dim,
                           const double *p, double t_transient,
                           double t_sample, double dt,
                           int n_samples, int var_idx)
{
    int n_trans = (int)(t_transient / dt);
    int n_samp  = (int)(t_sample / dt);
    int skip = n_samp / n_samples;
    if (skip < 1) skip = 1;

    double *x = (double*)chaos_calloc(dim, sizeof(double));
    double *x_next = (double*)chaos_calloc(dim, sizeof(double));
    double *result = (double*)chaos_calloc(n_samples, sizeof(double));
    int i, k;

    for (k = 0; k < dim; k++) x[k] = x0[k];

    /* 瞬态消除 */
    double t = 0.0;
    for (i = 0; i < n_trans; i++) {
        rk4_step(x_next, f, x, t, dt, dim, p);
        for (k = 0; k < dim; k++) x[k] = x_next[k];
        t += dt;
    }

    /* 渐近采样 */
    int count = 0;
    for (i = 0; i < n_samp && count < n_samples; i++) {
        rk4_step(x_next, f, x, t, dt, dim, p);
        for (k = 0; k < dim; k++) x[k] = x_next[k];
        if (i % skip == 0)
            result[count++] = x[var_idx];
        t += dt;
    }

    free(x); free(x_next);
    return result;
}
