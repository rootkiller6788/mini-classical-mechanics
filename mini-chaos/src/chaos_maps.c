/**
 * chaos_maps.c — 离散时间混沌映射实现
 *
 * 涵盖:
 *   L1: Logistic/Tent/Sine/Gauss映射定义
 *   L6: Henon/Standard/Ikeda/Circle/Arnold/Baker等经典映射
 *   L5: 迭代工具(瞬态消除+轨迹收集)
 *
 * 参考: May (1976), Henon (1976), Chirikov (1979), Strogatz Ch.10
 */

#include "chaos_maps.h"
#include <stdio.h>
#include <math.h>

/* ============================================================
 * L6: 一维离散映射
 * ============================================================ */

/**
 * Logistic映射: x_{n+1} = r*x_n*(1 - x_n)
 *
 * 定义域: x in [0,1], r in [0,4]
 *
 * 不动点(解析):
 *   x* = 0 (对所有r)
 *   x* = 1 - 1/r (r > 1)
 *
 * 稳定性(乘子): mu = r*(1-2x*)
 *   x*=0: mu=r     -> 稳定当 r<1
 *   x*=1-1/r: mu=2-r -> 稳定当 1<r<3
 *
 * 周期倍增点: r1=3, r2=1+sqrt(6)≈3.449, r3≈3.544, ...
 * Feigenbaum点: r_inf≈3.569946
 *
 * 完全混沌(r=4): 不变测度 rho(x)=1/(pi*sqrt(x*(1-x)))
 * Lyapunov指数=lambda=log(2)≈0.693
 *
 * 定理(Li-Yorke 1975): 周期3存在 => 所有周期存在 => 混沌
 * 周期3窗口: r≈3.828
 *
 * 参考: May (1976) Nature, Strogatz §10.1-10.4
 * MIT 8.012 / Caltech Ph 106
 */
double logistic_map(double x, double r)
{
    return r * x * (1.0 - x);
}

/**
 * 立方映射: x_{n+1} = r*x_n - x_n^3
 *
 * 不动点: x*=0 和 x*=±sqrt(r-1) (r>1)
 * 叉形分岔在r=1
 *
 * 知识点: 对称性破缺, 双稳态, 超临界叉形分岔
 * 参考: Strogatz §3.4
 */
double cubic_map(double x, double r)
{
    return r * x - x * x * x;
}

/**
 * Sine映射: x_{n+1} = r*sin(pi*x_n)
 *
 * 与Logistic映射拓扑共轭(通过h(x)=sin^2(pi*x/2))
 * 分岔序列相同, Feigenbaum常数相同
 *
 * 知识点: 拓扑共轭, 普适性
 * 参考: Strogatz §10.5
 */
double sine_map(double x, double r)
{
    return r * sin(M_PI * x);
}

/**
 * Tent映射: x_{n+1} = mu*min(x, 1-x)
 *
 * mu=2: Lyapunov=log(2), 完全混沌
 * mu<2: 有稳定周期轨道
 * mu>2: 逃逸到负无穷
 *
 * 拓扑共轭于Bernoulli移位(通过共轭函数)
 * 不变测度(μ=2): 均匀分布 rho(x)=1
 *
 * 知识点: 分段线性, Markov分割, 符号动力学, 共轭
 * 参考: Devaney "Chaotic Dynamical Systems" §1.4
 */
double tent_map(double x, double mu)
{
    return mu * (x < 0.5 ? x : 1.0 - x);
}

/**
 * Gauss映射(连分数映射): x_{n+1} = 1/x_n mod 1
 *
 * 即: x_{n+1} = 1/x_n - floor(1/x_n), x≠0
 *
 * 不变测度(Gauss测度): rho(x) = 1/((1+x)*log(2))
 * 各态历经: int_0^1 f(x)*rho(x)dx = lim (1/n)*sum f(x_k)
 *
 * 知识点: 遍历理论, 连分数, 不变测度
 * 参考: Khinchin "Continued Fractions" (1964)
 * Princeton PHY 505: 遍历性
 */
double gauss_map(double x)
{
    if (fabs(x) < 1e-15) return 0.0;
    double inv = 1.0 / x;
    return inv - floor(inv);
}

/**
 * Bernoulli移位: x_{n+1} = 2*x_n mod 1
 *
 * 二进制展开: x = 0.b1b2b3..., 移位操作相当于左移并丢弃整数部分
 * 符号动力学: 序列{b_k}完全决定轨道
 *
 * 知识点: 决定论混沌=符号序列的伪随机性
 * Lyapunov指数=lambda=log(2)
 *
 * 定理: Bernoulli移位是混沌的(拓扑传递, 周期轨道稠密, 敏感依赖)
 * 参考: Devaney §1.3
 */
double bernoulli_shift(double x)
{
    double y = 2.0 * x;
    return y - floor(y);
}

/**
 * 透镜映射: x_{n+1} = (a*x_n) mod 1, a>1
 *
 * Lyapunov指数 = log(a) > 0 => 混沌
 *
 * 知识点: 扩张映射, 一致双曲性
 */
double lens_map(double x, double a)
{
    double y = a * x;
    return y - floor(y);
}

/* ============================================================
 * L6: 二维离散映射
 * ============================================================ */

/**
 * Henon映射: x_{n+1} = 1 - a*x_n^2 + y_n
 *            y_{n+1} = b*x_n
 *
 * 经典参数: a=1.4, b=0.3
 *
 * Jacobian: det(J) = -b => 面积收缩因子|b|
 * 不动点(解析): x* = [-(1-b) ± sqrt((1-b)^2+4a)]/(2a), y*=b*x*
 *
 * Smale马蹄: 在适当参数下存在横截面马蹄 -> 混沌严格证明
 *
 * 知识点: 奇怪吸引子, Smale马蹄, 同宿缠结
 * 参考: Henon (1976), Strogatz §12.2
 * Berkeley PHYS 231 / Cambridge Part III
 */
void henon_map(double *xy_out, double x, double y, double a, double b)
{
    xy_out[0] = 1.0 - a * x * x + y;
    xy_out[1] = b * x;
}

/**
 * Lozi映射: x_{n+1} = 1 - a*|x_n| + y_n
 *           y_{n+1} = b*x_n
 *
 * 分段线性 => 可严格证明混沌(Misiurewicz 1980)
 *
 * 知识点: 分段线性奇怪吸引子, 严格混沌证明
 */
void lozi_map(double *xy_out, double x, double y, double a, double b)
{
    xy_out[0] = 1.0 - a * fabs(x) + y;
    xy_out[1] = b * x;
}

/**
 * 标准映射(Chirikov-Taylor映射)
 *
 * p_{n+1} = p_n + K*sin(theta_n)
 * theta_{n+1} = theta_n + p_{n+1}  (mod 2*pi)
 *
 * Hamiltonian: H = p^2/2 + K*cos(theta)*sum delta(t-n)
 *
 * K<<1: 可积(不变环面)
 * K≈1: 最后一个KAM环面破坏(临界K_c≈0.9716)
 * K>>1: 全局混沌
 *
 * 知识点: Hamiltonian混沌, KAM定理, 共振重叠准则(Chirikov判据)
 * 参考: Chirikov (1979), Lichtenberg & Lieberman (1992)
 * Goldstein §11.4 / Caltech Ph 205
 */
void standard_map(double *thetap_out, double theta, double p, double K)
{
    double p_new = p + K * sin(theta);
    double theta_new = theta + p_new;

    /* 映射到[0, 2*pi) */
    theta_new = fmod(theta_new, M_2PI);
    if (theta_new < 0.0) theta_new += M_2PI;
    p_new = fmod(p_new + M_PI, M_2PI) - M_PI;

    thetap_out[0] = theta_new;
    thetap_out[1] = p_new;
}

/**
 * Ikeda映射(光学环形腔)
 *
 * 复形式: E_{n+1} = A + B*E_n*exp(i*(phi - p/(1+|E|^2)))
 *
 * 实形式:
 *   x_{n+1} = A + B*(x_n*cos(tau_n) - y_n*sin(tau_n))
 *   y_{n+1} = B*(x_n*sin(tau_n) + y_n*cos(tau_n))
 *   其中 tau_n = phi - p/(1 + x_n^2 + y_n^2)
 *
 * 物理背景: 光学双稳态腔, 延迟反馈
 * 知识点: 光学混沌, 延迟微分方程近似
 * 参考: Ikeda (1979), Ikeda, Daido, Akimoto (1980)
 */
void ikeda_map(double *xy_out, double x, double y,
               double A, double B, double phi, double p)
{
    double r2 = x * x + y * y;
    double tau = phi - p / (1.0 + r2);
    double cs = cos(tau), sn = sin(tau);
    xy_out[0] = A + B * (x * cs - y * sn);
    xy_out[1] = B * (x * sn + y * cs);
}

/**
 * 圆映射: theta_{n+1} = theta_n + Omega - (K/(2*pi))*sin(2*pi*theta_n)  (mod 1)
 *
 * 卷绕数: W = lim (theta_n - theta_0)/n
 *
 * 锁频: 在有理W处形成Arnold舌头
 * 临界K=1: 魔鬼阶梯具有完整测度
 * K>1: 混沌
 *
 * 知识点: 锁频, Arnold舌头, 魔鬼阶梯, 拟周期->混沌
 * 参考: Arnold (1965), Strogatz §8.6
 * Oxford CMT: 频率同步
 */
double circle_map(double theta, double Omega, double K)
{
    double t = theta + Omega - (K / M_2PI) * sin(M_2PI * theta);
    t = fmod(t, 1.0);
    if (t < 0.0) t += 1.0;
    return t;
}

/**
 * Arnold猫映射: x_{n+1} = 2*x_n + y_n (mod 1)
 *               y_{n+1} = x_n + y_n (mod 1)
 *
 * 矩阵表示: |2 1|
 *          |1 1|  det=1 (保面积)
 *
 * 特征值: lambda1,2 = (3 ± sqrt(5))/2 ≈ 2.618, 0.382
 * Lyapunov指数: lambda1=log(2.618)≈0.962, lambda2=log(0.382)≈-0.962
 *
 * 知识点: Anosov微分同胚, 遍历性, 混合性, 双曲环面自同构
 * 参考: Arnold & Avez (1968)
 */
void arnold_cat_map(double *xy_out, double x, double y)
{
    double xn = 2.0 * x + y;
    double yn = x + y;
    xy_out[0] = xn - floor(xn);
    xy_out[1] = yn - floor(yn);
}

/**
 * Baker映射(面包师变换): x_{n+1} = 2*x_n mod 1
 *                        y_{n+1} = 0.5*y_n + 0.5*floor(2*x_n)
 *
 * 相当于面的拉伸和折叠(面包师揉面)
 *
 * 知识点: Bernoulli系统的二维嵌入, 拉伸-折叠范式
 * 定量混沌: Lyapunov lambda1=log(2), lambda2=-log(2)
 * KS熵 = log(2)
 */
void bakers_map(double *xy_out, double x, double y)
{
    double sx = (x < 0.5) ? 0.0 : 1.0;
    xy_out[0] = 2.0 * x - floor(2.0 * x);
    xy_out[1] = 0.5 * y + 0.5 * sx;
}

/**
 * Gingerbreadman映射: x_{n+1} = 1 - y_n + |x_n|
 *                     y_{n+1} = x_n
 *
 * 保面积: Jacobian行列式=1
 * 特征: 六边形周期岛, 混沌海
 *
 * 知识点: 保面积映射中的混合结构, KAM岛屿
 */
void gingerbreadman_map(double *xy_out, double x, double y)
{
    xy_out[0] = 1.0 - y + fabs(x);
    xy_out[1] = x;
}

/* ============================================================
 * L5: 迭代工具
 * ============================================================ */

/**
 * 迭代1D映射+瞬态消除+轨迹收集
 *
 * 知识点: 瞬态消除的重要性 — 避免初始条件带来的非渐进行为
 *
 * @param n_transient  丢弃的前n个点
 * @param n_iter       收集的点数
 * @return 长度为n_iter的数组(caller free)
 */
double* iterate_map1d(map1d_t f, double x0, double param,
                      int n_transient, int n_iter)
{
    double *traj = (double*)chaos_calloc(n_iter, sizeof(double));
    double x = x0;
    int i;

    /* 瞬态消除 */
    for (i = 0; i < n_transient; i++)
        x = f(x, param);

    /* 轨迹收集 */
    for (i = 0; i < n_iter; i++) {
        x = f(x, param);
        traj[i] = x;
    }
    return traj;
}

/**
 * 迭代2D映射+瞬态消除+轨迹收集
 *
 * @param n_transient  丢弃的迭代次数
 * @param n_iter       收集的迭代次数
 * @return 长度为2*n_iter的数组[x0,y0,x1,y1,...] (caller free)
 */
double* iterate_map2d(map2d_t f, double x0, double y0,
                      double p1, double p2,
                      int n_transient, int n_iter)
{
    double *traj = (double*)chaos_calloc((size_t)n_iter * 2, sizeof(double));
    double x = x0, y = y0;
    double xy[2];
    int i;

    /* 瞬态消除 */
    for (i = 0; i < n_transient; i++) {
        f(xy, x, y, p1, p2);
        x = xy[0]; y = xy[1];
    }

    /* 轨迹收集 */
    for (i = 0; i < n_iter; i++) {
        f(xy, x, y, p1, p2);
        x = xy[0]; y = xy[1];
        traj[2 * i]     = x;
        traj[2 * i + 1] = y;
    }
    return traj;
}
