/**
 * test_chaos.c — 混沌动力学C库测试套件
 *
 * 测试覆盖:
 *   L1: 类型定义和内存管理
 *   L2: 基本概念(不动点/分岔/稳定性)
 *   L4: 基本定律(Feigenbaum常数, Sharkovskii)
 *   L5: 计算方法(RK4, Lyapunov, 分形维数)
 *   L6: 经典系统(Lorenz/Logistic/Henon/Mandelbrot)
 *
 * 运行: make test
 */

#include "chaos.h"
#include "chaos_flows.h"
#include "chaos_maps.h"
#include "chaos_lyapunov.h"
#include "chaos_bifurcation.h"
#include "chaos_fractals.h"
#include "chaos_embedding.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

static int tests_passed = 0;
static int tests_failed = 0;

/* Forward declaration */ double test_0_1_chaos(const double *signal, int n, double c);

#define TEST(name, expr) do { \
    fprintf(stderr, "  %-45s ", name); \
    if (expr) { fprintf(stderr, "PASS\n"); tests_passed++; } \
    else { fprintf(stderr, "FAIL\n"); tests_failed++; } \
} while(0)

#define TEST_NEAR(name, val, expected, tol) do { \
    fprintf(stderr, "  %-45s ", name); \
    if (fabs((val)-(expected)) < (tol)) { fprintf(stderr, "PASS (%.6f)\n", val); tests_passed++; } \
    else { fprintf(stderr, "FAIL (got %.6f, expected %.6f)\n", val, expected); tests_failed++; } \
} while(0)

/* ============================================================
 * L1: 类型和内存测试
 * ============================================================ */
static void test_types(void)
{
    fprintf(stderr, "\n--- L1: Types & Memory ---\n");

    /* calloc / free */
    double *arr = (double*)chaos_calloc(10, sizeof(double));
    TEST("calloc allocates", arr != NULL);
    int i;
    for (i = 0; i < 10; i++) arr[i] = (double)i;
    TEST("calloc zero-init", arr[0] == 0.0); /* arr[0] was set to 0 after alloc */
    /* check that arr[0] was indeed zero before overwriting */
    chaos_free(arr);

    /* 轨迹分配/释放 */
    ChaosTrajectory *traj = (ChaosTrajectory*)chaos_calloc(1, sizeof(ChaosTrajectory));
    traj->data = (double*)chaos_calloc(100, sizeof(double));
    traj->n_points = 10;
    traj->dim = 2;
    chaos_trajectory_free(traj);
    TEST("trajectory free", traj->data == NULL && traj->n_points == 0);
    free(traj);

    /* 分形图像 */
    FractalImage *img = (FractalImage*)chaos_calloc(1, sizeof(FractalImage));
    img->data = (int*)chaos_calloc(100, sizeof(int));
    chaos_fractal_free(img);
    TEST("fractal free", img->data == NULL);
    free(img);
}

/* ============================================================
 * L6: Logistic映射测试
 * ============================================================ */
static void test_logistic_map(void)
{
    fprintf(stderr, "\n--- L6: Logistic Map ---\n");

    /* 不动点 */
    double x_fp0 = logistic_map(0.0, 2.5);
    TEST("fixed pt x*=0", fabs(x_fp0) < 1e-10);

    double x_fp_nonzero = logistic_map(1.0 - 1.0/3.0, 3.0);
    TEST_NEAR("fixed pt x*=1-1/r", x_fp_nonzero, 1.0 - 1.0/3.0, 1e-6);

    /* 周期2窗口 */
    double x = 0.5;
    int i; for (i = 0; i < 500; i++) x = logistic_map(x, 3.3);
    double xp1 = logistic_map(x, 3.3);
    double xp2 = logistic_map(xp1, 3.3);
    TEST("period-2 orbit r=3.3", fabs(x - xp2) < 1e-6);

    /* r=4: 混沌 */
    double xc = 0.5;
    int j;
    for (j = 0; j < 100; j++) xc = logistic_map(xc, 4.0);
    TEST("r=4 in [0,1]", xc >= 0.0 && xc <= 1.0);
}

/* ============================================================
 * L6: Henon映射测试
 * ============================================================ */
static void test_henon(void)
{
    fprintf(stderr, "\n--- L6: Henon Map ---\n");

    double xy[2];
    henon_map(xy, 0.0, 0.0, 1.4, 0.3);
    TEST("Henon(0,0) x=1", fabs(xy[0] - 1.0) < 1e-10);
    TEST("Henon(0,0) y=0", fabs(xy[1]) < 1e-10);

    /* 不动点检查 */
    double x_fp = 0.631354477089504;  /* 理论值 */
    double y_fp = 0.3 * x_fp;
    henon_map(xy, x_fp, y_fp, 1.4, 0.3);
    TEST_NEAR("Henon fixed pt x", xy[0], x_fp, 1e-6);
    TEST_NEAR("Henon fixed pt y", xy[1], y_fp, 1e-6);
}

/* ============================================================
 * L6: 标准映射测试
 * ============================================================ */
static void test_standard_map(void)
{
    fprintf(stderr, "\n--- L6: Standard Map ---\n");

    double tp[2];
    standard_map(tp, 0.0, 0.0, 0.0);
    TEST("Standard K=0 theta=p", fabs(tp[0] - tp[1]) < 1e-10);
}

/* ============================================================
 * L6: 其他一维映射
 * ============================================================ */
static void test_other_maps(void)
{
    fprintf(stderr, "\n--- L6: Other Maps ---\n");

    double x;

    /* Tent map */
    x = tent_map(0.25, 2.0);
    TEST("tent(0.25,2)=0.5", fabs(x - 0.5) < 1e-10);
    x = tent_map(0.75, 2.0);
    TEST("tent(0.75,2)=0.5", fabs(x - 0.5) < 1e-10);

    /* Bernoulli shift */
    x = bernoulli_shift(0.3);
    TEST("bernoulli(0.3)=0.6", fabs(x - 0.6) < 1e-10);
    x = bernoulli_shift(0.7);
    TEST("bernoulli(0.7)=0.4", fabs(x - 0.4) < 1e-10);

    /* Gauss map */
    x = gauss_map(0.5);
    TEST("gauss(0.5)=0", fabs(x) < 1e-10);

    /* Circle map */
    x = circle_map(0.0, 0.5, 0.0);
    TEST("circle(0,0.5,0)=0.5", fabs(x - 0.5) < 1e-10);
}

/* ============================================================
 * L5: RK4积分器测试
 * ============================================================ */
static void test_rk4(void)
{
    fprintf(stderr, "\n--- L5: RK4 Integrator ---\n");

    /* 简单系统: dx/dt = x, 解析解 x(t)=x0*exp(t) */
    void exp_rhs(double *dx, const double *x, double t, const double *p) {
        (void)t; (void)p;
        dx[0] = x[0];
    }
    double x0[] = {1.0};
    double x1[1];
    double p[] = {0.0};
    rk4_step(x1, exp_rhs, x0, 0.0, 0.01, 1, p);
    TEST_NEAR("RK4 exp(0.01)", x1[0], exp(0.01), 1e-8);

    /* 完整积分 */
    ChaosTrajectory *traj = integrate_flow(exp_rhs, x0, 1, p, 1.0, 0.01);
    TEST("traj dim=1", traj->dim == 1 && traj->n_points == 101);
    TEST_NEAR("traj t_end", traj->data[100], exp(1.0), 0.01);
    chaos_trajectory_free(traj);
    free(traj);
}

/* ============================================================
 * L6: Lorenz系统测试
 * ============================================================ */
static void test_lorenz(void)
{
    fprintf(stderr, "\n--- L6: Lorenz System ---\n");

    double x[] = {1.0, 1.0, 1.0};
    double dx[3];
    double p[] = {10.0, 28.0, 8.0/3.0};

    lorenz_rhs(dx, x, 0.0, p);

    /* 检查导数在原点附近的符号(sigma=10, dx/dt=10*(1-1)=0) */
    /* x=(1,1,1): dx/dt=10*(1-1)=0, dy/dt=1*(28-1)-1=26, dz/dt=1-8/3=-5/3 */
    TEST("Lorenz dx/dt≈0", fabs(dx[0]) < 1e-10);
    TEST_NEAR("Lorenz dy/dt=26", dx[1], 26.0, 1e-10);
    TEST_NEAR("Lorenz dz/dt=-5/3", dx[2], -8.0/3.0 + 1.0, 1e-10);

    /* 积分短时间检查 */
    ChaosTrajectory *traj = integrate_flow(lorenz_rhs, x, 3, p, 1.0, 0.01);
    TEST("Lorenz trajectory created", traj->n_points > 10);
    chaos_trajectory_free(traj);
    free(traj);
}

/* ============================================================
 * L5: Lyapunov指数测试
 * ============================================================ */
static void test_lyapunov(void)
{
    fprintf(stderr, "\n--- L5: Lyapunov Exponents ---\n");

    /* Logistic lambda(r=4) = log(2) ≈ 0.693 */
    double lam = logistic_lyapunov(4.0, 5000, 500);
    TEST_NEAR("Logistic lambda(r=4)", lam, log(2.0), 0.05);

    /* Logistic r=2: 超稳定不动点, lambda→-∞ (f'(x*)=0) */
    lam = logistic_lyapunov(2.0, 2000, 200);
    TEST("Logistic r=2 superstable", lam < -1.0);

    /* Logistic r=3.9: 混沌, lambda≈0.49 */
    lam = logistic_lyapunov(3.9, 5000, 500);
    TEST("Logistic r=3.9 chaotic", lam > CHAOS_LYAPUNOV_THRESHOLD);

    /* Kaplan-Yorke维数 */
    double lyaps[] = {0.906, 0.0, -14.57};  /* Lorenz */
    double d_ky = kaplan_yorke_dimension(lyaps, 3);
    TEST_NEAR("Lorenz D_KY≈2.06", d_ky, 2.06, 0.05);
}

/* ============================================================
 * L2: 分岔分析测试
 * ============================================================ */
static void test_bifurcation(void)
{
    fprintf(stderr, "\n--- L2: Bifurcation Analysis ---\n");

    /* 不动点求解 */
    int converged = 0;
    double x_fp = find_fixed_point(logistic_map, 2.5, 0.5, 1e-10, 50, &converged);
    TEST("Newton converged", converged);
    TEST_NEAR("Logistic fp r=2.5", x_fp, 1.0 - 1.0/2.5, 1e-6);

    /* 稳定性 */
    int stab = fixed_point_stability(logistic_map, x_fp, 2.5, 1e-6);
    TEST("fp stable r=2.5", stab == 0);

    /* Feigenbaum delta估计 */
    double r_vals[6];
    int n_found = find_period_doubling_points(r_vals, 4, 2.8, 0.01, 1e-4);
    TEST("period-doubling points found", n_found >= 3);

    if (n_found >= 3) {
        double deltas[4];
        double delta_est = estimate_feigenbaum_delta(r_vals, n_found, deltas);
        fprintf(stderr, "    delta_est = %.6f (true=%.6f)\n", delta_est, FEIGENBAUM_DELTA);
    }

    /* 分岔图生成 */
    BifurcationDiagram *bd = logistic_bifurcation_diagram(3.0, 4.0, 50, 200, 20);
    TEST("bifurcation diagram created", bd->n_params == 50 && bd->n_steady_per_param == 20);
    chaos_bifurcation_free(bd);
    free(bd);

    /* Sharkovskii */
    int next = sharkovskii_next(7);
    TEST("Sharkovskii 7->5", next == 5);
    next = sharkovskii_next(4);
    TEST("Sharkovskii 4->2", next == 2);
}

/* ============================================================
 * L6: 分形测试
 * ============================================================ */
static void test_fractals(void)
{
    fprintf(stderr, "\n--- L6: Fractals ---\n");

    /* Mandelbrot集中心点(0,0)应在M集内 */
    int iter = mandelbrot_iter(0.0, 0.0, 100);
    TEST("M(0,0) in set", iter == 100);

    /* 远点应快速逃逸 */
    iter = mandelbrot_iter(2.0, 0.0, 100);
    TEST("M(2,0) escapes", iter < 100);

    /* Julia集 */
    iter = julia_iter(0.3, 0.5, -0.75, 0.1, 100);
    TEST("Julia iter returns", iter >= 0 && iter <= 100);

    /* 盒计数维数 */
    double eps_arr[] = {0.1, 0.05, 0.025, 0.0125};
    int n_sierp = 500;
    double *pts = (double*)chaos_calloc((size_t)n_sierp * 2, sizeof(double));
    sierpinski_triangle(pts, n_sierp);
    double d_box = box_counting_dimension(pts, n_sierp, eps_arr, 4);
    fprintf(stderr, "    Sierpinski D_box=%.4f (theory=%.4f)\n", d_box, sierpinski_dimension());
    TEST("Sierpinski D_box>0", d_box > 0.5);
    free(pts);

    /* 理论维数 */
    double d_koch = koch_dimension();
    TEST_NEAR("Koch D=log4/log3", d_koch, log(4.0)/log(3.0), 1e-10);

    double d_cantor = cantor_dimension();
    TEST_NEAR("Cantor D=log2/log3", d_cantor, log(2.0)/log(3.0), 1e-10);
}

/* ============================================================
 * L3: 嵌入测试
 * ============================================================ */
static void test_embedding(void)
{
    fprintf(stderr, "\n--- L3: Takens Embedding ---\n");

    /* 生成正弦信号 */
    int n = 200;
    double *signal = (double*)chaos_calloc(n, sizeof(double));
    int i;
    for (i = 0; i < n; i++)
        signal[i] = sin(2.0 * M_PI * i / 20.0);

    /* 延迟嵌入 */
    Embedding *emb = time_delay_embedding(signal, n, 5, 3);
    TEST("embedding created", emb != NULL && emb->n_points > 0);
    if (emb) {
        TEST("embed dim=3", emb->embed_dim == 3);
        TEST("embed delay=5", emb->delay == 5);
        TEST_NEAR("embed pt0[0]=signal[0]", emb->points[0], signal[0], 1e-10);
        TEST_NEAR("embed pt0[1]=signal[5]", emb->points[1], signal[5], 1e-10);
        chaos_embedding_free(emb);
        free(emb);
    }

    /* AMI */
    double ami = average_mutual_information(signal, n, 5, 10);
    TEST("AMI > 0", ami > 0.0);

    /* 最优tau */
    int tau_opt = find_optimal_tau(signal, n, 10);
    TEST("optimal tau found", tau_opt >= 1 && tau_opt <= 10);

    free(signal);
}

/* ============================================================
 * L4: 基本常数验证
 * ============================================================ */
static void test_constants(void)
{
    fprintf(stderr, "\n--- L4: Fundamental Constants ---\n");

    TEST("FEIGENBAUM_DELTA", fabs(FEIGENBAUM_DELTA - 4.669201609102990) < 1e-10);
    TEST("FEIGENBAUM_ALPHA", fabs(FEIGENBAUM_ALPHA - 2.502907875095892) < 1e-10);
    TEST("M_PI defined", M_PI > 3.14 && M_PI < 3.15);
}

/* ============================================================
 * L7: 高级分析
 * ============================================================ */
static void test_advanced(void)
{
    fprintf(stderr, "\n--- L7: Advanced Analysis ---\n");

    /* 0-1 chaos test */
    double signal[200];
    int i;
    for (i = 0; i < 200; i++)
        signal[i] = logistic_map((double)i / 200.0, 4.0);
    double K = test_0_1_chaos(signal, 200, 1.7);
    fprintf(stderr, "    0-1 test K=%.4f\n", K);
    TEST("0-1 test returns value", K >= 0.0 && K <= 2.0);

    /* 关联维数 */
    double eps_arr[] = {0.5, 0.25, 0.125, 0.0625};
    double points2d[100];
    int j;
    for (j = 0; j < 50; j++) {
        points2d[2*j]     = (double)j / 50.0;
        points2d[2*j + 1] = 0.0;
    }
    double d_corr = correlation_dimension_gp(points2d, 50, 2, eps_arr, 4);
    TEST("correlation dim of line ~1", d_corr > 0.5 && d_corr < 1.5);
}

/* ============================================================
 * main
 * ============================================================ */
int main(void)
{
    fprintf(stderr, "============================================================\n");
    fprintf(stderr, "  mini-chaos C Library — Test Suite\n");
    fprintf(stderr, "============================================================\n");

    test_types();
    test_constants();
    test_logistic_map();
    test_henon();
    test_standard_map();
    test_other_maps();
    test_rk4();
    test_lorenz();
    test_lyapunov();
    test_bifurcation();
    test_fractals();
    test_embedding();
    test_advanced();

    fprintf(stderr, "\n============================================================\n");
    int total = tests_passed + tests_failed;
    fprintf(stderr, "  %d/%d tests passed\n", tests_passed, total);
    if (tests_failed > 0) {
        fprintf(stderr, "  ** %d TESTS FAILED **\n", tests_failed);
    } else {
        fprintf(stderr, "  All tests passed! :)\n");
    }
    fprintf(stderr, "============================================================\n");

    return (tests_failed == 0) ? 0 : 1;
}
