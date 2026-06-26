/**
 * example_lorenz.c — Lorenz吸引子模拟
 *
 * 演示:
 *   - Lorenz系统积分
 *   - Poincare截面提取
 *   - Lyapunov指数计算
 *
 * 编译: make examples
 * 运行: ./examples_c/example_lorenz
 *
 * 参考: Lorenz (1963) J.Atmos.Sci. 20:130-141
 */

#include "chaos.h"
#include "chaos_flows.h"
#include "chaos_lyapunov.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("============================================================\n");
    printf("  Lorenz Attractor Simulation\n");
    printf("  Lorenz (1963) J.Atmos.Sci. 20:130-141\n");
    printf("============================================================\n\n");

    /* 经典参数: sigma=10, rho=28, beta=8/3 */
    double p[] = {10.0, 28.0, 8.0/3.0};
    double x0[] = {1.0, 1.0, 1.0};

    /* 积分轨迹 */
    printf("[1] Integrating Lorenz system (t_end=50, dt=0.01)...\n");
    ChaosTrajectory *traj = integrate_flow(lorenz_rhs, x0, 3, p, 50.0, 0.01);
    printf("    Trajectory: %d points, dim=%d\n", traj->n_points, traj->dim);

    /* 打印最后几个点 */
    int i, k;
    int n = traj->n_points;
    printf("    Last 3 points:\n");
    for (i = n - 3; i < n && i >= 0; i++) {
        printf("      t≈%.2f: (", i * traj->dt);
        for (k = 0; k < 3; k++)
            printf("% .6f%s", traj->data[i * 3 + k], k < 2 ? ", " : "");
        printf(")\n");
    }

    /* 统计x分量的范围 */
    double x_min = traj->data[0], x_max = traj->data[0];
    double z_min = traj->data[2], z_max = traj->data[2];
    for (i = 1; i < n; i++) {
        if (traj->data[i * 3] < x_min) x_min = traj->data[i * 3];
        if (traj->data[i * 3] > x_max) x_max = traj->data[i * 3];
        if (traj->data[i * 3 + 2] < z_min) z_min = traj->data[i * 3 + 2];
        if (traj->data[i * 3 + 2] > z_max) z_max = traj->data[i * 3 + 2];
    }
    printf("    x range: [%.3f, %.3f]\n", x_min, x_max);
    printf("    z range: [%.3f, %.3f]\n\n", z_min, z_max);

    /* Poincare截面: z = rho-1 = 27 (过Lorenz映射的典型截面) */
    printf("[2] Poincare section at z = 27...\n");
    PoincareSection *ps = poincare_section_extract(
        lorenz_rhs, x0, 3, p, 50.0, 0.01, 2, 27.0);
    printf("    Found %d section crossings\n", ps->n_sections);

    if (ps->n_sections > 0) {
        printf("    First 5 crossing x-coordinates (Lorenz map):\n");
        for (i = 0; i < 5 && i < ps->n_sections; i++)
            printf("      x[%d] = %.6f\n", i, ps->points[i * 3]);
    }
    chaos_poincare_free(ps);
    free(ps);

    /* Lyapunov指数 */
    printf("\n[3] Computing Lyapunov spectrum (Benettin method)...\n");
    printf("    This may take a moment...\n");

    LyapunovSpectrum spect;
    int ret = lyapunov_spectrum_benettin(&spect, lorenz_rhs, x0, 3, p, 100.0, 0.01, 1e-8, 500, 10);
    if (ret == 0) {
        printf("    Lyapunov exponents:\n");
        for (i = 0; i < spect.n_exponents; i++)
            printf("      lambda[%d] = %+.6f\n", i, spect.exponents[i]);
        printf("    KS entropy bound: %.6f\n", spect.sum_positive);
        printf("    Kaplan-Yorke dimension: %.4f\n", spect.kaplan_yorke);
        printf("    (Literature: D_KY ≈ 2.06 for Lorenz)\n");
        chaos_lyapunov_free(&spect);
    } else {
        printf("    Lyapunov computation failed.\n");
    }

    /* 对初始条件的敏感性 */
    printf("\n[4] Sensitivity to initial conditions...\n");
    double xA[] = {1.0, 1.0, 1.0};
    double xB[] = {1.0001, 1.0, 1.0};  /* 微扰 1e-4 */
    double *xA_cur = (double*)chaos_calloc(3, sizeof(double));
    double *xB_cur = (double*)chaos_calloc(3, sizeof(double));
    double *xA_next = (double*)chaos_calloc(3, sizeof(double));
    double *xB_next = (double*)chaos_calloc(3, sizeof(double));
    for (k = 0; k < 3; k++) { xA_cur[k] = xA[k]; xB_cur[k] = xB[k]; }

    for (i = 0; i <= 1000; i++) {
        double t = i * 0.01;
        rk4_step(xA_next, lorenz_rhs, xA_cur, t, 0.01, 3, p);
        rk4_step(xB_next, lorenz_rhs, xB_cur, t, 0.01, 3, p);
        for (k = 0; k < 3; k++) {
            xA_cur[k] = xA_next[k];
            xB_cur[k] = xB_next[k];
        }

        if (i % 200 == 0) {
            double diff = 0.0;
            for (k = 0; k < 3; k++) {
                double d = xA_cur[k] - xB_cur[k];
                diff += d * d;
            }
            diff = sqrt(diff);
            printf("    t=%.1f: |xA-xB| = %.6e\n", i * 0.01, diff);
        }
    }
    free(xA_cur); free(xB_cur); free(xA_next); free(xB_next);

    chaos_trajectory_free(traj);
    free(traj);
    printf("\nDone.\n");
    return 0;
}
