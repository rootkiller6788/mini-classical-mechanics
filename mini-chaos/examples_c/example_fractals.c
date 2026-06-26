/**
 * example_fractals.c — 分形生成演示
 *
 * 演示:
 *   - Mandelbrot集迭代计数
 *   - Julia集图像
 *   - Sierpinski三角
 *   - 盒计数维数
 *
 * 编译: make examples
 * 运行: ./examples_c/example_fractals
 *
 * 参考: Mandelbrot (1983), Peitgen & Richter (1986)
 */

#include "chaos.h"
#include "chaos_fractals.h"
#include <stdio.h>
#include <math.h>

static void print_ascii_mandelbrot(int width, int height)
{
    double x_min = -2.0, x_max = 0.5;
    double y_min = -1.2, y_max = 1.2;
    double dx = (x_max - x_min) / (width - 1);
    double dy = (y_max - y_min) / (height - 1);
    const char *chars = " .:-=+*#%@";
    int py, px;

    for (py = 0; py < height; py++) {
        double cy = y_min + py * dy;
        for (px = 0; px < width; px++) {
            double cx = x_min + px * dx;
            int iter = mandelbrot_iter(cx, cy, 30);
            int idx = (iter * 9) / 30;
            if (idx > 9) idx = 9;
            putchar(chars[idx]);
        }
        putchar('\n');
    }
}

int main(void)
{
    printf("============================================================\n");
    printf("  Fractal Gallery\n");
    printf("  Mandelbrot / Julia / Sierpinski / Newton\n");
    printf("============================================================\n\n");

    /* 1. Mandelbrot集ASCII艺术 */
    printf("[1] Mandelbrot Set (ASCII art, 80x40):\n\n");
    print_ascii_mandelbrot(80, 40);

    /* 2. Mandelbrot集关键点的迭代数 */
    printf("\n[2] Iteration counts at key points:\n\n");
    printf("    Point (c)       Iterations  In Set?\n");
    printf("    --------------  ----------  -------\n");
    double cx[] = {0.0, -0.75, -1.0, -1.5, 0.3, -0.5, -2.0};
    double cy[] = {0.0,  0.1,   0.0,  0.0,  0.5,  0.6,  0.0};
    int n_pts = sizeof(cx) / sizeof(cx[0]);
    int i;

    for (i = 0; i < n_pts; i++) {
        int iter = mandelbrot_iter(cx[i], cy[i], 200);
        printf("    (% .2f, % .2f)      %3d         %s\n",
               cx[i], cy[i], iter,
               iter == 200 ? "YES" : "NO");
    }

    /* 3. Julia集选集 */
    printf("\n[3] Julia sets (iteration count at origin):\n\n");
    double jc[] = {-0.75, 0.1, -0.1, 0.65, -1.25};
    double jcy[] = {0.0, 0.0, 0.8, 0.0, 0.0};
    const char *jname[] = {"San Marco", "Dendrite", "Spiral", "Rabbit", "Airplane"};
    int nj = sizeof(jc) / sizeof(jc[0]);

    for (i = 0; i < nj; i++) {
        int iter = julia_iter(0.0, 0.0, jc[i], jcy[i], 200);
        printf("    %-10s  c=(% .2f,% .2f)  origin iter=%3d  %s\n",
               jname[i], jc[i], jcy[i], iter,
               iter == 200 ? "bounded" : "escaped");
    }

    /* 4. Newton分形 */
    printf("\n[4] Newton fractal for z^3-1=0:\n\n");
    printf("    Roots: 1, exp(2*pi*i/3), exp(4*pi*i/3)\n");
    printf("    Basin classification at test points:\n");

    double nx[] = {0.5, 0.0, -0.8, -0.3,  0.0, 0.2};
    double ny[] = {0.0, 0.5, -0.2, -0.5, -0.5, 0.3};
    int ntest = sizeof(nx) / sizeof(nx[0]);

    for (i = 0; i < ntest; i++) {
        int root = newton_fractal(nx[i], ny[i], 50, 1e-6);
        const char *root_name[] = {"none", "1", "exp(2pi*i/3)", "exp(4pi*i/3)"};
        printf("    z0=(% .2f,% .2f) -> root %d (%s)\n",
               nx[i], ny[i], root,
               root >= 0 && root <= 3 ? root_name[root] : "?");
    }

    /* 5. Sierpinski三角和维数 */
    printf("\n[5] Sierpinski triangle (box-counting dimension):\n\n");
    int n_sierp = 2000;
    double *pts = (double*)chaos_calloc((size_t)n_sierp * 2, sizeof(double));
    sierpinski_triangle(pts, n_sierp);

    double eps_arr[] = {0.2, 0.1, 0.05, 0.025, 0.0125};
    int n_eps = sizeof(eps_arr) / sizeof(eps_arr[0]);
    double d_box = box_counting_dimension(pts, n_sierp, eps_arr, n_eps);

    printf("    Box-counting dimension: %.4f\n", d_box);
    printf("    Theoretical value:      %.4f  (= log(3)/log(2))\n",
           sierpinski_dimension());
    printf("    Points generated: %d\n", n_sierp);

    /* 打印前面几个点的坐标 */
    printf("    First 5 points: ");
    int pi;
    for (pi = 0; pi < 5 && pi < n_sierp; pi++)
        printf("(%.4f,%.4f) ", pts[2*pi], pts[2*pi+1]);
    printf("\n");

    free(pts);

    /* 6. Koch雪花 */
    printf("\n[6] Koch snowflake (3 iterations):\n\n");
    int max_segs = 200;
    double *kx = (double*)chaos_calloc(max_segs, sizeof(double));
    double *ky = (double*)chaos_calloc(max_segs, sizeof(double));
    int n_segs = koch_snowflake(kx, ky, max_segs, 3);

    printf("    Koch D = log(4)/log(3) = %.4f\n", koch_dimension());
    printf("    Segments after 3 iterations: %d\n", n_segs - 1);
    printf("    Perimeter growth: 3 -> 12 -> 48 -> 192 segments\n");

    free(kx); free(ky);

    printf("\nDone.\n");
    return 0;
}
