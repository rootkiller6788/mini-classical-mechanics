/**
 * example_logistic.c - Logistic map analysis demo
 * May (1976) Nature 261:459
 */
#include "chaos.h"
#include "chaos_maps.h"
#include "chaos_lyapunov.h"
#include "chaos_bifurcation.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("============================================================\n");
    printf("  Logistic Map Analysis  (May 1976)\n");
    printf("============================================================\n\n");

    printf("[1] Asymptotic behavior at various r:\n\n");
    double r_vals[] = {2.5, 3.2, 3.5, 3.9, 4.0};
    int n_r = sizeof(r_vals) / sizeof(r_vals[0]);
    int ri;

    for (ri = 0; ri < n_r; ri++) {
        double r = r_vals[ri];
        double x = 0.5;
        int i;
        for (i = 0; i < 500; i++) x = logistic_map(x, r);
        printf("    r=%.2f: ", r);
        int j;
        for (j = 0; j < 8; j++) {
            x = logistic_map(x, r);
            printf("%.4f ", x);
        }
        printf("\n");
    }

    printf("\n[2] Lyapunov exponents vs r:\n\n");
    printf("    r        lambda       Type\n");
    printf("    ------   ----------   ----\n");
    double r_check[] = {2.0, 2.5, 3.0, 3.57, 3.8, 3.9, 4.0};
    int nc = sizeof(r_check) / sizeof(r_check[0]);
    int si;

    for (si = 0; si < nc; si++) {
        double lam = logistic_lyapunov(r_check[si], 5000, 500);
        const char *type;
        if (lam < -0.01) type = "stable fixed pt";
        else if (lam < 0.01) type = "bifurcation/pt";
        else type = "chaotic";
        printf("    %.4f    %+10.6f    %s\n", r_check[si], lam, type);
    }

    printf("\n[3] Period-3 window detection (r ~ 3.828):\n");
    int found = 0;
    double r;
    for (r = 3.82; r <= 3.84; r += 0.0001) {
        int p = find_period_at_param(logistic_map, r, 0.5, 1000, 2000, 1e-5);
        if (p == 3 && !found) {
            printf("    Period-3 at r = %.6f  (Li-Yorke: => chaos!)\n", r);
            found = 1;
            break;
        }
    }

    printf("\n[4] Feigenbaum delta estimation:\n");
    double r_bif[8];
    int nf = find_period_doubling_points(r_bif, 5, 2.8, 0.005, 1e-4);
    int i;
    for (i = 0; i < nf; i++)
        printf("    r_%d = %.8f (period %d)\n", i+1, r_bif[i], 1 << (i+1));

    if (nf >= 3) {
        double deltas[6];
        double de = estimate_feigenbaum_delta(r_bif, nf, deltas);
        printf("    Estimated delta = %.6f (true = %.6f)\n", de, FEIGENBAUM_DELTA);
    }

    printf("\nDone.\n");
    return 0;
}
