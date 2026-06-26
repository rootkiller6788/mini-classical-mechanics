#include "continuum_types.h"
#include "continuum_elasticity.h"
#include "continuum_waves.h"
#include "continuum_beams.h"
#include "continuum_failure.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

int main(void) {
    int passed = 0, failed = 0;
#define TEST(name, cond) do { \
    if (cond) { printf("  * %s\n", name); passed++; } \
    else { printf("  X %s\n", name); failed++; } \
} while(0)

    printf("=== mini-continuum C Tests ===\n");

    /* StressTensor */
    StressTensor s = stress_create(100e6, 50e6, 0, 25e6, 0, 0);
    StressInvariants inv = stress_invariants(s);
    TEST("stress invariants I1", fabs(inv.I1 - 1.5e8) < 1e3);
    TEST("stress von Mises > 0", stress_von_mises(s) > 0);
    TEST("stress Tresca > 0", stress_tresca(s) > 0);
    double p[3]; stress_principal(s, p);
    TEST("principal ordering", p[0] >= p[1] && p[1] >= p[2]);

    /* StrainTensor */
    StrainTensor e = strain_from_eng(0.001, 0, 0, 0.002, 0, 0);
    TEST("strain dilatation", fabs(strain_dilatation(e) - 0.001) < 1e-10);

    /* Elasticity */
    double lam, mu;
    lame_constants(200e9, 0.3, &lam, &mu);
    TEST("Lame constants > 0", lam > 0 && mu > 0);

    StressTensor sh = hookes_law_isotropic(e, lam, mu);
    TEST("Hooke returns stress", sh.xx > 0);

    /* Waves */
    double vp = p_wave_speed(lam, mu, 7800);
    double vs = s_wave_speed(mu, 7800);
    TEST("P > S wave speed", vp > vs);

    /* Beams */
    double EI = 1e6, L = 10.0;
    double Pcr = euler_buckling_load(EI, L, 1, 1.0);
    TEST("Euler buckling > 0", Pcr > 0);

    double wmax = beam_max_deflection_uniform(EI, L, 1000.0);
    TEST("beam deflection > 0", wmax > 0);

    /* Failure */
    double sigma_vm = stress_von_mises(s);
    double sf = safety_factor_static(250e6, sigma_vm);
    TEST("safety factor > 0", sf > 0);

    printf("\n%d passed, %d failed\n", passed, failed);
    return failed > 0 ? 1 : 0;
}
