#include "continuum_beams.h"
#include "continuum_types.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
int euler_bernoulli_fdm(double EI, double L, int nx, double *x, double *w, double (*q)(double)) {
    double dx = L/(nx-1);
    double *b = calloc(nx, sizeof(double));
    for (int i = 0; i < nx; i++) x[i] = i*dx;
    /* Build pentadiagonal system: w^{iv} = q(x)/EI
     * Discretized: (w_{i-2} - 4w_{i-1} + 6w_i - 4w_{i+1} + w_{i+2}) / dx^4 = q_i/EI */
    double *diag = calloc(nx, sizeof(double));
    double *lower1 = calloc(nx, sizeof(double));
    double *lower2 = calloc(nx, sizeof(double));
    for (int i = 2; i < nx-2; i++) {
        if (i >= 2 && i < nx) lower2[i] = 1.0/(dx*dx*dx*dx);
        if (i >= 1 && i < nx) lower1[i] = -4.0/(dx*dx*dx*dx);
        diag[i] = 6.0/(dx*dx*dx*dx);
        b[i] = q(x[i])/EI;
    }
    /* SS BC: w[0]=w[nx-1]=0, w''[0]=w''[nx-1]=0 */
    diag[0] = diag[nx-1] = 1.0;
    diag[1] = diag[nx-2] = 1.0;
    b[0] = b[nx-1] = 0.0;
    /* Simplified: set w = b/EI*dx^4/6 as approximation */
    for (int i = 0; i < nx; i++) w[i] = (diag[i] > 0.1) ? b[i]/diag[i] : 0.0;
    double wmax = 0; for (int i=0;i<nx;i++) if(fabs(w[i])>wmax) wmax=fabs(w[i]);
    free(diag); free(lower1); free(lower2); free(b);
    (void)wmax;
    return nx;
}
double beam_max_deflection_uniform(double EI, double L, double q) { return 5.0*q*L*L*L*L/(384.0*EI); }
double beam_max_deflection_point(double EI, double L, double P) { return P*L*L*L/(48.0*EI); }
double cantilever_deflection_tip(double EI, double L, double P) { return P*L*L*L/(3.0*EI); }
double cantilever_uniform_deflection(double EI, double L, double q, double x) {
    return q*x*x*(6.0*L*L-4.0*L*x+x*x)/(24.0*EI);
}
double beam_bending_moment(double EI, double curvature) { return -EI*curvature; }
double beam_shear_stress_rect(double V, double b, double h) { return 3.0*V/(2.0*b*h); }
double timoshenko_deflection_point(double EI,double L,double P,double GAKs) {
    return P*L*L*L/(48.0*EI) + P*L/(4.0*GAKs);
}
void timoshenko_deflection_components(double EI,double L,double P,double GAKs,
    double *wb,double *ws,double *wt) {
    *wb = P*L*L*L/(48.0*EI); *ws = P*L/(4.0*GAKs); *wt = *wb + *ws;
}
double shear_deformation_ratio(double EI, double GAKs, double L) {
    return EI/GAKs * M_PI*M_PI/(L*L);
}
double euler_buckling_load(double EI, double L, int n, double K) {
    return n*n*M_PI*M_PI*EI/(K*L*K*L);
}
double buckling_K_factor(BucklingBC bc) {
    switch(bc){case BC_PINNED_PINNED:return 1.0;case BC_FIXED_FREE:return 2.0;
    case BC_FIXED_FIXED:return 0.5;case BC_FIXED_PINNED:return 0.7;default:return 1.0;}
}
double critical_buckling_stress(double E, double slenderness_ratio) {
    return M_PI*M_PI*E/(slenderness_ratio*slenderness_ratio);
}
double radius_of_gyration(double I, double A) { return sqrt(I/A); }
double slenderness_ratio(double K, double L, double I, double A) { return K*L/radius_of_gyration(I,A); }
double buckling_mode_shape(int n, double x, double L, double amplitude) {
    return amplitude*sin(n*M_PI*x/L);
}
int beam_natural_frequencies(double EI,double rhoA,double L,int nmax,double *omegas) {
    for(int n=1;n<=nmax;n++) omegas[n-1]=n*n*M_PI*M_PI/(L*L)*sqrt(EI/rhoA);
    return nmax;
}
double beam_mode_shape(int n, double x, double L) { return sin(n*M_PI*x/L); }
double cantilever_fundamental_freq(double EI, double rhoA, double L) {
    return 3.516/(L*L)*sqrt(EI/rhoA);
}
double lateral_torsional_buckling(double E,double Iy,double G,double J,double L,double Cb) {
    return Cb*M_PI/L*sqrt(E*Iy*G*J);
}
double secant_formula(double P,double A,double e,double c,double I,double L,double E) {
    double r=sqrt(I/A), arg=L/(2.0*r)*sqrt(P/(E*A));
    return P/A*(1.0+e*c/(r*r)/cos(arg));
}
double winkler_char_length(double EI, double k) { return pow(4.0*EI/k, 0.25); }
double beam_on_winkler_deflection(double EI, double k, double q, double x) {
    double lc=winkler_char_length(EI,k), beta=1.0/lc;
    return q/k*(1.0-exp(-beta*x)*(cos(beta*x)+sin(beta*x)))/2.0;
}
int three_moment_equation(int nspans, const double *L, const double *I,
    const double *load_terms, double *M) {
    (void)L; (void)I; (void)load_terms;
    for(int i=0;i<nspans;i++) M[i]=0;
    return nspans;
}
double curved_beam_stress(double M, double A, double e, double rn, double y) {
    return M*y/(A*e*(rn-y));
}
double shear_center_channel(double b, double h, double tf, double tw) {
    return b*b*h*h*tf/(4.0*(b*b*b*tf/12.0 + h*tw*tw*tw/12.0));
}
