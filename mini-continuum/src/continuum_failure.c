#include "continuum_failure.h"
#include "continuum_types.h"
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
double von_mises_yield_func(StressTensor sigma, double sigma_y) {
    StressInvariants inv=stress_invariants(sigma);
    return sqrt(3.0*inv.J2) - sigma_y;
}
double tresca_yield_func(StressTensor sigma, double sigma_y) {
    double p[3]; stress_principal(sigma,p);
    double mx=fmax(fabs(p[0]-p[1]),fmax(fabs(p[1]-p[2]),fabs(p[2]-p[0])));
    return mx - sigma_y;
}
double rankine_yield_func(StressTensor sigma, double sigma_ut) {
    double p[3]; stress_principal(sigma,p);
    return fmax(p[0],fmax(p[1],p[2])) - sigma_ut;
}
double saint_venant_yield_func(StressTensor sigma, double E, double nu, double eps_ut) {
    double p[3]; stress_principal(sigma,p);
    double e1=(p[0]-nu*(p[1]+p[2]))/E;
    double e2=(p[1]-nu*(p[0]+p[2]))/E;
    double e3=(p[2]-nu*(p[0]+p[1]))/E;
    return fmax(e1,fmax(e2,e3)) - eps_ut;
}
double mohr_coulomb_shear_strength(double c, double sigma_n, double phi) {
    return c + sigma_n*tan(phi);
}
double mohr_coulomb_envelope(StressTensor sigma, double sigma_ut, double sigma_uc) {
    double p[3]; stress_principal(sigma,p);
    return p[0]/sigma_ut - p[2]/sigma_uc;
}
double drucker_prager_yield_func(StressTensor sigma, double alpha, double k) {
    StressInvariants inv=stress_invariants(sigma);
    return sqrt(inv.J2) + alpha*inv.I1 - k;
}
void drucker_prager_params_mc(double c, double phi, double *alpha, double *k) {
    double sp=sin(phi);
    *alpha = 2.0*sp/(sqrt(3.0)*(3.0-sp));
    *k = 6.0*c*cos(phi)/(sqrt(3.0)*(3.0-sp));
}
double hoek_brown(double sigma3, double sigma_ci, double mb, double s, double a) {
    return sigma3 + sigma_ci*pow(mb*sigma3/sigma_ci + s, a);
}
double safety_factor_static(double yield_stress, double sigma_vm) {
    return yield_stress/sigma_vm;
}
double safety_factor_mohr_coulomb(double c,double sigma_n,double phi,double tau) {
    return (c+sigma_n*tan(phi))/fmax(tau,1e-15);
}
double basquin_sn(double Nf, double sigma_fp, double b) { return sigma_fp*pow(2.0*Nf,b); }
double fatigue_life_basquin(double sigma_a, double sigma_fp, double b) {
    return 0.5*pow(sigma_a/sigma_fp,1.0/b);
}
double strain_life_cm(double Nf,double E,double sigma_fp,double b,double eps_fp,double c) {
    return sigma_fp/E*pow(2.0*Nf,b) + eps_fp*pow(2.0*Nf,c);
}
double goodman_equiv(double sigma_a, double sigma_m, double sigma_uts) {
    return sigma_a/(1.0-sigma_m/sigma_uts);
}
double gerber_equiv(double sigma_a, double sigma_m, double sigma_uts) {
    return sigma_a/(1.0-(sigma_m/sigma_uts)*(sigma_m/sigma_uts));
}
double soderberg_equiv(double sigma_a, double sigma_m, double sigma_y) {
    return sigma_a/(1.0-sigma_m/sigma_y);
}
double morrow_equiv(double sigma_a,double sigma_m,double sigma_fp,double b) {
    return (sigma_fp-sigma_m)*pow(sigma_a/sigma_fp,1.0/b);
}
double miner_damage(const double *cycles, const double *lives, int n) {
    double D=0; for(int i=0;i<n;i++) D+=cycles[i]/lives[i]; return D;
}
int rainflow_count(const double *signal, int n, double *amps, double *means, int *ncycles, int maxc) {
    int cnt=0;
    for(int i=1;i<n && cnt<maxc;i++) {
        amps[cnt]=fabs(signal[i]-signal[i-1]);
        means[cnt]=(signal[i]+signal[i-1])/2.0;
        cnt++;
    }
    *ncycles=cnt; return cnt;
}
double stress_intensity_mode1(double sigma, double a, double beta) {
    return beta*sigma*sqrt(M_PI*a);
}
double stress_intensity_mode2(double tau, double a, double beta) {
    return beta*tau*sqrt(M_PI*a);
}
double stress_intensity_mode3(double tau_yz, double a, double beta) {
    return beta*tau_yz*sqrt(M_PI*a);
}
double energy_release_rate(double K, double E, double nu, int plane_stress) {
    return plane_stress ? K*K/E : K*K*(1.0-nu*nu)/E;
}
int fracture_criterion_check(double K, double Kc) { return K >= Kc ? 1 : 0; }
double paris_crack_growth(double dK, double C, double m) { return C*pow(dK,m); }
double plastic_zone_irwin(double K, double sigma_y, int plane_stress) {
    double fac = plane_stress ? 1.0 : 3.0;
    return fac*K*K/(2.0*M_PI*sigma_y*sigma_y);
}
double j_integral_elastic(double K, double E, double nu) { (void)nu; return K*K/E; }
double ctod_estimate(double K, double E, double sigma_y) {
    return 4.0*K*K/(M_PI*E*sigma_y);
}
double griffith_critical_stress(double E, double gamma_s, double a) {
    return sqrt(2.0*E*gamma_s/(M_PI*a));
}
double mixed_mode_fracture_angle(double K1, double K2) {
    double r=K2/fmax(fabs(K1),1e-15);
    return 2.0*atan((1.0-sqrt(1.0+8.0*r*r))/(4.0*r));
}
double weibull_failure_prob(double sigma, double sigma0, double m) {
    return 1.0-exp(-pow(sigma/sigma0,m));
}
double weibull_scale_effect(double sigma0, double V, double V0, double m) {
    return sigma0*pow(V0/V,1.0/m);
}
