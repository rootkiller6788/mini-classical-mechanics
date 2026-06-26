#ifndef CONTINUUM_BEAMS_H
#define CONTINUUM_BEAMS_H
#include "continuum_types.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Euler-Bernoulli beam equation: d^2/dx^2 (EI d^2w/dx^2) = q(x)
 * Assumptions: plane sections remain plane, no shear deformation,
 * small deflections. Reference: Timoshenko & Gere Ch.1-4 */
int euler_bernoulli_fdm(double EI, double L, int nx,
    double *x, double *w, double (*q)(double));
double beam_max_deflection_uniform(double EI, double L, double q);
double beam_max_deflection_point(double EI, double L, double P);
double cantilever_deflection_tip(double EI, double L, double P);
double cantilever_uniform_deflection(double EI, double L, double q, double x);
double beam_bending_moment(double EI, double curvature);
double beam_shear_stress_rect(double V, double b, double h);
/* Timoshenko beam: accounts for shear deformation
 * w_total = w_bending + w_shear
 * w_bend = PL^3/(48EI), w_shear = PL/(4 kappa G A) */
double timoshenko_deflection_point(double EI,double L,double P,double GAKs);
void timoshenko_deflection_components(double EI,double L,double P,double GAKs,
    double *w_bend,double *w_shear,double *w_total);
double shear_deformation_ratio(double EI, double GAKs, double L);
/* Euler buckling: P_cr = n^2 * pi^2 * EI / (K*L)^2
 * K: effective length factor. n: buckling mode number.
 * Reference: Timoshenko & Gere "Theory of Elastic Stability" */
double euler_buckling_load(double EI, double L, int n, double K);
typedef enum { BC_PINNED_PINNED, BC_FIXED_FREE, BC_FIXED_FIXED,
    BC_FIXED_PINNED, BC_CANTILEVER } BucklingBC;
double buckling_K_factor(BucklingBC bc);
double critical_buckling_stress(double E, double slenderness_ratio);
double radius_of_gyration(double I, double A);
double slenderness_ratio(double K, double L, double I, double A);
double buckling_mode_shape(int n, double x, double L, double amplitude);
/* Beam vibration: omega_n = (n*pi/L)^2 * sqrt(EI/(rho*A))
 * Mode shape: phi_n(x) = sin(n*pi*x/L) for simply supported */
int beam_natural_frequencies(double EI,double rhoA,double L,int nmax,double *omegas);
double beam_mode_shape(int n, double x, double L);
double cantilever_fundamental_freq(double EI, double rhoA, double L);
/* Lateral-torsional buckling: M_cr = C_b * pi/L * sqrt(E*Iy*G*J)
 * C_b: moment gradient factor (=1.0 for uniform moment) */
double lateral_torsional_buckling(double E,double Iy,double G,double J,
    double L,double Cb);
/* Secant formula for eccentrically loaded column:
 * sigma_max = P/A * (1 + e*c/r^2 * sec(L/(2r)*sqrt(P/(E*A)))) */
double secant_formula(double P,double A,double e,double c,double I,
    double L,double E);
/* Beam on elastic foundation (Winkler): EI d^4w/dx^4 + k*w = q(x)
 * Characteristic length: l_c = (4*EI/k)^(1/4) */
double winkler_char_length(double EI, double k);
double beam_on_winkler_deflection(double EI, double k, double q, double x);
/* Continuous beam: 3-moment equation (Clapeyron)
 * M_i*L_i/I_i + 2*M_{i+1}*(L_i/I_i+L_{i+1}/I_{i+1}) + M_{i+2}*L_{i+1}/I_{i+1}
 * = -6*(A_i*a_i/(I_i*L_i) + A_{i+1}*b_{i+1}/(I_{i+1}*L_{i+1})) */
int three_moment_equation(int nspans, const double *L, const double *I,
    const double *load_terms, double *M);
/* Curved beam: circumferential stress (Winkler-Bach)
 * sigma = M*y/(A*e*(r_n - y)), r_n = A / integral(dA/r) */
double curved_beam_stress(double M, double A, double e, double rn, double y);
/* Shear center: point where transverse load must be applied
 * to produce bending without twist. For channel section. */
double shear_center_channel(double b, double h, double tf, double tw);
#ifdef __cplusplus
}
#endif
#endif
