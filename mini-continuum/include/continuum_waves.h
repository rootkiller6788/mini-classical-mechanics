#ifndef CONTINUUM_WAVES_H
#define CONTINUUM_WAVES_H
#include "continuum_types.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Elastic wave speeds in isotropic unbounded medium
 * P-wave (dilatational): c_p = sqrt((lambda+2*mu)/rho)
 * S-wave (shear):        c_s = sqrt(mu/rho)
 * Ratio: c_p/c_s = sqrt(2(1-nu)/(1-2nu)) >= sqrt(2)
 * Reference: Achenbach "Wave Propagation in Elastic Solids" */
double p_wave_speed(double lambda, double mu, double rho);
double s_wave_speed(double mu, double rho);
double p_wave_speed_eng(double E, double nu, double rho);
double s_wave_speed_eng(double E, double nu, double rho);
/* Rayleigh surface wave speed (approximation by Bergmann)
 * Exact: satisfies (2-k^2)^2 - 4*sqrt(1-k^2)*sqrt(1-(c_s/c_p)^2*k^2)=0
 * where k=c_R/c_s. Approximation: c_R/c_s = (0.862+1.14*nu)/(1+nu) */
double rayleigh_wave_speed_approx(double lambda, double mu, double rho);
double rayleigh_wave_speed_exact(double lambda, double mu, double rho);
/* Love wave: SH surface wave in layered half-space
 * Exists when c_s1 < c < c_s2. Dispersion: tan(kH*sqrt(c^2/c_s1^2-1))
 * = (mu2/mu1)*sqrt(1-c^2/c_s2^2)/sqrt(c^2/c_s1^2-1) */
void love_wave_speed_range(double mu1,double rho1,double mu2,double rho2,
    double *c_min, double *c_max);
double love_wave_dispersion_residual(double k,double H,double mu1,double mu2,
    double rho1,double rho2,double c);
/* Stoneley wave: interface wave between two elastic half-spaces.
 * Exists for limited range of material contrasts. */
int stoneley_wave_exists(double lam1,double mu1,double rho1,
    double lam2,double mu2,double rho2);
double stoneley_wave_speed(double lam1,double mu1,double rho1,
    double lam2,double mu2,double rho2);
/* Lamb waves in plates: symmetric (S0,S1,...) and antisymmetric (A0,A1,...)
 * Rayleigh-Lamb frequency equation:
 *   tan(beta*h/2)/tan(alpha*h/2) = -(4*k^2*alpha*beta)/(k^2-beta^2)^2 (sym)
 * where alpha^2=omega^2/c_p^2-k^2, beta^2=omega^2/c_s^2-k^2 */
double lamb_wave_disp_symmetric(double omega,double k,double h,
    double cp,double cs);
double lamb_wave_disp_antisymmetric(double omega,double k,double h,
    double cp,double cs);
/* Acoustic impedance: Z = rho * c */
double acoustic_impedance(double rho, double c);
/* Reflection coefficient at normal incidence: R = (Z2-Z1)/(Z2+Z1) */
double reflection_coefficient(double Z1, double Z2);
/* Transmission coefficient at normal incidence: T = 2*Z2/(Z2+Z1) */
double transmission_coefficient(double Z1, double Z2);
/* Snell's law: sin(theta2) = (c2/c1)*sin(theta1).
 * Returns theta2 in radians, or NAN if total internal reflection. */
double snell_angle(double theta1_rad, double c1, double c2);
double critical_angle(double c_slow, double c_fast);
/* Free surface P-SV reflection coefficients (Achenbach §5.4)
 * Incident P-wave at angle theta_p, reflected P and SV.
 * Returns (Rpp, Rps) — amplitude ratios. */
void free_surface_reflection(double cp,double cs,double theta_p,
    double *Rpp, double *Rps);
/* Zoeppritz equations for P-SV at welded interface between two elastic media.
 * Incident P-wave from medium 1. Returns (Rpp,Rps,Tpp,Tps) */
void zoeppritz_pp(double vp1,double vs1,double rho1,
    double vp2,double vs2,double rho2,double theta,
    double *Rpp,double *Rps,double *Tpp,double *Tps);
/* 1D wave equation FDTD solver: d^2u/dt^2 = c^2 d^2u/dx^2
 * Domain [0,L], time [0,T]. Returns number of time steps.
 * u must be pre-allocated [nx][nt+1]. */
int wave_1d_fdtd(double c, double L, double T, int nx, double *x,
    double *t, double **u, double (*u0)(double), double (*v0)(double));
/* d'Alembert solution: u(x,t) = f(x-ct) + g(x+ct),
 * where f=forward wave, g=backward wave, c=wave speed */
double dalembert_1d(double (*f)(double), double (*g)(double),
    double x, double t, double c);
/* Dispersion relation for uniform medium: omega = c * k
 * Phase velocity: v_phase = omega/k = c
 * Group velocity: v_group = d omega/dk = c */
void dispersion_uniform(double k, double c, double *omega, double *v_phase);
/* Group velocity via central finite difference */
double group_velocity(double (*omega)(double), double k, double dk);
/* Timoshenko beam dispersion: accounts for rotary inertia and shear
 * omega^4 * rho*I*rho*A/(kappa*G) - omega^2 * (rho*A + rho*I*k^2*(1+E/(kappa*G)))
 * + E*I*k^4 = 0 */
double timoshenko_dispersion(double k,double E,double G,double rho,
    double nu,double h,double kappa);
/* Attenuation coefficient for viscoelastic waves:
 * alpha(omega) = omega/(2*c*Q), where Q is quality factor */
double attenuation_coefficient(double omega, double c, double Q);
/* Geometric spreading: amplitude ~ 1/sqrt(r) in 2D, ~1/r in 3D */
double geometric_spreading_2d(double r, double r0, double A0);
double geometric_spreading_3d(double r, double r0, double A0);
/* Wave energy density: E = 0.5*rho*v^2 + 0.5*stress*strain */
double wave_energy_density(double rho, double v, double sigma, double eps);
/* Poynting vector for elastic waves: P_i = -sigma_ij * v_j */
void elastic_poynting_vector(StressTensor sigma, Vector3 velocity, Vector3 *P);
#ifdef __cplusplus
}
#endif
#endif
