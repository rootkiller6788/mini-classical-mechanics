#ifndef CONTINUUM_FAILURE_H
#define CONTINUUM_FAILURE_H
#include "continuum_types.h"
#ifdef __cplusplus
extern "C" {
#endif
/* YIELD CRITERIA
 * von Mises: sigma_vm = sqrt(3*J2) <= sigma_y
 *   Physical interpretation: distortional strain energy
 * Tresca: max(|sigma1-sigma2|,|sigma2-sigma3|,|sigma3-sigma1|) <= sigma_y
 *   Physical interpretation: maximum shear stress
 * Reference: Hill "Mathematical Theory of Plasticity" */
double von_mises_yield_func(StressTensor sigma, double sigma_y);
double tresca_yield_func(StressTensor sigma, double sigma_y);
/* Rankine (max principal stress): max(sigma1,sigma2,sigma3) <= sigma_ut */
double rankine_yield_func(StressTensor sigma, double sigma_ut);
/* Saint-Venant (max principal strain) */
double saint_venant_yield_func(StressTensor sigma, double E, double nu, double eps_ut);
/* Mohr-Coulomb: tau = c + sigma_n*tan(phi)
 * Principal form: sigma1/sigma_ut - sigma3/sigma_uc <= 1 */
double mohr_coulomb_shear_strength(double c, double sigma_n, double phi_rad);
double mohr_coulomb_envelope(StressTensor sigma, double sigma_ut, double sigma_uc);
/* Drucker-Prager: sqrt(J2) + alpha*I1 <= k
 * alpha = 2*sin(phi)/(sqrt(3)*(3-sin(phi))), k = 6*c*cos(phi)/(sqrt(3)*(3-sin(phi)))
 * Matches Mohr-Coulomb outer apices for 3D stress states */
double drucker_prager_yield_func(StressTensor sigma, double alpha, double k);
void drucker_prager_params_mc(double c, double phi_rad, double *alpha, double *k);
/* Hoek-Brown for rock: sigma1 = sigma3 + sigma_ci*(m_b*sigma3/sigma_ci + s)^a */
double hoek_brown(double sigma3, double sigma_ci, double mb, double s, double a);
/* SAFETY FACTORS */
double safety_factor_static(double yield_stress, double sigma_vm);
double safety_factor_mohr_coulomb(double c,double sigma_n,double phi_rad,double tau);
/* FATIGUE
 * Basquin (high-cycle): sigma_a = sigma_f' * (2*Nf)^b
 * Coffin-Manson (low-cycle): eps_a = (sigma_f'/E)*(2*Nf)^b + eps_f'*(2*Nf)^c
 * Reference: Dowling "Mechanical Behavior of Materials", Suresh "Fatigue" */
double basquin_sn(double Nf, double sigma_fp, double b);
double fatigue_life_basquin(double sigma_a, double sigma_fp, double b);
double strain_life_cm(double Nf, double E, double sigma_fp, double b,
    double eps_fp, double c);
/* Mean stress corrections:
 * Goodman: sigma_a/sigma_ar + sigma_m/sigma_uts = 1
 * Gerber:  sigma_a/sigma_ar + (sigma_m/sigma_uts)^2 = 1
 * Soderberg: sigma_a/sigma_ar + sigma_m/sigma_y = 1
 * Morrow: sigma_ar = (sigma_f' - sigma_m)*(2*Nf)^b */
double goodman_equiv(double sigma_a, double sigma_m, double sigma_uts);
double gerber_equiv(double sigma_a, double sigma_m, double sigma_uts);
double soderberg_equiv(double sigma_a, double sigma_m, double sigma_y);
double morrow_equiv(double sigma_a,double sigma_m,double sigma_fp,double b);
/* Palmgren-Miner linear damage accumulation: D = sum(n_i/N_fi). D>=1 = failure */
double miner_damage(const double *cycles, const double *lives, int n);
/* Rainflow counting (simplified 4-point algorithm)
 * Input: stress/strain time history, length n
 * Output: cycles array (amplitudes), means array, ncycles set to count */
int rainflow_count(const double *signal, int n,
    double *amplitudes, double *means, int *ncycles, int max_cycles);
/* FRACTURE MECHANICS
 * Stress intensity factor: K_I = beta * sigma * sqrt(pi*a)
 * Energy release rate (plane stress): G = K_I^2/E
 * Paris law: da/dN = C*(Delta K)^m
 * Reference: Anderson "Fracture Mechanics" */
double stress_intensity_mode1(double sigma, double a, double beta);
double stress_intensity_mode2(double tau, double a, double beta);
double stress_intensity_mode3(double tau_yz, double a, double beta);
double energy_release_rate(double K, double E, double nu, int plane_stress);
int fracture_criterion_check(double K, double Kc);
double paris_crack_growth(double dK, double C, double m);
double plastic_zone_irwin(double K, double sigma_y, int plane_stress);
double j_integral_elastic(double K, double E, double nu);
double ctod_estimate(double K, double E, double sigma_y);
double griffith_critical_stress(double E, double gamma_s, double a);
double mixed_mode_fracture_angle(double K1, double K2);
/* Weibull weakest-link statistics: Pf = 1-exp(-(sigma/sigma0)^m)
 * Scale effect: sigma0_V = sigma0_0 * (V0/V)^(1/m) */
double weibull_failure_prob(double sigma, double sigma0, double m);
double weibull_scale_effect(double sigma0, double V, double V0, double m);
#ifdef __cplusplus
}
#endif
#endif
