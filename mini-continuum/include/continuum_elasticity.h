#ifndef CONTINUUM_ELASTICITY_H
#define CONTINUUM_ELASTICITY_H
#include "continuum_types.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Hooke's law isotropic: sigma = lambda*tr(eps)*I + 2*mu*eps
 * Reference: Landau Vol.7 §4, Timoshenko & Goodier Ch.1
 * Knowledge: L4 Fundamental Law — generalized Hooke's law */
StressTensor hookes_law_isotropic(StrainTensor eps, double lambda, double mu);
StrainTensor hookes_law_inverse(StressTensor sigma, double lambda, double mu);
void lame_constants(double E, double nu, double *lambda, double *mu);
void engineering_from_lame(double lambda, double mu, double *E, double *nu);
void elastic_constants_table(double E, double nu,
    double *lam, double *mu, double *G, double *K);
void from_bulk_shear(double K, double G, double *E, double *nu);
/* Plane stress stiffness: Q = E/(1-nu^2)*[[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]]
 * Plane strain stiffness: D uses lambda+2mu */
void plane_stress_stiffness(double E, double nu, double Q[3][3]);
void plane_strain_stiffness(double E, double nu, double D[3][3]);
void plane_stress_to_strain(double E, double nu, double *Ep, double *nup);
/* Orthotropic compliance matrix [S] 6x6 from 9 engineering constants
 * Reference: Ting "Anisotropic Elasticity", Jones "Mechanics of Composites"
 * Symmetry: S_ij = S_ji, 9 independent constants */
void orthotropic_compliance(double E1,double E2,double E3,
    double nu12,double nu13,double nu23,double G12,double G13,double G23,
    double S[6][6]);
void transversely_isotropic_stiffness(double E1,double E3,
    double nu12,double nu31,double G31, double C[6][6]);
/* Strain energy density: W = 0.5*lambda*(tr eps)^2 + mu * sum(eps_ij^2)
 * Complementary energy: W* = (1/2E)*sum(sigma_ij^2) - (nu/2E)*(tr sigma)^2 */
double strain_energy_density(StrainTensor eps, double lambda, double mu);
double complementary_energy_density(StressTensor sigma, double E, double nu);
/* Thermoelastic stress: sigma_th = -E*alpha*dT/(1-2*nu) * I (volumetric)
 * Full: sigma = lambda*tr(eps)*I + 2*mu*eps - (3*lambda+2*mu)*alpha*dT*I */
StressTensor thermoelastic_stress(StrainTensor eps, double lambda, double mu,
    double alpha, double dT);
/* Generalized Hooke's law 6x6: sigma_i = C_ij * eps_j (Voigt)
 * For isotropic: 2 independent constants C11, C12, C44=(C11-C12)/2 */
void stiffness_isotropic_6x6(double lambda, double mu, double C[6][6]);
void compliance_isotropic_6x6(double E, double nu, double S[6][6]);
/* Anisotropic: rotate stiffness by angle theta about z-axis
 * C'_ijkl = R_ip R_jq R_kr R_ls C_pqrs using Bond transformation */
void rotate_stiffness_6x6(const double C[6][6], double theta, double Cp[6][6]);
/* Compute stiffness from compliance: C = S^{-1} (6x6 matrix inverse) */
int stiffness_from_compliance(const double S[6][6], double C[6][6]);
/* Saint-Venant compatibility conditions (2D plane strain):
 * d^2 eps_xx/dy^2 + d^2 eps_yy/dx^2 = 2 d^2 eps_xy/dxdy
 * Checked via finite difference at point (i,j) on grid */
double compatibility_residual_2d(const double *eps_xx, const double *eps_yy,
    const double *eps_xy, int nx, int ny, int i, int j, double dx, double dy);
/* Airy stress function phi(x,y): sigma_xx=d^2phi/dy^2, sigma_yy=d^2phi/dx^2,
 * sigma_xy=-d^2phi/dxdy. Automatically satisfies equilibrium. */
void airy_to_stresses(double d2phi_dx2, double d2phi_dy2, double d2phi_dxdy,
    double *sxx, double *syy, double *sxy);
/* Biharmonic equation for Airy function: nabla^4 phi = 0 (no body forces)
 * In polar coords: (d^2/dr^2+1/r*d/dr+1/r^2*d^2/dtheta^2)^2 phi = 0 */
double biharmonic_polar_residual(double (*phi)(double,double),
    double r, double theta, double dr, double dtheta);
#ifdef __cplusplus
}
#endif
#endif
