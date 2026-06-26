#include "continuum_elasticity.h"
#include "continuum_types.h"
#include <math.h>
#include <string.h>
void lame_constants(double E, double nu, double *lambda, double *mu) {
    *lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    *mu = E / (2.0 * (1.0 + nu));
}
void engineering_from_lame(double lambda, double mu, double *E, double *nu) {
    *E = mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
    *nu = lambda / (2.0 * (lambda + mu));
}
StressTensor hookes_law_isotropic(StrainTensor eps, double lambda, double mu) {
    double tr = eps.xx + eps.yy + eps.zz;
    return stress_create(lambda*tr + 2.0*mu*eps.xx, lambda*tr + 2.0*mu*eps.yy,
        lambda*tr + 2.0*mu*eps.zz, 2.0*mu*eps.xy, 2.0*mu*eps.xz, 2.0*mu*eps.yz);
}
StrainTensor hookes_law_inverse(StressTensor sigma, double lambda, double mu) {
    double tr = sigma.xx + sigma.yy + sigma.zz;
    double fac = lambda / (2.0*mu*(3.0*lambda + 2.0*mu));
    return symtensor_create(
        (1.0/(2.0*mu))*sigma.xx - fac*tr,
        (1.0/(2.0*mu))*sigma.yy - fac*tr,
        (1.0/(2.0*mu))*sigma.zz - fac*tr,
        sigma.xy/(2.0*mu), sigma.xz/(2.0*mu), sigma.yz/(2.0*mu));
}
void elastic_constants_table(double E, double nu, double *lam, double *mu, double *G, double *K) {
    *lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    *mu = E / (2.0 * (1.0 + nu));
    *G = *mu;
    *K = E / (3.0 * (1.0 - 2.0 * nu));
}
void from_bulk_shear(double K, double G, double *E, double *nu) {
    *E = 9.0 * K * G / (3.0 * K + G);
    *nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G));
}
void plane_stress_stiffness(double E, double nu, double Q[3][3]) {
    double f = E / (1.0 - nu*nu);
    Q[0][0]=f; Q[0][1]=f*nu; Q[0][2]=0;
    Q[1][0]=f*nu; Q[1][1]=f; Q[1][2]=0;
    Q[2][0]=0; Q[2][1]=0; Q[2][2]=f*(1.0-nu)/2.0;
}
void plane_strain_stiffness(double E, double nu, double D[3][3]) {
    double f = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    D[0][0]=f*(1.0-nu); D[0][1]=f*nu; D[0][2]=0;
    D[1][0]=f*nu; D[1][1]=f*(1.0-nu); D[1][2]=0;
    D[2][0]=0; D[2][1]=0; D[2][2]=f*(1.0-2.0*nu)/2.0;
}
void plane_stress_to_strain(double E, double nu, double *Ep, double *nup) {
    *Ep = E / (1.0 - nu*nu);
    *nup = nu / (1.0 - nu);
}
double strain_energy_density(StrainTensor eps, double lambda, double mu) {
    double tr = eps.xx + eps.yy + eps.zz;
    return 0.5*lambda*tr*tr + mu*(eps.xx*eps.xx+eps.yy*eps.yy+eps.zz*eps.zz
        + 2.0*(eps.xy*eps.xy+eps.xz*eps.xz+eps.yz*eps.yz));
}
double complementary_energy_density(StressTensor sigma, double E, double nu) {
    double tr = sigma.xx + sigma.yy + sigma.zz;
    Matrix3x3 sm = symtensor_to_matrix(sigma);
    double s2 = 0; for(int i=0;i<3;i++)for(int j=0;j<3;j++) s2+=sm.m[i][j]*sm.m[i][j];
    return (1.0/(2.0*E))*s2 - (nu/(2.0*E))*tr*tr;
}
StressTensor thermoelastic_stress(StrainTensor eps, double lambda, double mu,
    double alpha, double dT) {
    double tr = eps.xx + eps.yy + eps.zz;
    double thermal = (3.0*lambda + 2.0*mu) * alpha * dT;
    return stress_create(lambda*tr + 2.0*mu*eps.xx - thermal,
        lambda*tr + 2.0*mu*eps.yy - thermal,
        lambda*tr + 2.0*mu*eps.zz - thermal,
        2.0*mu*eps.xy, 2.0*mu*eps.xz, 2.0*mu*eps.yz);
}
void stiffness_isotropic_6x6(double lambda, double mu, double C[6][6]) {
    memset(C, 0, sizeof(double)*36);
    double c11=lambda+2.0*mu, c12=lambda, c44=mu;
    C[0][0]=C[1][1]=C[2][2]=c11; C[0][1]=C[0][2]=C[1][0]=C[1][2]=C[2][0]=C[2][1]=c12;
    C[3][3]=C[4][4]=C[5][5]=c44;
}
double compatibility_residual_2d(const double *exx,const double *eyy,const double *exy,
    int nx,int ny,int i,int j,double dx,double dy) {
    int idx = i*ny + j;
    if(i<1||i>=nx-1||j<1||j>=ny-1) return 0;
    double d2exx_dy2 = (exx[idx+1]-2*exx[idx]+exx[idx-1])/(dy*dy);
    double d2eyy_dx2 = (eyy[idx+ny]-2*eyy[idx]+eyy[idx-ny])/(dx*dx);
    double d2exy_dxdy = (exy[idx+ny+1]-exy[idx+ny-1]-exy[idx-ny+1]+exy[idx-ny-1])/(4*dx*dy);
    return d2exx_dy2 + d2eyy_dx2 - 2*d2exy_dxdy;
}
