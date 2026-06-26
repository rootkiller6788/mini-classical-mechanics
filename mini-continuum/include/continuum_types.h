#ifndef CONTINUUM_TYPES_H
#define CONTINUUM_TYPES_H
#include <stddef.h>
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct { double x, y, z; } Vector3;
Vector3 vec3_add(Vector3 a, Vector3 b);
Vector3 vec3_sub(Vector3 a, Vector3 b);
Vector3 vec3_scale(Vector3 v, double s);
double vec3_dot(Vector3 a, Vector3 b);
Vector3 vec3_cross(Vector3 a, Vector3 b);
double vec3_norm(Vector3 v);
Vector3 vec3_normalize(Vector3 v);
typedef struct { double m[3][3]; } Matrix3x3;
Matrix3x3 mat3_identity(void);
Matrix3x3 mat3_zero(void);
Matrix3x3 mat3_add(Matrix3x3 a, Matrix3x3 b);
Matrix3x3 mat3_sub(Matrix3x3 a, Matrix3x3 b);
Matrix3x3 mat3_scale(Matrix3x3 a, double s);
Matrix3x3 mat3_multiply(Matrix3x3 a, Matrix3x3 b);
Matrix3x3 mat3_transpose(Matrix3x3 a);
double mat3_trace(Matrix3x3 a);
double mat3_determinant(Matrix3x3 a);
Matrix3x3 mat3_inverse(Matrix3x3 a);
void mat3_eigenvalues_3x3(Matrix3x3 a, double e[3]);
Matrix3x3 mat3_from_voigt(double v[6]);
typedef struct { double xx, yy, zz, xy, xz, yz; } SymTensor2;
SymTensor2 symtensor_zero(void);
SymTensor2 symtensor_create(double xx, double yy, double zz, double xy, double xz, double yz);
SymTensor2 symtensor_from_matrix(Matrix3x3 m);
Matrix3x3 symtensor_to_matrix(SymTensor2 t);
SymTensor2 symtensor_add(SymTensor2 a, SymTensor2 b);
SymTensor2 symtensor_sub(SymTensor2 a, SymTensor2 b);
SymTensor2 symtensor_scale(SymTensor2 t, double s);
double symtensor_ddot(SymTensor2 a, SymTensor2 b);
double symtensor_norm(SymTensor2 t);
double symtensor_equiv_norm(SymTensor2 t);
typedef SymTensor2 StressTensor;
typedef struct { double I1, J2, J3; } StressInvariants;
StressTensor stress_create(double xx,double yy,double zz,double xy,double xz,double yz);
StressInvariants stress_invariants(StressTensor sigma);
double stress_hydrostatic(StressTensor sigma);
StressTensor stress_deviatoric(StressTensor sigma);
double stress_von_mises(StressTensor sigma);
double stress_tresca(StressTensor sigma);
void stress_principal(StressTensor sigma, double p[3]);
void stress_principal_dirs(StressTensor sigma, Vector3 dirs[3]);
Vector3 stress_traction(StressTensor sigma, Vector3 n);
void stress_octahedral(StressTensor sigma, double *so, double *to);
double stress_lode(StressTensor sigma);
typedef SymTensor2 StrainTensor;
StrainTensor strain_from_disp_grad(Matrix3x3 grad_u);
double strain_dilatation(StrainTensor eps);
StrainTensor strain_deviatoric(StrainTensor eps);
void strain_principal(StrainTensor eps, double p[3]);
double strain_equivalent(StrainTensor eps);
StrainTensor strain_from_eng(double ex,double ey,double ez,double gxy,double gxz,double gyz);
typedef struct { double F11,F12,F13,F21,F22,F23,F31,F32,F33; } DeformationGradient;
DeformationGradient defgrad_identity(void);
Matrix3x3 defgrad_to_matrix(DeformationGradient F);
double defgrad_det(DeformationGradient F);
Matrix3x3 defgrad_right_cg(DeformationGradient F);
Matrix3x3 defgrad_left_cg(DeformationGradient F);
StrainTensor defgrad_green_lagrange(DeformationGradient F);
void defgrad_polar(DeformationGradient F, Matrix3x3 *R, Matrix3x3 *U);
double defgrad_vol_ratio(DeformationGradient F);
typedef struct { double E,nu,rho,G,K,lambda,mu; } ElasticMaterial;
ElasticMaterial elastic_create(double E, double nu, double rho);
ElasticMaterial elastic_from_bulk_shear(double K, double G, double rho);
ElasticMaterial elastic_from_lame(double lam, double mu, double rho);
typedef enum { SYM_ISO, SYM_TRANS_ISO, SYM_ORTHO, SYM_MONO, SYM_TRIC } MaterialSymmetry;
int symmetry_nconst(MaterialSymmetry sym);
typedef enum { PLANE_STRESS, PLANE_STRAIN, AXISYMM } PlaneProblemType;
#ifdef __cplusplus
}
#endif
#endif
