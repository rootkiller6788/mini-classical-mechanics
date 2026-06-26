#include "continuum_types.h"
#include <math.h>
#include <string.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Vector3 operations */
Vector3 vec3_add(Vector3 a, Vector3 b) { Vector3 r={a.x+b.x,a.y+b.y,a.z+b.z}; return r; }
Vector3 vec3_sub(Vector3 a, Vector3 b) { Vector3 r={a.x-b.x,a.y-b.y,a.z-b.z}; return r; }
Vector3 vec3_scale(Vector3 v, double s) { Vector3 r={v.x*s,v.y*s,v.z*s}; return r; }
double vec3_dot(Vector3 a, Vector3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
Vector3 vec3_cross(Vector3 a, Vector3 b) { Vector3 r={a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; return r; }
double vec3_norm(Vector3 v) { return sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }
Vector3 vec3_normalize(Vector3 v) { double n=vec3_norm(v); if(n<1e-30)return v; Vector3 r={v.x/n,v.y/n,v.z/n}; return r; }

/* Matrix3x3 operations */
Matrix3x3 mat3_identity(void) {
    Matrix3x3 r = {{{1,0,0},{0,1,0},{0,0,1}}}; return r;
}
Matrix3x3 mat3_zero(void) { Matrix3x3 r; memset(&r,0,sizeof(r)); return r; }
Matrix3x3 mat3_add(Matrix3x3 a, Matrix3x3 b) {
    Matrix3x3 r; int i,j; for(i=0;i<3;i++)for(j=0;j<3;j++)r.m[i][j]=a.m[i][j]+b.m[i][j]; return r;
}
Matrix3x3 mat3_sub(Matrix3x3 a, Matrix3x3 b) {
    Matrix3x3 r; int i,j; for(i=0;i<3;i++)for(j=0;j<3;j++)r.m[i][j]=a.m[i][j]-b.m[i][j]; return r;
}
Matrix3x3 mat3_scale(Matrix3x3 a, double s) {
    Matrix3x3 r; int i,j; for(i=0;i<3;i++)for(j=0;j<3;j++)r.m[i][j]=a.m[i][j]*s; return r;
}
Matrix3x3 mat3_multiply(Matrix3x3 a, Matrix3x3 b) {
    Matrix3x3 r; int i,j,k; for(i=0;i<3;i++)for(j=0;j<3;j++){r.m[i][j]=0;for(k=0;k<3;k++)r.m[i][j]+=a.m[i][k]*b.m[k][j];} return r;
}
Matrix3x3 mat3_transpose(Matrix3x3 a) {
    Matrix3x3 r; int i,j; for(i=0;i<3;i++)for(j=0;j<3;j++)r.m[i][j]=a.m[j][i]; return r;
}
double mat3_trace(Matrix3x3 a) { return a.m[0][0]+a.m[1][1]+a.m[2][2]; }
double mat3_determinant(Matrix3x3 a) {
    return a.m[0][0]*(a.m[1][1]*a.m[2][2]-a.m[1][2]*a.m[2][1])
          -a.m[0][1]*(a.m[1][0]*a.m[2][2]-a.m[1][2]*a.m[2][0])
          +a.m[0][2]*(a.m[1][0]*a.m[2][1]-a.m[1][1]*a.m[2][0]);
}
Matrix3x3 mat3_inverse(Matrix3x3 a) {
    double det=mat3_determinant(a), invdet=1.0/det; Matrix3x3 r;
    r.m[0][0]=(a.m[1][1]*a.m[2][2]-a.m[1][2]*a.m[2][1])*invdet;
    r.m[0][1]=(a.m[0][2]*a.m[2][1]-a.m[0][1]*a.m[2][2])*invdet;
    r.m[0][2]=(a.m[0][1]*a.m[1][2]-a.m[0][2]*a.m[1][1])*invdet;
    r.m[1][0]=(a.m[1][2]*a.m[2][0]-a.m[1][0]*a.m[2][2])*invdet;
    r.m[1][1]=(a.m[0][0]*a.m[2][2]-a.m[0][2]*a.m[2][0])*invdet;
    r.m[1][2]=(a.m[1][0]*a.m[0][2]-a.m[0][0]*a.m[1][2])*invdet;
    r.m[2][0]=(a.m[1][0]*a.m[2][1]-a.m[2][0]*a.m[1][1])*invdet;
    r.m[2][1]=(a.m[2][0]*a.m[0][1]-a.m[0][0]*a.m[2][1])*invdet;
    r.m[2][2]=(a.m[0][0]*a.m[1][1]-a.m[0][1]*a.m[1][0])*invdet;
    return r;
}
void mat3_eigenvalues_3x3(Matrix3x3 a, double e[3]) {
    double I1=mat3_trace(a); Matrix3x3 a2=mat3_multiply(a,a);
    double I2=(I1*I1-mat3_trace(a2))/2.0, I3=mat3_determinant(a);
    double Q=(I1*I1-3.0*I2)/9.0, R=(2.0*I1*I1*I1-9.0*I1*I2+27.0*I3)/54.0;
    if(Q<0)Q=0;
    double ratio=fmin(fmax(R/sqrt(Q*Q*Q),-1.0),1.0), theta=acos(ratio);
    double sqQ=2.0*sqrt(Q), t0=I1/3.0;
    e[0]=sqQ*cos(theta/3.0)+t0; e[1]=sqQ*cos((theta+2.0*M_PI)/3.0)+t0;
    e[2]=sqQ*cos((theta+4.0*M_PI)/3.0)+t0;
    if(e[0]<e[1]){double t=e[0];e[0]=e[1];e[1]=t;}
    if(e[0]<e[2]){double t=e[0];e[0]=e[2];e[2]=t;}
    if(e[1]<e[2]){double t=e[1];e[1]=e[2];e[2]=t;}
}
Matrix3x3 mat3_from_voigt(double v[6]) {
    Matrix3x3 r={{{v[0],v[3],v[4]},{v[3],v[1],v[5]},{v[4],v[5],v[2]}}}; return r;
}

/* SymTensor2 operations */
SymTensor2 symtensor_zero(void) { SymTensor2 z={0,0,0,0,0,0}; return z; }
SymTensor2 symtensor_create(double xx,double yy,double zz,double xy,double xz,double yz) {
    SymTensor2 t={xx,yy,zz,xy,xz,yz}; return t;
}
SymTensor2 symtensor_from_matrix(Matrix3x3 m) {
    SymTensor2 t={m.m[0][0],m.m[1][1],m.m[2][2],m.m[0][1],m.m[0][2],m.m[1][2]}; return t;
}
Matrix3x3 symtensor_to_matrix(SymTensor2 t) {
    Matrix3x3 r={{{t.xx,t.xy,t.xz},{t.xy,t.yy,t.yz},{t.xz,t.yz,t.zz}}}; return r;
}
SymTensor2 symtensor_add(SymTensor2 a,SymTensor2 b) {
    SymTensor2 t={a.xx+b.xx,a.yy+b.yy,a.zz+b.zz,a.xy+b.xy,a.xz+b.xz,a.yz+b.yz}; return t;
}
SymTensor2 symtensor_sub(SymTensor2 a,SymTensor2 b) {
    SymTensor2 t={a.xx-b.xx,a.yy-b.yy,a.zz-b.zz,a.xy-b.xy,a.xz-b.xz,a.yz-b.yz}; return t;
}
SymTensor2 symtensor_scale(SymTensor2 t,double s) {
    SymTensor2 r={t.xx*s,t.yy*s,t.zz*s,t.xy*s,t.xz*s,t.yz*s}; return r;
}
double symtensor_ddot(SymTensor2 a,SymTensor2 b) {
    return a.xx*b.xx+a.yy*b.yy+a.zz*b.zz+2.0*(a.xy*b.xy+a.xz*b.xz+a.yz*b.yz);
}
double symtensor_norm(SymTensor2 t) { return sqrt(symtensor_ddot(t,t)); }
double symtensor_equiv_norm(SymTensor2 t) {
    double tr=t.xx+t.yy+t.zz,dx=t.xx-tr/3,dy=t.yy-tr/3,dz=t.zz-tr/3;
    return sqrt(2.0/3.0*(dx*dx+dy*dy+dz*dz+2.0*(t.xy*t.xy+t.xz*t.xz+t.yz*t.yz)));
}

/* StressTensor operations */
StressTensor stress_create(double xx,double yy,double zz,double xy,double xz,double yz) {
    return symtensor_create(xx,yy,zz,xy,xz,yz);
}
StressInvariants stress_invariants(StressTensor sigma) {
    StressInvariants inv; inv.I1=sigma.xx+sigma.yy+sigma.zz;
    double p3=inv.I1/3.0,sx=sigma.xx-p3,sy=sigma.yy-p3,sz=sigma.zz-p3;
    inv.J2=0.5*(sx*sx+sy*sy+sz*sz)+sigma.xy*sigma.xy+sigma.xz*sigma.xz+sigma.yz*sigma.yz;
    inv.J3=sx*(sy*sz-sigma.yz*sigma.yz)-sigma.xy*(sigma.xy*sz-sigma.yz*sigma.xz)
           +sigma.xz*(sigma.xy*sigma.yz-sy*sigma.xz);
    return inv;
}
double stress_hydrostatic(StressTensor sigma) { return (sigma.xx+sigma.yy+sigma.zz)/3.0; }
StressTensor stress_deviatoric(StressTensor sigma) {
    double p=stress_hydrostatic(sigma);
    return symtensor_create(sigma.xx-p,sigma.yy-p,sigma.zz-p,sigma.xy,sigma.xz,sigma.yz);
}
double stress_von_mises(StressTensor sigma) {
    StressInvariants inv=stress_invariants(sigma); return sqrt(3.0*inv.J2);
}
double stress_tresca(StressTensor sigma) {
    double p[3]; stress_principal(sigma,p);
    double d12=fabs(p[0]-p[1]),d23=fabs(p[1]-p[2]),d31=fabs(p[2]-p[0]);
    return fmax(d12,fmax(d23,d31));
}
void stress_principal(StressTensor sigma,double p[3]) {
    Matrix3x3 m=symtensor_to_matrix(sigma); mat3_eigenvalues_3x3(m,p);
}
void stress_principal_dirs(StressTensor sigma,Vector3 dirs[3]) {
    double p[3]; stress_principal(sigma,p);
    Matrix3x3 m=symtensor_to_matrix(sigma); int k;
    for(k=0;k<3;k++) {
        double a00=m.m[0][0]-p[k],a01=m.m[0][1],a02=m.m[0][2];
        double a11=m.m[1][1]-p[k],a12=m.m[1][2],a22=m.m[2][2]-p[k];
        double c00=a11*a22-a12*a12,c01=a02*a12-a01*a22,c02=a01*a12-a02*a11;
        double c11=a00*a22-a02*a02,c12=a01*a02-a00*a12,c22=a00*a11-a01*a01;
        double n0=c00*c00+c01*c01+c02*c02,n1=c01*c01+c11*c11+c12*c12,n2=c02*c02+c12*c12+c22*c22;
        Vector3 v;
        if(n0>=n1&&n0>=n2){v.x=c00;v.y=c01;v.z=c02;}
        else if(n1>=n2){v.x=c01;v.y=c11;v.z=c12;}
        else{v.x=c02;v.y=c12;v.z=c22;}
        dirs[k]=vec3_normalize(v);
    }
}
Vector3 stress_traction(StressTensor sigma,Vector3 n) {
    Matrix3x3 m=symtensor_to_matrix(sigma);
    Vector3 r={m.m[0][0]*n.x+m.m[0][1]*n.y+m.m[0][2]*n.z,
               m.m[1][0]*n.x+m.m[1][1]*n.y+m.m[1][2]*n.z,
               m.m[2][0]*n.x+m.m[2][1]*n.y+m.m[2][2]*n.z};
    return r;
}
void stress_octahedral(StressTensor sigma,double *so,double *to) {
    StressInvariants inv=stress_invariants(sigma);
    *so=inv.I1/3.0; *to=sqrt(2.0*inv.J2/3.0);
}
double stress_lode(StressTensor sigma) {
    double p[3]; stress_principal(sigma,p);
    double denom=p[0]-p[2]; if(fabs(denom)<1e-30)return 0.0;
    return (2.0*p[1]-p[0]-p[2])/denom;
}

/* StrainTensor operations */
StrainTensor strain_from_disp_grad(Matrix3x3 grad_u) {
    Matrix3x3 gt=mat3_transpose(grad_u);
    Matrix3x3 eps_mat=mat3_scale(mat3_add(grad_u,gt),0.5);
    return symtensor_from_matrix(eps_mat);
}
double strain_dilatation(StrainTensor eps) { return eps.xx+eps.yy+eps.zz; }
StrainTensor strain_deviatoric(StrainTensor eps) {
    double th=strain_dilatation(eps)/3.0;
    return symtensor_create(eps.xx-th,eps.yy-th,eps.zz-th,eps.xy,eps.xz,eps.yz);
}
void strain_principal(StrainTensor eps,double p[3]) {
    Matrix3x3 m=symtensor_to_matrix(eps); mat3_eigenvalues_3x3(m,p);
}
double strain_equivalent(StrainTensor eps) {
    StrainTensor ed=strain_deviatoric(eps);
    return sqrt(2.0/3.0*symtensor_ddot(ed,ed));
}
StrainTensor strain_from_eng(double ex,double ey,double ez,double gxy,double gxz,double gyz) {
    return symtensor_create(ex,ey,ez,gxy/2.0,gxz/2.0,gyz/2.0);
}

/* DeformationGradient */
DeformationGradient defgrad_identity(void) { DeformationGradient F={1,0,0,0,1,0,0,0,1}; return F; }
Matrix3x3 defgrad_to_matrix(DeformationGradient F) {
    Matrix3x3 r={{{F.F11,F.F12,F.F13},{F.F21,F.F22,F.F23},{F.F31,F.F32,F.F33}}}; return r;
}
double defgrad_det(DeformationGradient F) { return mat3_determinant(defgrad_to_matrix(F)); }
Matrix3x3 defgrad_right_cg(DeformationGradient F) {
    Matrix3x3 m=defgrad_to_matrix(F); return mat3_multiply(mat3_transpose(m),m);
}
Matrix3x3 defgrad_left_cg(DeformationGradient F) {
    Matrix3x3 m=defgrad_to_matrix(F); return mat3_multiply(m,mat3_transpose(m));
}
StrainTensor defgrad_green_lagrange(DeformationGradient F) {
    Matrix3x3 C=defgrad_right_cg(F),I3=mat3_identity();
    Matrix3x3 E=mat3_scale(mat3_sub(C,I3),0.5);
    return symtensor_from_matrix(E);
}
void defgrad_polar(DeformationGradient F,Matrix3x3 *R,Matrix3x3 *U) {
    Matrix3x3 C=defgrad_right_cg(F),Uk=mat3_identity(); int iter;
    for(iter=0;iter<15;iter++) {
        Matrix3x3 Uk_inv=mat3_inverse(Uk),Uk_invT=mat3_transpose(Uk_inv);
        Matrix3x3 CUki=mat3_multiply(C,Uk_invT);
        Uk=mat3_scale(mat3_add(Uk,CUki),0.5);
    }
    *U=Uk; *R=mat3_multiply(defgrad_to_matrix(F),mat3_inverse(Uk));
}
double defgrad_vol_ratio(DeformationGradient F) { return defgrad_det(F); }

/* ElasticMaterial */
ElasticMaterial elastic_create(double E,double nu,double rho) {
    ElasticMaterial m; m.E=E; m.nu=nu; m.rho=rho;
    m.G=E/(2.0*(1.0+nu)); m.K=E/(3.0*(1.0-2.0*nu));
    m.lambda=E*nu/((1.0+nu)*(1.0-2.0*nu)); m.mu=m.G; return m;
}
ElasticMaterial elastic_from_bulk_shear(double K,double G,double rho) {
    double E=9.0*K*G/(3.0*K+G), nu=(3.0*K-2.0*G)/(2.0*(3.0*K+G));
    return elastic_create(E,nu,rho);
}
ElasticMaterial elastic_from_lame(double lam,double mu,double rho) {
    double E=mu*(3.0*lam+2.0*mu)/(lam+mu), nu=lam/(2.0*(lam+mu));
    ElasticMaterial m; m.E=E; m.nu=nu; m.rho=rho;
    m.G=mu; m.K=lam+2.0*mu/3.0; m.lambda=lam; m.mu=mu; return m;
}
int symmetry_nconst(MaterialSymmetry sym) {
    switch(sym){case SYM_ISO:return 2;case SYM_TRANS_ISO:return 5;
    case SYM_ORTHO:return 9;case SYM_MONO:return 13;case SYM_TRIC:return 21;default:return 0;}
}
