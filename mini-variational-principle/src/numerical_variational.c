/**
 * numerical_variational.c -- Numerical Methods for Variational Problems (L5)
 * FEM, Rayleigh-Ritz, spectral Galerkin, VFD, obstacle problem,
 * error estimation, adaptive refinement, Newton-Kantorovich.
 */
#include "numerical_variational.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void fem_stiffness_1d(const double *xn, int n, double (*D)(double), double *K)
{
    if (!xn||!D||!K||n<2) return;
    for (int i=0;i<n*n;i++) K[i]=0.0;
    for (int e=0;e<n-1;e++){
        double h=xn[e+1]-xn[e],xm=xn[e]+0.5*h,Dv=D(xm);
        double Ke[4]={Dv/h,-Dv/h,-Dv/h,Dv/h};
        K[e*n+e]+=Ke[0]; K[e*n+e+1]+=Ke[1];
        K[(e+1)*n+e]+=Ke[2]; K[(e+1)*n+e+1]+=Ke[3];
    }
}

void fem_mass_1d(const double *xn, int n, double (*rho)(double), double *M)
{
    if (!xn||!rho||!M||n<2) return;
    for (int i=0;i<n*n;i++) M[i]=0.0;
    for (int e=0;e<n-1;e++){
        double h=xn[e+1]-xn[e],xm=xn[e]+0.5*h,rv=rho(xm);
        double Me[4]={rv*h/3.0,rv*h/6.0,rv*h/6.0,rv*h/3.0};
        M[e*n+e]+=Me[0]; M[e*n+e+1]+=Me[1];
        M[(e+1)*n+e]+=Me[2]; M[(e+1)*n+e+1]+=Me[3];
    }
}

double integrate_over_element(double (*f)(double), double h, int ng)
{
    if (!f||ng<1) return 0.0;
    if (ng>=2){
        static const double g2[2]={-0.5773502691896257,0.5773502691896257};
        static const double w2[2]={1.0,1.0};
        double r=0.0;
        for (int k=0;k<2;k++){
            double xm=0.5*h*(g2[k]+1.0);
            r+=w2[k]*f(xm);
        }
        return 0.5*h*r;
    }
    return h*f(0.5*h);
}

void rayleigh_ritz_solve(const double *K, const double *fv, int n,
    const int *bc, const double *bcv, int nbc, double *u)
{
    if (!K||!u||n<1) return;
    /* Simple iteration with Dirichlet BC */
    for (int i=0;i<n;i++) u[i]=0.0;
    for (int b=0;b<nbc;b++) u[bc[b]]=bcv[b];
    for (int it=0;it<1000;it++){
        double md=0.0;
        for (int i=0;i<n;i++){
            int is_bc=0;
            for (int b=0;b<nbc;b++) if (bc[b]==i){ is_bc=1; break; }
            if (is_bc) continue;
            double s=0.0,diag=K[i*n+i];
            if (fabs(diag)<1e-15) diag=1e-15;
            for (int j=0;j<n;j++) if (j!=i) s+=K[i*n+j]*u[j];
            double nu=((fv?fv[i]:0.0)-s)/diag;
            double d=fabs(nu-u[i]); if (d>md) md=d; u[i]=nu;
        }
        if (md<1e-12) break;
    }
}

int generalized_eigenvalue(const double *K, const double *M, int n,
    int nm, double *ev, double *evc)
{
    if (!K||!M||!ev||!evc||n<2||nm<1) return 0;
    /* Power iteration for largest eigenvalue of M^{-1}*K */
    double *v=malloc((size_t)n*sizeof(double));
    if (!v) return 0;
    for (int i=0;i<n;i++) v[i]=1.0/sqrt((double)n);
    for (int mode=0;mode<nm;mode++){
        double lambda=0.0;
        for (int it=0;it<100;it++){
            double *w=malloc((size_t)n*sizeof(double));
            if (!w) break;
            /* w = K*v, then v = M^{-1}*w (simplified: diagonal M) */
            for (int i=0;i<n;i++){
                double s=0.0;
                for (int j=0;j<n;j++) s+=K[i*n+j]*v[j];
                w[i]=s/(M[i*n+i]>1e-15?M[i*n+i]:1e-15);
            }
            double nrm=0.0;
            for (int i=0;i<n;i++) nrm+=w[i]*w[i];
            nrm=sqrt(nrm); if (nrm<1e-15) nrm=1e-15;
            double l2=0.0;
            for (int i=0;i<n;i++){ v[i]=w[i]/nrm; l2+=v[i]*w[i]; }
            double diff=fabs(l2-lambda);
            lambda=l2;
            free(w);
            if (diff<1e-10) break;
        }
        ev[mode]=lambda;
        for (int i=0;i<n;i++) evc[mode*n+i]=v[i];
    }
    free(v);
    return nm;
}

void spectral_galerkin_fourier(double (*op)(int m, double x, double L),
    int nm, double L, double *Km)
{
    if (!op||!Km||nm<1) return;
    for (int i=0;i<nm;i++) for (int j=0;j<nm;j++){
        double s=0.0; int nq=200; double dx=L/(double)(nq-1);
        for (int k=0;k<nq;k++){
            double x=(double)k*dx, w=(k==0||k==nq-1)?0.5:1.0;
            s+=w*sin(M_PI*(i+1)*x/L)*op(j,x,L)*dx;
        }
        Km[i*nm+j]=s;
    }
}

void variational_fd_laplacian(int n, double dx, double *K)
{
    if (!K||n<3) return;
    for (int i=0;i<n*n;i++) K[i]=0.0;
    double id2=1.0/(dx*dx);
    K[0]=id2; K[1]=-id2;
    for (int i=1;i<n-1;i++){
        K[i*n+i-1]=-id2; K[i*n+i]=2.0*id2; K[i*n+i+1]=-id2;
    }
    K[(n-1)*n+n-2]=-id2; K[(n-1)*n+n-1]=id2;
}

void obstacle_problem_solve(const double *K, const double *fv,
    const double *psi, int n, int ni, double *u)
{
    if (!K||!psi||!u||n<2) return;
    for (int i=0;i<n;i++) u[i]=0.0;
    for (int it=0;it<ni;it++){
        for (int i=0;i<n;i++){
            double s=0.0;
            for (int j=0;j<n;j++) if (j!=i) s+=K[i*n+j]*u[j];
            double ui=((fv?fv[i]:0.0)-s)/K[i*n+i];
            if (ui<psi[i]) ui=psi[i];
            u[i]=ui;
        }
    }
}

double zienkiewicz_zhu_error_estimate(const double *fe, const double *sr, int n)
{
    if (!fe||!sr||n<1) return 0.0;
    double eta=0.0,U=0.0;
    for (int i=0;i<n;i++){ double d=fe[i]-sr[i]; eta+=d*d; U+=sr[i]*sr[i]; }
    return sqrt(eta)/sqrt(U>1e-15?U:1e-15);
}

double h_refinement_criterion(double e, double tol, double h, int p)
{
    (void)p; return (e>tol)?h*0.5:h;
}

double vms_stabilization_parameter(double h, double kappa)
{ return h*h/(4.0*kappa>1e-15?4.0*kappa:1e-15); }

double vms_advection_diffusion(double Pe, double h)
{ return h/(2.0*(fabs(Pe)>1e-15?fabs(Pe):1e-15)); }

int morse_index_estimate(const double *H, int n, double tol)
{
    if (!H||n<1) return 0;
    /* Count negative diagonal entries as rough estimate */
    int count=0;
    for (int i=0;i<n;i++) if (H[i*n+i]<-tol) count++;
    return count;
}

void newton_kantorovich_step(
    void (*F)(const double *u, double *r),
    void (*J)(const double *u, double *jac),
    double *u, int n, double tol, int mi)
{
    if (!F||!J||!u||n<1) return;
    double *r=malloc((size_t)n*sizeof(double));
    double *jac=malloc((size_t)(n*n)*sizeof(double));
    double *du=malloc((size_t)n*sizeof(double));
    if (!r||!jac||!du){ free(r); free(jac); free(du); return; }
    for (int it=0;it<mi;it++){
        F(u,r); J(u,jac);
        double nr=0.0;
        for (int i=0;i<n;i++){ nr+=r[i]*r[i]; du[i]=r[i]/jac[i*n+i]; }
        if (sqrt(nr)<tol) break;
        for (int i=0;i<n;i++) u[i]-=du[i];
    }
    free(r); free(jac); free(du);
}

double adaptive_variational_stepsize(double eest, double tol, double h)
{
    double sf=0.9,rn=tol/(eest>1e-15?eest:1e-15);
    double hn=h*sf*pow(rn,1.0/3.0);
    if (hn<0.1*h) hn=0.1*h;
    if (hn>2.0*h) hn=2.0*h;
    return hn;
}
