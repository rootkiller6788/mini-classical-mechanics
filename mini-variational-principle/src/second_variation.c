/**
 * second_variation.c -- Second Variation Theory (L4-L6)
 * Legendre, Jacobi, Weierstrass conditions, Hilbert integral,
 * sufficient conditions, Ritz/Galerkin methods, Morse theory.
 */
#include "second_variation.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

bool legendre_condition(double v) { return v >= -1e-12; }
bool strong_legendre(double v) { return v > 1e-12; }

void verify_legendre_along_extremal(
    LagrangianDensity d2Fyp2, const double *x, const double *y,
    const double *yp, int n, bool *ok, double *mn, double *mx, double *pv)
{
    if (!d2Fyp2||!x||!y||!yp||n<2) { if(ok)*ok=false; if(mn)*mn=0; if(mx)*mx=0; return; }
    bool p=true; double a=DBL_MAX,b=-DBL_MAX;
    for (int i=0;i<n;i++) {
        double P=d2Fyp2(x[i],y[i],yp[i]); if(pv) pv[i]=P;
        if(P<a)a=P; if(P>b)b=P; if(P<-1e-10) p=false;
    }
    if(ok)*ok=p; if(mn)*mn=a; if(mx)*mx=b;
}

int find_conjugate_points(double a, double b,
    LagrangianDensity P_func, LagrangianDensity Q_func,
    const double *xe, const double *ye, const double *ype,
    int ng, double *cout, int mxp)
{
    if (!P_func||!Q_func||!xe||!ye||!ype||ng<3) return 0;
    double dx=(b-a)/(double)(ng-1);
    double *u=malloc((size_t)ng*sizeof(double));
    double *up=malloc((size_t)ng*sizeof(double));
    if (!u||!up) { free(u); free(up); return 0; }
    u[0]=0.0; up[0]=1.0;
    for (int i=0;i<ng-1;i++) {
        double xm=xe[i]+0.5*dx, P=P_func(xm,0.0,0.0), Q=Q_func(xm,0.0,0.0);
        if (fabs(P)<1e-15) P=1e-15;
        up[i+1]=up[i]+dx*Q*u[i]/P; u[i+1]=u[i]+dx*up[i];
    }
    int np=0;
    for (int i=1;i<ng&&np<mxp;i++) {
        if (u[i-1]*u[i]<0.0&&cout) {
            double t=fabs(u[i-1])/(fabs(u[i-1])+fabs(u[i]));
            cout[np++]=xe[i-1]+t*dx;
        }
    }
    free(u); free(up); return np;
}

bool jacobi_condition(const double *cp, int n, double a, double b) {
    for (int i=0;i<n;i++) if (cp[i]>a&&cp[i]<b) return false;
    return true;
}
double first_conjugate_point(const double *cp, int n, double a) {
    double f=INFINITY;
    for (int i=0;i<n;i++) if (cp[i]>a&&cp[i]<f) f=cp[i];
    return f;
}
void solve_jacobi_equation(double a, double b,
    LagrangianDensity Pf, LagrangianDensity Qf,
    const double *xe, const double *ye, const double *ype,
    int ng, double u0p, double *uo, double *upo)
{
    if (!Pf||!Qf||!xe||!uo||!upo||ng<2) return;
    double dx=(b-a)/(double)(ng-1); uo[0]=0.0; upo[0]=u0p;
    for (int i=0;i<ng-1;i++) {
        double xm=xe[i]+0.5*dx, P=Pf(xm,0.0,0.0), Q=Qf(xm,0.0,0.0);
        if (fabs(P)<1e-15) P=1e-15;
        upo[i+1]=upo[i]+dx*Q*uo[i]/P; uo[i+1]=uo[i]+dx*upo[i];
    }
}

double weierstrass_excess(LagrangianDensity F, LagrangianDensity dFdp,
    double x, double y, double po, double p) {
    if (!F||!dFdp) return 0.0;
    return F(x,y,p)-F(x,y,po)-(p-po)*dFdp(x,y,po);
}
bool weierstrass_condition(LagrangianDensity F, LagrangianDensity dFdp,
    double x, double y, double po, const double *pr, int np) {
    if (!F||!dFdp||!pr||np<1) return false;
    for (int i=0;i<np;i++)
        if (weierstrass_excess(F,dFdp,x,y,po,pr[i])<-1e-10) return false;
    return true;
}
bool verify_weierstrass_along_extremal(LagrangianDensity F,
    LagrangianDensity dFdp, const double *x, const double *y,
    const double *yp, int nn, const double *pr, int np) {
    if (!F||!dFdp||!x||!y||!yp||nn<1) return false;
    for (int i=0;i<nn;i++)
        if (!weierstrass_condition(F,dFdp,x[i],y[i],yp[i],pr,np)) return false;
    return true;
}

double hilbert_invariant_integral(LagrangianDensity F, LagrangianDensity dFdp,
    const double *x, const double *y, const double *yp,
    double (*pf)(double,double), int n) {
    if (!F||!dFdp||!x||!y||!yp||!pf||n<2) return 0.0;
    double I=0.0,dx=x[1]-x[0];
    for (int i=0;i<n;i++) {
        double p=pf(x[i],y[i]), w=(i==0||i==n-1)?0.5:1.0;
        I+=w*(F(x[i],y[i],p)+(yp[i]-p)*dFdp(x[i],y[i],p))*dx;
    }
    return I;
}
double hilbert_integral_test(LagrangianDensity F, LagrangianDensity dFdp,
    double (*pf)(double,double), const double *x,
    const double *y1, const double *yp1,
    const double *y2, const double *yp2, int n) {
    double I1=hilbert_invariant_integral(F,dFdp,x,y1,yp1,pf,n);
    double I2=hilbert_invariant_integral(F,dFdp,x,y2,yp2,pf,n);
    return fabs(I1-I2);
}

SufficientConditions check_sufficient_conditions(
    LagrangianDensity F, LagrangianDensity dFdp, LagrangianDensity d2Fyp2,
    LagrangianDensity Pf, LagrangianDensity Qf,
    const double *x, const double *y, const double *yp,
    int n, const double *pr, int np)
{
    SufficientConditions sc; memset(&sc,0,sizeof(sc));
    bool lp=false; double mn,mx;
    double *pv=malloc((size_t)n*sizeof(double));
    if (pv) {
        verify_legendre_along_extremal(d2Fyp2,x,y,yp,n,&lp,&mn,&mx,pv);
        sc.legendre_ok=lp&&(mn>1e-12); free(pv);
    }
    double cp[100];
    int nc=find_conjugate_points(x[0],x[n-1],Pf,Qf,x,y,yp,n,cp,100);
    sc.jacobi_ok=jacobi_condition(cp,nc,x[0],x[n-1]);
    sc.weierstrass_ok=verify_weierstrass_along_extremal(F,dFdp,x,y,yp,n,pr,np);
    sc.sufficient=sc.legendre_ok&&sc.jacobi_ok&&sc.weierstrass_ok;
    return sc;
}

int ritz_method(LagrangianDensity d2Fyp2, LagrangianDensity d2Fy2,
    double a, double b, int nb, double (**ba)(double),
    double (**bd)(double), int nq, double *co)
{
    if (!d2Fyp2||!d2Fy2||!ba||!bd||!co||nb<1||nq<2) return -1;
    double *A=calloc((size_t)(nb*nb),sizeof(double));
    double *fv=calloc((size_t)nb,sizeof(double));
    if (!A||!fv) { free(A); free(fv); return -1; }
    double dx=(b-a)/(double)(nq-1);
    for (int k=0;k<nq;k++) {
        double x=a+(double)k*dx, w=(k==0||k==nq-1)?0.5:1.0;
        double P=d2Fyp2(x,0.0,0.0), Q=d2Fy2(x,0.0,0.0);
        for (int i=0;i<nb;i++)
            for (int j=0;j<nb;j++)
                A[i*nb+j]+=w*(P*bd[i](x)*bd[j](x)+Q*ba[i](x)*ba[j](x))*dx;
    }
    for (int iter=0;iter<1000;iter++) {
        double md=0.0;
        for (int i=0;i<nb;i++) {
            double s=0.0;
            for (int j=0;j<nb;j++) if (j!=i) s+=A[i*nb+j]*co[j];
            double nc=-s/A[i*nb+i];
            double d=fabs(nc-co[i]); if (d>md) md=d; co[i]=nc;
        }
        if (md<1e-12) break;
    }
    free(A); free(fv); return 0;
}

int galerkin_method(double (*el_op)(double(*u)(double),double x),
    double a, double b, int nb, double (**ba)(double),
    int nq, double *co)
{
    if (!el_op||!ba||!co||nb<1||nq<2) return -1;
    double *A=calloc((size_t)(nb*nb),sizeof(double));
    if (!A) return -1;
    double dx=(b-a)/(double)(nq-1);
    for (int k=0;k<nq;k++) {
        double x=a+(double)k*dx, w=(k==0||k==nq-1)?0.5:1.0;
        for (int i=0;i<nb;i++) {
            double ev=el_op(ba[i],x);
            for (int j=0;j<nb;j++) A[i*nb+j]+=w*ev*ba[j](x)*dx;
        }
    }
    for (int iter=0;iter<500;iter++) {
        double md=0.0;
        for (int i=0;i<nb;i++) {
            double s=0.0;
            for (int j=0;j<nb;j++) if (j!=i) s+=A[i*nb+j]*co[j];
            double nc=-s/A[i*nb+i]; double d=fabs(nc-co[i]);
            if (d>md) md=d; co[i]=nc;
        }
        if (md<1e-10) break;
    }
    free(A); return 0;
}

int morse_index(const double *cp, int n, double a, double b) {
    int c=0; for (int i=0;i<n;i++) if (cp[i]>a&&cp[i]<b) c++; return c;
}
int count_jacobi_zeros(const double *u, int n) {
    if (!u||n<2) return 0;
    int c=0; for (int i=1;i<n;i++) if (u[i-1]*u[i]<0.0) c++; return c;
}
