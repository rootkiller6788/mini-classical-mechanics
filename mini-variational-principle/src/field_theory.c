/**
 * field_theory.c -- Field Theory Variational Principles (L6-L8)
 * Ginzburg-Landau, phase field (Allen-Cahn, Cahn-Hilliard),
 * elastic energy, electromagnetic action.
 */
#include "field_theory.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double ginzburg_landau_free_energy(const double *psi, const double *dp,
    double a, double b, double c, double x0, double x1, int n)
{
    if (!psi||!dp||n<2) return 0.0;
    double F=0.0,dx=(x1-x0)/(double)(n-1);
    for (int i=0;i<n;i++){
        double p=psi[i],d=dp[i],w=(i==0||i==n-1)?0.5:1.0;
        F+=w*(a*p*p+0.5*b*p*p*p*p+c*d*d)*dx;
    }
    return F;
}

void solve_ginzburg_landau(double a, double b, double c, double L,
    double bcL, double bcR, int nx, int ni, double omega,
    double *x, double *psi)
{
    if (!x||!psi||nx<3) return;
    double dx=L/(double)(nx-1);
    for (int i=0;i<nx;i++){ x[i]=(double)i*dx; psi[i]=0.5*(bcL+bcR); }
    psi[0]=bcL; psi[nx-1]=bcR;
    for (int it=0;it<ni;it++){
        double *po=malloc((size_t)nx*sizeof(double));
        if (!po) return;
        memcpy(po,psi,(size_t)nx*sizeof(double));
        double md=0.0;
        for (int i=1;i<nx-1;i++){
            double lap=(psi[i+1]-2.0*psi[i]+psi[i-1])/(dx*dx);
            double rhs=(2.0*psi[i]+a*psi[i]+b*psi[i]*psi[i]*psi[i])/c;
            double denom=2.0/(dx*dx)+a/c; if(denom<1e-10) denom=1e-10;
            double np=(1.0-omega)*psi[i]+omega*(lap+rhs)/denom;
            double d=fabs(np-psi[i]); if (d>md) md=d; psi[i]=np;
        }
        free(po);
        if (md<1e-10) break;
    }
}

double gl_coherence_length(double a, double c)
{ return (fabs(a)<1e-15)?INFINITY:sqrt(c/fabs(a)); }

double gl_critical_field(double a, double b)
{ return (b<1e-15)?INFINITY:sqrt(4.0*M_PI*a*a/b); }

void allen_cahn_rhs(const double *phi, double eps, double M,
    double dx, int n, double *dp)
{
    if (!phi||!dp||n<3) return;
    for (int i=1;i<n-1;i++){
        double lap=(phi[i+1]-2.0*phi[i]+phi[i-1])/(dx*dx);
        double fp=phi[i]*(phi[i]*phi[i]-1.0);
        dp[i]=M*(eps*eps*lap-fp);
    }
    dp[0]=0.0; dp[n-1]=0.0;
}

void cahn_hilliard_rhs(const double *phi, double eps, double M,
    double dx, int n, double *dp)
{
    if (!phi||!dp||n<5) return;
    double *mu=malloc((size_t)n*sizeof(double));
    if (!mu) return;
    for (int i=1;i<n-1;i++){
        double lap=-(phi[i+1]-2.0*phi[i]+phi[i-1])/(dx*dx);
        double fp=phi[i]*(phi[i]*phi[i]-1.0);
        mu[i]=eps*eps*lap+fp;
    }
    mu[0]=mu[1]; mu[n-1]=mu[n-2];
    for (int i=2;i<n-2;i++)
        dp[i]=M*(mu[i+1]-2.0*mu[i]+mu[i-1])/(dx*dx);
    dp[0]=dp[1]=dp[2]; dp[n-1]=dp[n-2]=dp[n-3];
    free(mu);
}

double interface_energy(const double *phi, double eps, double dx, int n)
{
    if (!phi||n<2) return 0.0;
    double E=0.0;
    for (int i=1;i<n;i++){
        double g=(phi[i]-phi[i-1])/dx, pm=0.5*(phi[i]+phi[i-1]);
        double f=(pm*pm-1.0)*(pm*pm-1.0)*0.25;
        E+=(0.5*eps*eps*g*g+f)*dx;
    }
    return E;
}

double helmholtz_free_energy(const double *strain, const double *C,
    double alpha, double dT, double lam, double mu)
{
    if (!strain||!C) return 0.0;
    double W=0.0;
    for (int i=0;i<6;i++) for (int j=0;j<6;j++) W+=0.5*strain[i]*C[i*6+j]*strain[j];
    double trE=strain[0]+strain[1]+strain[2];
    return W-(3.0*lam+2.0*mu)*alpha*dT*trE;
}

double gibbs_free_energy(const double *stress, const double *S)
{
    if (!stress||!S) return 0.0;
    double G=0.0;
    for (int i=0;i<6;i++) for (int j=0;j<6;j++) G-=0.5*stress[i]*S[i*6+j]*stress[j];
    return G;
}

double em_action_density(const double *E, const double *B,
    double rho, const double *J, const double *A, double phi)
{
    double F2=0.0;
    for (int i=0;i<3;i++){ F2+=E[i]*E[i]-B[i]*B[i]; }
    double JdotA=0.0;
    if (J&&A) for (int i=0;i<3;i++) JdotA+=J[i]*A[i];
    return -0.25*F2-rho*phi+JdotA;
}

void stress_energy_tensor(
    double (*Ld)(const double *f, const double *g, int n),
    const double *f, const double *g, int n, double *T)
{
    if (!Ld||!f||!g||!T||n<1) return;
    double L0=Ld(f,g,n);
    for (int mu=0;mu<4;mu++){
        for (int nu=0;nu<4;nu++){
            double sum=0.0;
            for (int i=0;i<n;i++)
                sum+=g[mu*n+i]*g[nu*n+i]; /* simplified */
            T[mu*4+nu]=sum-(mu==nu?L0:0.0);
        }
    }
}

double topological_charge_density(const double *phi, double dx, int n)
{
    if (!phi||n<3) return 0.0;
    double Q=0.0;
    for (int i=1;i<n;i++){
        double dphi=(phi[i]-phi[i-1])/dx;
        Q+=dphi*dx;
    }
    return Q/(2.0*M_PI);
}
