/**
 * hamilton_principle.c -- Hamilton's Principle (L4-L6)
 */
#include "hamilton_principle.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

double phase_space_action(
    double (*H)(const double *q, const double *p, int n),
    const double *t, const double *q, const double *qdot,
    const double *p, int nd, int ns)
{
    if (!H||!t||!q||!qdot||!p||nd<1||ns<2) return 0.0;
    double S=0.0, dt=t[1]-t[0];
    for (int i=0;i<ns;i++){
        double pq=0.0;
        for (int d=0;d<nd;d++) pq+=p[i*nd+d]*qdot[i*nd+d];
        double w=(i==0||i==ns-1)?0.5:1.0;
        S+=w*(pq-H(q+i*nd,p+i*nd,nd))*dt;
    }
    return S;
}

void phase_space_el_residual(
    double (*H)(const double *q, const double *p, int n),
    void (*gHq)(const double *q, const double *p, int n, double *o),
    void (*gHp)(const double *q, const double *p, int n, double *o),
    const double *q, const double *p,
    const double *qdot, const double *pdot,
    int n, double *res)
{
    (void)H; /* Hamiltonian used via its gradient functions instead */
    if (!gHq||!gHp||!res) return;
    double *gq=malloc((size_t)n*sizeof(double));
    double *gp=malloc((size_t)n*sizeof(double));
    if (!gq||!gp){ free(gq); free(gp); return; }
    gHq(q,p,n,gq); gHp(q,p,n,gp);
    for (int i=0;i<n;i++){ res[i]=qdot[i]-gp[i]; res[n+i]=pdot[i]+gq[i]; }
    free(gq); free(gp);
}

void discrete_hamilton_principle(
    double (*Ld)(const double *qk, const double *qk1, double h),
    const double *qs, int ns, int nd, double h, double *res)
{
    if (!Ld||!qs||!res||ns<3) return;
    double *qk=malloc((size_t)nd*sizeof(double));
    double *qkm=malloc((size_t)nd*sizeof(double));
    if (!qk||!qkm){ free(qk); free(qkm); return; }
    for (int k=1;k<ns-1;k++){
        double eps=1e-6;
        for (int d=0;d<nd;d++){
            memcpy(qk,qs+(k+1)*nd,nd*sizeof(double));
            memcpy(qkm,qs+(k-1)*nd,nd*sizeof(double));
            double orig=qk[d];
            qk[d]+=eps;
            double Lp=Ld(qs+k*nd,qk,h)+Ld(qkm,qs+k*nd,h);
            qk[d]-=2.0*eps;
            double Lm=Ld(qs+k*nd,qk,h)+Ld(qkm,qs+k*nd,h);
            qk[d]=orig;
            res[(k-1)*nd+d]=(Lp-Lm)/(2.0*eps);
        }
    }
    free(qk); free(qkm);
}

double symplectic_2form(const double *dxi, const double *deta, int n)
{
    if (!dxi||!deta||n<1) return 0.0;
    double r=0.0; int nh=n/2;
    for (int i=0;i<nh;i++)
        r+=dxi[i]*deta[nh+i]-dxi[nh+i]*deta[i];
    return r;
}

double momentum_map(const double *G, const double *q, const double *p, int n)
{
    if (!G||!q||!p||n<1) return 0.0;
    double J=0.0;
    for (int i=0;i<n;i++){
        double Giq=0.0;
        for (int j=0;j<n;j++) Giq+=G[i*n+j]*q[j];
        J+=p[i]*Giq;
    }
    return J;
}

bool symplectic_reduction(double J, double mu, double tol)
{ return fabs(J-mu)<tol; }

void hamiltonian_vector_field(
    void (*gH)(const double *q, const double *p, int n, double *o),
    const double *q, const double *p, int n, double *out)
{
    if (!gH||!out) return;
    double *grad=malloc((size_t)(2*n)*sizeof(double));
    if (!grad) return;
    gH(q,p,n,grad);
    int nh=n/2;
    for (int i=0;i<nh;i++){ out[i]=grad[nh+i]; out[nh+i]=-grad[i]; }
    free(grad);
}

double poisson_bracket(
    double (*F)(const double *q, const double *p, int n),
    double (*G)(const double *q, const double *p, int n),
    const double *q, const double *p, int n, double eps)
{
    if (!F||!G) return 0.0;
    double pb=0.0;
    double *qp=malloc((size_t)n*sizeof(double));
    double *pp=malloc((size_t)n*sizeof(double));
    if (!qp||!pp){ free(qp); free(pp); return 0.0; }
    for (int i=0;i<n;i++){
        memcpy(qp,q,n*sizeof(double)); memcpy(pp,p,n*sizeof(double));
        qp[i]+=eps; double Fqp=F(qp,pp,n); qp[i]-=2.0*eps; double Fqm=F(qp,pp,n);
        qp[i]=q[i];
        pp[i]+=eps; double Fpp=F(qp,pp,n); pp[i]-=2.0*eps; double Fpm=F(qp,pp,n);
        pp[i]=p[i];
        memcpy(qp,q,n*sizeof(double)); memcpy(pp,p,n*sizeof(double));
        qp[i]+=eps; double Gqp=G(qp,pp,n); qp[i]-=2.0*eps; double Gqm=G(qp,pp,n);
        qp[i]=q[i];
        pp[i]+=eps; double Gpp=G(qp,pp,n); pp[i]-=2.0*eps; double Gpm=G(qp,pp,n);
        pp[i]=p[i];
        pb+=((Fqp-Fqm)/(2.0*eps))*((Gpp-Gpm)/(2.0*eps))
           -((Fpp-Fpm)/(2.0*eps))*((Gqp-Gqm)/(2.0*eps));
    }
    free(qp); free(pp); return pb;
}

double casimir_check(
    double (*C)(const double *q, const double *p, int n),
    double (*H)(const double *q, const double *p, int n),
    const double *q, const double *p, int n, double eps)
{ return fabs(poisson_bracket(C,H,q,p,n,eps)); }

void multisymplectic_residual(
    const double *K, const double *M,
    const double *z, const double *zt, const double *zx,
    int n, double *res)
{
    if (!K||!M||!z||!zt||!zx||!res) return;
    for (int i=0;i<n;i++){
        double s=0.0;
        for (int j=0;j<n;j++) s+=K[i*n+j]*zt[j]+M[i*n+j]*zx[j];
        res[i]=s;
    }
}

void least_action_path(
    double (*L)(double t, const double *q, const double *qdot, int n),
    void (*gqL)(double t, const double *q, const double *qdot, int n, double *o),
    void (*gqdL)(double t, const double *q, const double *qdot, int n, double *o),
    const double *q0, const double *qT,
    double t0, double tF, int nd, int ns,
    double *qp, int mi, double lr)
{
    if (!L||!gqL||!gqdL||!q0||!qT||!qp||nd<1||ns<3) return;
    double dt=(tF-t0)/(double)(ns-1);
    /* Initialize with linear path */
    for (int i=0;i<ns;i++){
        double t=t0+(double)i*dt, s=(t-t0)/(tF-t0);
        for (int d=0;d<nd;d++) qp[i*nd+d]=q0[d]+s*(qT[d]-q0[d]);
    }
    /* Gradient descent */
    double *g=malloc((size_t)(ns*nd)*sizeof(double));
    if (!g) return;
    for (int it=0;it<mi;it++){
        memset(g,0,(size_t)(ns*nd)*sizeof(double));
        for (int i=1;i<ns-1;i++){
            double t=t0+(double)i*dt;
            double *qd=malloc((size_t)nd*sizeof(double));
            if (!qd) continue;
            /* qdot via central difference */
            for (int d=0;d<nd;d++)
                qd[d]=(qp[(i+1)*nd+d]-qp[(i-1)*nd+d])/(2.0*dt);
            double *gq=malloc((size_t)nd*sizeof(double));
            double *gqd=malloc((size_t)nd*sizeof(double));
            if (gq&&gqd){
                gqL(t,qp+i*nd,qd,nd,gq); gqdL(t,qp+i*nd,qd,nd,gqd);
                for (int d=0;d<nd;d++)
                    g[i*nd+d]=gq[d]-(gqd[(i+1)*nd+d]-gqd[(i-1)*nd+d])/(2.0*dt);
            }
            free(qd); free(gq); free(gqd);
        }
        double mg=0.0;
        for (int i=1;i<ns-1;i++)
            for (int d=0;d<nd;d++){
                qp[i*nd+d]-=lr*g[i*nd+d];
                double ag=fabs(g[i*nd+d]); if(ag>mg) mg=ag;
            }
        if (mg<1e-10) break;
    }
    free(g);
}

double noether_charge_energy(
    double (*L)(double t, const double *q, const double *qdot, int n),
    double t, const double *q, const double *qdot, int n, double eps)
{
    if (!L||!q||!qdot||n<1) return 0.0;
    double L0=L(t,q,qdot,n);
    double *qd=malloc((size_t)n*sizeof(double));
    if (!qd) return 0.0;
    double E=0.0;
    for (int i=0;i<n;i++){
        memcpy(qd,qdot,n*sizeof(double));
        qd[i]+=eps; double Lp=L(t,q,qd,n);
        qd[i]-=2.0*eps; double Lm=L(t,q,qd,n);
        double dLdqdot=(Lp-Lm)/(2.0*eps);
        E+=dLdqdot*qdot[i];
    }
    free(qd);
    return E-L0;
}
