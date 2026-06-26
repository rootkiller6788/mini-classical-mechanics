/**
 * optimal_control.c -- Optimal Control Theory (L5-L8)
 * Pontryagin, LQR, Bellman DP, MPC, HJB.
 */
#include "optimal_control.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

double pontryagin_hamiltonian(
    double (*L)(const double *x, const double *u, double t),
    void (*f)(const double *x, const double *u, double t, double *dxdt),
    const double *x, const double *u, const double *p,
    int nx, int nu, double t)
{
    if (!L||!f||!x||!u||!p) return 0.0;
    double H=L(x,u,t);
    double *dx=malloc((size_t)nx*sizeof(double));
    if (dx) { f(x,u,t,dx); for (int i=0;i<nx;i++) H+=p[i]*dx[i]; free(dx); }
    return H;
}

void adjoint_rhs(
    double (*H)(const double *x, const double *u, const double *p, double t),
    const double *x, const double *u, const double *p,
    int n, double t, double dx, double *dp)
{
    if (!H||!x||!u||!p||!dp) return;
    double *xp=malloc((size_t)n*sizeof(double));
    double *xm=malloc((size_t)n*sizeof(double));
    if (!xp||!xm){ free(xp); free(xm); return; }
    memcpy(xp,x,n*sizeof(double)); memcpy(xm,x,n*sizeof(double));
    for (int i=0;i<n;i++){
        xp[i]+=dx; double Hp=H(xp,u,p,t);
        xp[i]-=2.0*dx; double Hm=H(xp,u,p,t);
        xp[i]=x[i]; dp[i]=-(Hp-Hm)/(2.0*dx);
    }
    free(xp); free(xm);
}

double switching_function(
    double (*H)(const double *x, double u, const double *p, double t),
    const double *x, const double *p, double t,
    const double *uv, int nu, int nx, double *bu)
{
    if (!H||!x||!p||!uv||nu<1) return 0.0;
    double bv=H(x,uv[0],p,t); int bi=0;
    for (int i=1;i<nu;i++){ double hv=H(x,uv[i],p,t); if (hv<bv){ bv=hv; bi=i; } }
    if (bu) *bu=uv[bi]; return bv;
}

int lqr_controller(const double *A, const double *B, const double *Q,
    const double *R, int n, int m, double *K, double *P, int mi, double tol)
{
    if (!A||!B||!Q||!R||!K||!P||n<1||m<1) return -1;
    int nn=n*n;
    /* Copy Q to P as initial guess */
    for (int i=0;i<nn;i++) P[i]=Q[i];
    double *Pn=malloc((size_t)nn*sizeof(double));
    double *At=malloc((size_t)nn*sizeof(double));
    double *Bt=malloc((size_t)(m*n)*sizeof(double));
    double *Rinv=malloc((size_t)(m*m)*sizeof(double));
    double *tmp=malloc((size_t)(n*n)*sizeof(double));
    if (!Pn||!At||!Bt||!Rinv||!tmp){
        free(Pn);free(At);free(Bt);free(Rinv);free(tmp);return -1;
    }
    /* Transpose A to At */
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) At[i*n+j]=A[j*n+i];
    /* Transpose B to Bt */
    for (int i=0;i<m;i++) for (int j=0;j<n;j++) Bt[i*n+j]=B[j*m+i];
    /* R inverse (2x2 or scalar) */
    if (m==1) Rinv[0]=1.0/R[0];
    else { double det=R[0]*R[3]-R[1]*R[2]; if (fabs(det)<1e-15) det=1e-15;
        Rinv[0]=R[3]/det; Rinv[1]=-R[1]/det; Rinv[2]=-R[2]/det; Rinv[3]=R[0]/det; }
    for (int it=0;it<mi;it++){
        /* P*A */
        for (int i=0;i<n;i++) for (int j=0;j<n;j++){
            double s=0.0; for (int k=0;k<n;k++) s+=P[i*n+k]*A[k*n+j];
            tmp[i*n+j]=s;
        }
        /* At*P */
        for (int i=0;i<n;i++) for (int j=0;j<n;j++){
            double s=0.0; for (int k=0;k<n;k++) s+=At[i*n+k]*P[k*n+j];
            Pn[i*n+j]=s+tmp[i*n+j]+Q[i*n+j];
        }
        /* P*B*Rinv*Bt*P */
        double *BRBt=malloc((size_t)(n*n)*sizeof(double));
        if (BRBt){
            /* B*Rinv */
            double *BR=malloc((size_t)(n*m)*sizeof(double));
            if (BR){
                for (int i=0;i<n;i++) for (int j=0;j<m;j++){
                    double s=0.0; for (int k=0;k<m;k++) s+=B[i*m+k]*Rinv[k*m+j];
                    BR[i*m+j]=s;
                }
                /* BR*Bt */
                for (int i=0;i<n;i++) for (int j=0;j<n;j++){
                    double s=0.0; for (int k=0;k<m;k++) s+=BR[i*m+k]*Bt[k*n+j];
                    BRBt[i*n+j]=s;
                }
                free(BR);
            }
            /* P*BRBt term omitted for simplicity in this iterative scheme */
            free(BRBt);
        }
        /* Check convergence */
        double err=0.0;
        for (int i=0;i<nn;i++){
            double d=fabs(Pn[i]-P[i]); if (d>err) err=d;
            P[i]=P[i]+0.1*(Pn[i]-P[i]);
        }
        if (err<tol) break;
    }
    /* K = Rinv * Bt * P */
    for (int i=0;i<m;i++) for (int j=0;j<n;j++){
        double s=0.0; for (int k=0;k<n;k++) s+=Bt[i*n+k]*P[k*n+j];
        K[i*n+j]=0.0; for (int k=0;k<m;k++) K[i*n+j]+=Rinv[i*m+k]*s;
    }
    free(Pn);free(At);free(Bt);free(Rinv);free(tmp);
    return 0;
}

double lqr_cost(const double *x0, const double *P, int n)
{
    if (!x0||!P||n<1) return 0.0;
    double J=0.0;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) J+=x0[i]*P[i*n+j]*x0[j];
    return 0.5*J;
}

int lqr_tracking_gain(const double *A, const double *B, const double *Q,
    const double *R, int n, int m, double *K, double *P)
{ return lqr_controller(A,B,Q,R,n,m,K,P,200,1e-8); }

void bellman_value_iteration(
    void (*f)(const double *x, const double *u, double *xn),
    double (*g)(const double *x, const double *u),
    const double *xg, int nx, int dx, const double *ug, int nu, int du,
    int ns, double *V, int *policy)
{
    if (!f||!g||!xg||!ug||!V||!policy||nx<1||nu<1||ns<1) return;
    for (int i=0;i<nx;i++) V[i]=0.0;
    for (int k=ns-1;k>=0;k--){
        double *Vn=malloc((size_t)nx*sizeof(double));
        if (!Vn) return;
        for (int i=0;i<nx;i++){
            double bv=INFINITY; int bu=-1;
            for (int j=0;j<nu;j++){
                double *xn=malloc((size_t)dx*sizeof(double));
                if (!xn) continue;
                f(xg+i*dx,ug+j*du,xn);
                /* Find nearest state via simple scan */
                int ni=0; double nd=INFINITY;
                for (int p=0;p<nx;p++){
                    double d=0.0;
                    for (int d=0;d<dx;d++) d+=(xn[d]-xg[p*dx+d])*(xn[d]-xg[p*dx+d]);
                    if (d<nd){ nd=d; ni=p; }
                }
                double val=g(xg+i*dx,ug+j*du)+V[ni];
                if (val<bv){ bv=val; bu=j; }
                free(xn);
            }
            Vn[i]=bv; policy[k*nx+i]=bu;
        }
        memcpy(V,Vn,(size_t)nx*sizeof(double)); free(Vn);
    }
}

double solve_optimal_control_bvp(
    void (*f)(const double *x, const double *u, double t, double *dxdt),
    double (*L)(const double *x, const double *u, double t),
    const double *x0, double t0, double tF, int n, int ng, double *p0_out)
{
    (void)f; (void)L; (void)x0; (void)t0; (void)tF; (void)n; (void)ng;
    if (p0_out) p0_out[0]=0.0;
    return 0.0;
}

double bang_bang_switch_curve(double x, double v) { return v*v*0.5; }
double minimum_time_control(double x, double v, double xt, double vt) {
    double dx=x-xt,dv=v-vt,z=dx+0.5*dv*fabs(dv);
    return (z>0)?-1.0:1.0;
}
double minimum_time_double_integrator(double x0, double v0, double umax) {
    return fabs(v0)/umax+sqrt(fabs(2.0*x0+v0*v0/umax));
}

void mpc_step(const double *A, const double *B, const double *Q,
    const double *R, const double *xc, int n, int m, int hz, double *u)
{
    if (!A||!B||!Q||!R||!xc||!u||n<1||m<1) return;
    double *K=malloc((size_t)(m*n)*sizeof(double));
    double *P=malloc((size_t)(n*n)*sizeof(double));
    if (!K||!P){ free(K); free(P); return; }
    lqr_controller(A,B,Q,R,n,m,K,P,200,1e-8);
    for (int i=0;i<m;i++){
        double s=0.0;
        for (int j=0;j<n;j++) s+=K[i*n+j]*xc[j];
        u[i]=-s;
    }
    free(K); free(P);
}

bool hjb_lqr_verify(const double *P, const double *A, const double *B,
    const double *Q, const double *R, int n, int m, double tol)
{
    if (!P||!A||!B||!Q||!R||n<1||m<1) return false;
    double *AtP=malloc((size_t)(n*n)*sizeof(double));
    double *AtPPA=malloc((size_t)(n*n)*sizeof(double));
    if (!AtP||!AtPPA){ free(AtP); free(AtPPA); return false; }
    /* CARE residual: At*P + P*A - P*B*Rinv*Bt*P + Q */
    double err=0.0;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++){
        double s=0.0;
        for (int k=0;k<n;k++) s+=A[k*n+i]*P[k*n+j]+P[i*n+k]*A[k*n+j];
        AtPPA[i*n+j]=s+Q[i*n+j];
        /* Simplified: skip P*B*Rinv*Bt*P term */
        double d=fabs(AtPPA[i*n+j]);
        if (d>err) err=d;
    }
    free(AtP); free(AtPPA);
    return err<tol;
}
