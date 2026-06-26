#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "functional.h"
#include "euler_lagrange.h"
#include "constrained.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int passed = 0, failed = 0;
#define T(n, e) do { if(e) passed++; else { failed++; printf("FAIL: %s\n",n); } } while(0)

static double freeF(double x, double y, double yp) { (void)x;(void)y; return 0.5*yp*yp; }
static double harm_dFdy(double x, double y, double yp) { (void)x;(void)yp; return -y; }
static double harm_dFdyp(double x, double y, double yp) { (void)x;(void)y; return yp; }
static double harm_d2Fdyp2(double x, double y, double yp) { (void)x;(void)y;(void)yp; return 1.0; }
static double linf(double x) { return 2.0*x+1.0; }
static double lind(double x) { (void)x; return 2.0; }

int main(void) {
    printf("=== mini-variational-principle Test Suite ===\n\n");
    int n=101;
    double *x=malloc((size_t)n*sizeof(double));
    double *y=malloc((size_t)n*sizeof(double));
    double *yp=malloc((size_t)n*sizeof(double));
    generate_uniform_nodes(0.0,1.0,n,x);
    evaluate_function_samples(linf,lind,x,n,y,yp);

    printf("--- L2: Functional Evaluation ---\n");
    double J1=evaluate_functional(freeF,x,y,yp,n);
    T("J_trap>0",J1>0.0);
    double J2=evaluate_functional_simpson(freeF,x,y,yp,n);
    T("J_simp>0",J2>0.0);
    T("quad_agree",fabs(J1-J2)<0.01);
    T("norm>0",functional_norm(y,yp,0.0,1.0,n)>0.0);
    T("L2>0",functional_L2_norm(y,0.0,1.0,n)>0.0);

    printf("--- L3: Derivatives ---\n");
    double *eta=malloc((size_t)n*sizeof(double));
    double *etap=malloc((size_t)n*sizeof(double));
    double *res=malloc((size_t)n*sizeof(double));
    for (int i=0;i<n;i++){ y[i]=x[i]; yp[i]=1.0; }
    for (int i=0;i<n;i++){ eta[i]=sin(x[i]); etap[i]=cos(x[i]); }
    double dJ=gateaux_derivative(freeF,x,y,yp,eta,etap,n,1e-4);
    T("gateaux_finite",isfinite(dJ));
    functional_derivative(harm_dFdy,harm_dFdyp,x,y,yp,n,res);
    double mr=0.0; for(int i=0;i<n;i++) if(fabs(res[i])>mr) mr=fabs(res[i]);
    T("fd_residual_ok",mr<100.0);
    T("frechet_finite",isfinite(frechet_derivative_test(freeF,x,y,yp,eta,etap,n)));

    printf("--- L4: Beltrami ---\n");
    for (int i=0;i<n;i++){ y[i]=3.0*x[i]; yp[i]=3.0; }
    double *vals=malloc((size_t)n*sizeof(double));
    double be=verify_beltrami(freeF,harm_dFdyp,x,y,yp,n,vals);
    T("beltrami_constant",be<1e-6);
    T("beltrami_val",isfinite(beltrami_identity(freeF,harm_dFdyp,0.0,0.0,3.0)));

    printf("--- L4: Variations ---\n");
    for (int i=0;i<n;i++){ y[i]=x[i]; yp[i]=1.0; eta[i]=x[i]*(1.0-x[i]); etap[i]=1.0-2.0*x[i]; }
    T("1st_var_finite",isfinite(first_variation_formula(harm_dFdy,harm_dFdyp,x,y,yp,eta,etap,n)));
    T("2nd_var_>=0",second_variation_formula(harm_d2Fdyp2,harm_dFdy,x,eta,etap,n)>=-1e-10);

    printf("--- L5: Utilities ---\n");
    double xu[10]; generate_uniform_nodes(0.0,1.0,10,xu);
    T("uni_start",fabs(xu[0])<1e-15); T("uni_end",fabs(xu[9]-1.0)<1e-15);
    double xc[10]; generate_chebyshev_nodes(-1.0,1.0,10,xc);
    T("cheb_start",fabs(xc[0]-1.0)<1e-15); T("cheb_end",fabs(xc[9]+1.0)<1e-15);

    printf("--- L6: Canonical Systems ---\n");
    T("brach_>0",brachistochrone_lagrangian(0.5,1.0,0.5,9.81)>0.0);
    T("catenary_>0",catenary_lagrangian(0.0,1.0,0.5)>0.0);
    T("harm_finite",isfinite(harmonic_oscillator_lagrangian(0.0,0.5,1.0,1.0,1.0)));
    double R; dido_solution(2.0*M_PI,100,NULL,NULL,&R);
    T("dido_R1",fabs(R-1.0)<1e-10);

    printf("--- L7: Constrained ---\n");
    T("iso_lambda",isfinite(isoperimetric_lagrange_multiplier(freeF,freeF,0.5,0.0,1.0,0.0,1.0,0.0)));
    T("nat_bc",isfinite(variable_left_endpoint_condition(freeF,harm_dFdyp,0.0,0.0,1.0)));

    printf("--- L3: Convexity ---\n");
    for (int i=0;i<n;i++){ y[i]=x[i]; yp[i]=1.0; eta[i]=2.0*x[i]; etap[i]=2.0; }
    T("convex",is_convex_functional(freeF,x,y,yp,eta,etap,n,5));
    T("quadratic",is_quadratic_functional(freeF,x,y,yp,n));

    free(x);free(y);free(yp);free(eta);free(etap);free(res);free(vals);
    printf("\n=== Results: %d passed, %d failed ===\n",passed,failed);
    return failed>0?1:0;
}
