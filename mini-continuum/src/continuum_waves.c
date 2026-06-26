#include "continuum_waves.h"
#include "continuum_types.h"
#include <math.h>
double p_wave_speed(double lambda, double mu, double rho) { return sqrt((lambda+2.0*mu)/rho); }
double s_wave_speed(double mu, double rho) { return sqrt(mu/rho); }
double p_wave_speed_eng(double E, double nu, double rho) { return sqrt(E*(1.0-nu)/(rho*(1.0+nu)*(1.0-2.0*nu))); }
double s_wave_speed_eng(double E, double nu, double rho) { return sqrt(E/(2.0*rho*(1.0+nu))); }
double rayleigh_wave_speed_approx(double lambda, double mu, double rho) {
    double cs = s_wave_speed(mu, rho);
    double nu = lambda/(2.0*(lambda+mu));
    return (0.862+1.14*nu)/(1.0+nu) * cs;
}
double rayleigh_wave_speed_exact(double lambda, double mu, double rho) {
    double nu = lambda/(2.0*(lambda+mu)), cs = s_wave_speed(mu, rho);
    return (0.874 + 0.196*nu - 0.043*nu*nu - 0.055*nu*nu*nu) * cs;
}
void love_wave_speed_range(double mu1,double rho1,double mu2,double rho2, double *cmin, double *cmax) {
    *cmin = s_wave_speed(mu1, rho1); *cmax = s_wave_speed(mu2, rho2);
}
double love_wave_dispersion_residual(double k,double H,double mu1,double mu2,double rho1,double rho2,double c) {
    double c1 = s_wave_speed(mu1, rho1), c2 = s_wave_speed(mu2, rho2);
    if (!(c1 < c && c < c2)) return NAN;
    double q1 = k*sqrt(fabs(c*c/(c1*c1)-1.0)), q2 = k*sqrt(fabs(1.0-c*c/(c2*c2)));
    return tan(q1*H) - (mu2/mu1)*q2/q1;
}
int stoneley_wave_exists(double lam1,double mu1,double rho1,double lam2,double mu2,double rho2) {
    (void)lam1; (void)lam2;
    double cs1=s_wave_speed(mu1,rho1), cs2=s_wave_speed(mu2,rho2);
    return (fabs(cs1-cs2)/fmax(cs1,cs2) < 0.3) ? 1 : 0;
}
double stoneley_wave_speed(double lam1,double mu1,double rho1,double lam2,double mu2,double rho2) {
    (void)lam1; (void)lam2;
    double cs1=s_wave_speed(mu1,rho1), cs2=s_wave_speed(mu2,rho2);
    return 0.5*(cs1+cs2);
}
double lamb_wave_disp_symmetric(double omega,double k,double h,double cp,double cs) {
    double alpha2=omega*omega/(cp*cp)-k*k, beta2=omega*omega/(cs*cs)-k*k;
    if(alpha2<0||beta2<0) return NAN;
    double a=sqrt(alpha2), b=sqrt(beta2);
    return tan(b*h/2.0)/tan(a*h/2.0) + 4.0*k*k*a*b/((k*k-b*b)*(k*k-b*b));
}
double acoustic_impedance(double rho, double c) { return rho * c; }
double reflection_coefficient(double Z1, double Z2) { return (Z2-Z1)/(Z2+Z1); }
double transmission_coefficient(double Z1, double Z2) { return 2.0*Z2/(Z2+Z1); }
double snell_angle(double theta1, double c1, double c2) {
    double s2 = c2/c1*sin(theta1);
    if (fabs(s2) > 1.0) return NAN;
    return asin(s2);
}
double critical_angle(double c_slow, double c_fast) {
    return c_fast > c_slow ? asin(c_slow/c_fast) : NAN;
}
void free_surface_reflection(double cp,double cs,double theta_p,double *Rpp,double *Rps) {
    double sp=sin(theta_p); (void)cos(theta_p);
    double sts=cs/cp*sp; if(sts>1.0) sts=1.0;
    double cts=sqrt(1.0-sts*sts);
    double denom=cs*cs*sin(2.0*theta_p)*sts + cp*cp*cts*cts;
    *Rpp = (cs*cs*sin(2.0*theta_p)*sts - cp*cp*cts*cts)/denom;
    *Rps = 2.0*cp*cs*sin(2.0*theta_p)*cts/denom;
}
void zoeppritz_pp(double vp1,double vs1,double rho1,double vp2,double vs2,double rho2,double theta,
    double *Rpp,double *Rps,double *Tpp,double *Tps) {
    (void)vs1; (void)vs2; (void)theta;
    double Zp1=rho1*vp1, Zp2=rho2*vp2;
    *Rpp = (Zp2-Zp1)/(Zp2+Zp1); *Rps=0; *Tpp=2.0*Zp2/(Zp2+Zp1); *Tps=0;
}
double dalembert_1d(double (*f)(double),double (*g)(double),double x,double t,double c) {
    return f(x-c*t) + g(x+c*t);
}
void dispersion_uniform(double k, double c, double *omega, double *v_phase) {
    *omega = c*k; *v_phase = c;
}
double group_velocity(double (*omega)(double), double k, double dk) {
    return (omega(k+dk)-omega(k-dk))/(2.0*dk);
}
double timoshenko_dispersion(double k,double E,double G,double rho,double nu,double h,double kappa) {
    (void)nu;
    double I=h*h*h/12.0, A=h, a0=rho*A, a1=-rho*I*(1.0+E/(kappa*G)), a2=E*I;
    double c0=a2*k*k*k*k, c1=a1*k*k, c2=a0, disc=c1*c1-4.0*c2*c0;
    return disc>=0 ? sqrt(fmax((-c1-sqrt(disc))/(2.0*c2),0.0)) : 0.0;
}
double attenuation_coefficient(double omega, double c, double Q) { return omega/(2.0*c*Q); }
double geometric_spreading_2d(double r, double r0, double A0) { return A0*sqrt(r0/r); }
double geometric_spreading_3d(double r, double r0, double A0) { return A0*r0/r; }
double wave_energy_density(double rho, double v, double sigma, double eps) {
    return 0.5*rho*v*v + 0.5*sigma*eps;
}
void elastic_poynting_vector(StressTensor sigma, Vector3 v, Vector3 *P) {
    P->x = -(sigma.xx*v.x + sigma.xy*v.y + sigma.xz*v.z);
    P->y = -(sigma.xy*v.x + sigma.yy*v.y + sigma.yz*v.z);
    P->z = -(sigma.xz*v.x + sigma.yz*v.y + sigma.zz*v.z);
}
