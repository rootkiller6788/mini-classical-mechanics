/**
 * rigid_body_inertia.c — Inertia Tensor: Computation, Principal Axes, Standard Shapes
 *
 * Goldstein §5.1-5.3: Moment of inertia tensor definitions and computations.
 * All formulas verified against Goldstein Table 5.1.
 *
 * Key implementation decisions:
 *   - Jacobi eigenvalue algorithm for 3×3 symmetric diagonalization
 *   - All standard shape formulas are analytically exact (no numerical integration error)
 *   - Parallel axis theorem: exact algebraic transform
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/rigid_body_types.h"
#include "../include/rigid_body_inertia.h"

/* ============================================================================
 * Inertia Tensor from Discrete Particles
 * ============================================================================ */

void inertia_tensor_from_particles(
    int n_particles, const double masses[], const vec3 positions[],
    InertiaTensor *out)
{
    double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
    for (int k = 0; k < n_particles; k++) {
        double m = masses[k];
        double x = positions[k].x, y = positions[k].y, z = positions[k].z;
        double r2 = x*x + y*y + z*z;
        Ixx += m * (r2 - x*x);
        Iyy += m * (r2 - y*y);
        Izz += m * (r2 - z*z);
        Ixy -= m * x * y;
        Ixz -= m * x * z;
        Iyz -= m * y * z;
    }
    *out = (InertiaTensor){Ixx, Iyy, Izz, Ixy, Ixz, Iyz};
}

void inertia_tensor_numerical(
    double (*density_func)(double x, double y, double z),
    const double xlim[2], const double ylim[2], const double zlim[2],
    int n_grid, InertiaTensor *out)
{
    double dx = (xlim[1] - xlim[0]) / n_grid;
    double dy = (ylim[1] - ylim[0]) / n_grid;
    double dz = (zlim[1] - zlim[0]) / n_grid;
    double dV = dx * dy * dz;

    double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
    double M = 0; /* Total mass */
    double cm_x = 0, cm_y = 0, cm_z = 0; /* Center of mass accumulators */

    /* First pass: total mass and CM */
    for (int i = 0; i < n_grid; i++) {
        double x = xlim[0] + (i + 0.5) * dx;
        for (int j = 0; j < n_grid; j++) {
            double y = ylim[0] + (j + 0.5) * dy;
            for (int k = 0; k < n_grid; k++) {
                double z = zlim[0] + (k + 0.5) * dz;
                double rho = density_func(x, y, z);
                double dm = rho * dV;
                M += dm;
                cm_x += dm * x;
                cm_y += dm * y;
                cm_z += dm * z;
            }
        }
    }

    if (M < 1e-15) {
        *out = INERTIA_ZERO;
        return;
    }
    cm_x /= M; cm_y /= M; cm_z /= M;

    /* Second pass: inertia about CM */
    for (int i = 0; i < n_grid; i++) {
        double x = xlim[0] + (i + 0.5) * dx - cm_x;
        for (int j = 0; j < n_grid; j++) {
            double y = ylim[0] + (j + 0.5) * dy - cm_y;
            for (int k = 0; k < n_grid; k++) {
                double z = zlim[0] + (k + 0.5) * dz - cm_z;
                double rho = density_func(x + cm_x, y + cm_y, z + cm_z);
                double dm = rho * dV;
                double r2 = x*x + y*y + z*z;
                Ixx += dm * (r2 - x*x);
                Iyy += dm * (r2 - y*y);
                Izz += dm * (r2 - z*z);
                Ixy -= dm * x * y;
                Ixz -= dm * x * z;
                Iyz -= dm * y * z;
            }
        }
    }

    *out = (InertiaTensor){Ixx, Iyy, Izz, Ixy, Ixz, Iyz};
}

/* ============================================================================
 * Parallel Axis Theorem
 * ============================================================================ */

void parallel_axis_theorem(const InertiaTensor *I_cm, double M, vec3 d, InertiaTensor *out) {
    double d2 = d.x*d.x + d.y*d.y + d.z*d.z;
    out->Ixx = I_cm->Ixx + M * (d2 - d.x*d.x);
    out->Iyy = I_cm->Iyy + M * (d2 - d.y*d.y);
    out->Izz = I_cm->Izz + M * (d2 - d.z*d.z);
    out->Ixy = I_cm->Ixy - M * d.x * d.y;
    out->Ixz = I_cm->Ixz - M * d.x * d.z;
    out->Iyz = I_cm->Iyz - M * d.y * d.z;
}

void parallel_axis_inverse(const InertiaTensor *I_p, double M, vec3 d, InertiaTensor *out) {
    double d2 = d.x*d.x + d.y*d.y + d.z*d.z;
    out->Ixx = I_p->Ixx - M * (d2 - d.x*d.x);
    out->Iyy = I_p->Iyy - M * (d2 - d.y*d.y);
    out->Izz = I_p->Izz - M * (d2 - d.z*d.z);
    out->Ixy = I_p->Ixy + M * d.x * d.y;
    out->Ixz = I_p->Ixz + M * d.x * d.z;
    out->Iyz = I_p->Iyz + M * d.y * d.z;
}

/* ============================================================================
 * Principal Axes via Jacobi Eigenvalue Algorithm for 3×3 Symmetric Matrix
 *
 * Jacobi rotation method converges for real symmetric matrices.
 * Each rotation zeros the largest off-diagonal element.
 * For 3×3, convergence is very fast (typically < 10 sweeps).
 * ============================================================================ */

static double jacobi_max_offdiag(const double A[9], int *p, int *q) {
    double max_val = 0.0;
    *p = 0; *q = 1;
    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            double val = fabs(A[i * 3 + j]);
            if (val > max_val) {
                max_val = val;
                *p = i;
                *q = j;
            }
        }
    }
    return max_val;
}

void principal_axes_decompose(const InertiaTensor *I, PrincipalAxes *out) {
    /* Build symmetric 3×3 matrix */
    double A[9];
    A[0] = I->Ixx; A[1] = I->Ixy; A[2] = I->Ixz;
    A[3] = I->Ixy; A[4] = I->Iyy; A[5] = I->Iyz;
    A[6] = I->Ixz; A[7] = I->Iyz; A[8] = I->Izz;

    /* Initialize eigenvector matrix to identity */
    double V[9] = {1,0,0, 0,1,0, 0,0,1};

    /* Jacobi iterations (max 50 sweeps for 3×3, typically converges in < 10) */
    for (int sweep = 0; sweep < 50; sweep++) {
        int p, q;
        double max_off = jacobi_max_offdiag(A, &p, &q);
        if (max_off < 1e-14) break;

        double theta;
        if (fabs(A[p*3+p] - A[q*3+q]) < 1e-15) {
            theta = M_PI / 4.0;
        } else {
            theta = 0.5 * atan2(2.0 * A[p*3+q], A[p*3+p] - A[q*3+q]);
        }

        double c = cos(theta), s = sin(theta);

        /* Apply Jacobi rotation to A: A' = J^T · A · J */
        /* Update rows/cols p and q */
        double App = A[p*3+p], Aqq = A[q*3+q], Apq = A[p*3+q];

        A[p*3+p] = c*c * App - 2*s*c * Apq + s*s * Aqq;
        A[q*3+q] = s*s * App + 2*s*c * Apq + c*c * Aqq;
        A[p*3+q] = 0.0;
        A[q*3+p] = 0.0;

        for (int i = 0; i < 3; i++) {
            if (i != p && i != q) {
                double Aip = A[i*3+p], Aiq = A[i*3+q];
                A[i*3+p] = c * Aip - s * Aiq;
                A[p*3+i] = A[i*3+p];
                A[i*3+q] = s * Aip + c * Aiq;
                A[q*3+i] = A[i*3+q];
            }
        }

        /* Update eigenvector matrix V: V' = V · J */
        for (int i = 0; i < 3; i++) {
            double Vip = V[i*3+p], Viq = V[i*3+q];
            V[i*3+p] = c * Vip - s * Viq;
            V[i*3+q] = s * Vip + c * Viq;
        }
    }

    /* Extract eigenvalues (diagonal of A) and sort descending */
    double evals[3] = {A[0], A[4], A[8]};
    int idx[3] = {0, 1, 2};

    /* Bubble sort (3 elements) by eigenvalue descending */
    for (int i = 0; i < 2; i++) {
        for (int j = i + 1; j < 3; j++) {
            if (evals[idx[i]] < evals[idx[j]]) {
                int tmp = idx[i]; idx[i] = idx[j]; idx[j] = tmp;
            }
        }
    }

    out->moments[0] = evals[idx[0]];
    out->moments[1] = evals[idx[1]];
    out->moments[2] = evals[idx[2]];

    /* Place eigenvectors in columns sorted by descending eigenvalue */
    for (int col = 0; col < 3; col++) {
        out->axes.m[0 + col] = V[0 + idx[col]];
        out->axes.m[3 + col] = V[3 + idx[col]];
        out->axes.m[6 + col] = V[6 + idx[col]];
    }
}

void principal_moments_only(const InertiaTensor *I, double moments[3]) {
    PrincipalAxes pa;
    principal_axes_decompose(I, &pa);
    moments[0] = pa.moments[0];
    moments[1] = pa.moments[1];
    moments[2] = pa.moments[2];
}

int is_inertia_diagonal(const InertiaTensor *I, double tol) {
    return (fabs(I->Ixy) < tol) && (fabs(I->Ixz) < tol) && (fabs(I->Iyz) < tol);
}

/* ============================================================================
 * Standard Homogeneous Shapes — Inertia Tensors about CM
 * ============================================================================ */

InertiaTensor sphere_inertia(double M, double R) {
    double I = 0.4 * M * R * R; /* (2/5) M R² */
    return (InertiaTensor){I, I, I, 0.0, 0.0, 0.0};
}

InertiaTensor spherical_shell_inertia(double M, double R) {
    double I = (2.0/3.0) * M * R * R;
    return (InertiaTensor){I, I, I, 0.0, 0.0, 0.0};
}

InertiaTensor cylinder_inertia(double M, double R, double H) {
    double Izz = 0.5 * M * R * R;
    double Iperp = M * (3.0 * R * R + H * H) / 12.0;
    return (InertiaTensor){Iperp, Iperp, Izz, 0.0, 0.0, 0.0};
}

InertiaTensor hollow_cylinder_inertia(double M, double R, double H) {
    double Izz = M * R * R;
    double Iperp = 0.5 * M * R * R + M * H * H / 12.0;
    return (InertiaTensor){Iperp, Iperp, Izz, 0.0, 0.0, 0.0};
}

InertiaTensor cuboid_inertia(double M, double a, double b, double c) {
    return (InertiaTensor){
        M * (b*b + c*c) / 12.0,
        M * (a*a + c*c) / 12.0,
        M * (a*a + b*b) / 12.0,
        0.0, 0.0, 0.0
    };
}

InertiaTensor rod_inertia(double M, double L) {
    double Iperp = M * L * L / 12.0;
    return (InertiaTensor){Iperp, Iperp, 0.0, 0.0, 0.0, 0.0};
}

InertiaTensor disk_inertia(double M, double R) {
    double Izz = 0.5 * M * R * R;
    double Iperp = 0.25 * M * R * R;
    return (InertiaTensor){Iperp, Iperp, Izz, 0.0, 0.0, 0.0};
}

InertiaTensor ellipsoid_inertia(double M, double a, double b, double c) {
    return (InertiaTensor){
        0.2 * M * (b*b + c*c),
        0.2 * M * (a*a + c*c),
        0.2 * M * (a*a + b*b),
        0.0, 0.0, 0.0
    };
}

InertiaTensor solid_cone_inertia(double M, double R, double H) {
    /* CM is at z = 3H/4 from apex. In CM frame: */
    double Izz_cm = 0.3 * M * R * R;
    double Iperp_cm = 0.15 * M * R * R + 0.0375 * M * H * H;
    /* 3/20 = 0.15, 3/80 = 0.0375 */
    return (InertiaTensor){Iperp_cm, Iperp_cm, Izz_cm, 0.0, 0.0, 0.0};
}

InertiaTensor torus_inertia(double M, double R, double r) {
    double Izz = M * (R*R + 0.75 * r*r);
    double Iperp = 0.5 * M * R*R + 0.625 * M * r*r;
    /* 5/8 = 0.625 */
    return (InertiaTensor){Iperp, Iperp, Izz, 0.0, 0.0, 0.0};
}

void inertia_compose(int n_parts, const InertiaTensor parts[], InertiaTensor *out) {
    double Ixx = 0, Iyy = 0, Izz = 0, Ixy = 0, Ixz = 0, Iyz = 0;
    for (int i = 0; i < n_parts; i++) {
        Ixx += parts[i].Ixx; Iyy += parts[i].Iyy; Izz += parts[i].Izz;
        Ixy += parts[i].Ixy; Ixz += parts[i].Ixz; Iyz += parts[i].Iyz;
    }
    *out = (InertiaTensor){Ixx, Iyy, Izz, Ixy, Ixz, Iyz};
}

/* ============================================================================
 * Inertia Ellipsoid
 * ============================================================================ */

void inertia_ellipsoid_semiaxes(const InertiaTensor *I, double *a, double *b, double *c) {
    double moments[3];
    principal_moments_only(I, moments);
    *a = 1.0 / sqrt(moments[0]); /* shortest: largest inertia */
    *b = 1.0 / sqrt(moments[1]);
    *c = 1.0 / sqrt(moments[2]); /* longest: smallest inertia */
}

void inertia_ellipsoid_points(const InertiaTensor *I, int n_theta, int n_phi, vec3 points[]) {
    PrincipalAxes pa;
    principal_axes_decompose(I, &pa);
    double a, b_ax, c_ax;
    inertia_ellipsoid_semiaxes(I, &a, &b_ax, &c_ax);

    int idx = 0;
    for (int it = 0; it < n_theta; it++) {
        double theta = M_PI * it / (n_theta - 1);
        for (int ip = 0; ip < n_phi; ip++) {
            double phi = M_2PI * ip / n_phi;
            double xl = a * sin(theta) * cos(phi);
            double yl = b_ax * sin(theta) * sin(phi);
            double zl = c_ax * cos(theta);
            /* Transform to world frame: r = R · r_local */
            points[idx].x = pa.axes.m[0]*xl + pa.axes.m[1]*yl + pa.axes.m[2]*zl;
            points[idx].y = pa.axes.m[3]*xl + pa.axes.m[4]*yl + pa.axes.m[5]*zl;
            points[idx].z = pa.axes.m[6]*xl + pa.axes.m[7]*yl + pa.axes.m[8]*zl;
            idx++;
        }
    }
}

double inertia_ellipsoid_eval(const InertiaTensor *I, vec3 rho) {
    double x = rho.x, y = rho.y, z = rho.z;
    return I->Ixx*x*x + I->Iyy*y*y + I->Izz*z*z
         + 2.0*I->Ixy*x*y + 2.0*I->Ixz*x*z + 2.0*I->Iyz*y*z;
}

double radius_of_gyration(const InertiaTensor *I, double M, vec3 axis) {
    /* I_axis = n̂^T · I · n̂ */
    double nx = axis.x, ny = axis.y, nz = axis.z;
    double I_axis = I->Ixx*nx*nx + I->Iyy*ny*ny + I->Izz*nz*nz
                  + 2.0*I->Ixy*nx*ny + 2.0*I->Ixz*nx*nz + 2.0*I->Iyz*ny*nz;
    if (M < 1e-15) return 0.0;
    return sqrt(I_axis / M);
}
