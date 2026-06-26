// nbody_accel.cpp — N体引力加速度热内核
// 编译: g++ -std=c++17 -O3 -march=native -shared -fPIC -o libnbody.so nbody_accel.cpp
// Julia 调用: ccall((:nbody_accel_direct, "libnbody"), Cvoid, ...)

#include "nbody_accel.hpp"
#include <cmath>
#include <cstring>
#include <algorithm>

// 内联快速倒数平方根（可选启用）
#ifdef USE_FAST_RSQRT
static inline float fast_rsqrt(float x) {
    // Quake III 算法
    float xhalf = 0.5f * x;
    int i = *(int*)&x;
    i = 0x5f3759df - (i >> 1);
    x = *(float*)&i;
    return x * (1.5f - xhalf * x * x);
}
#endif

void nbody_accel_direct(
    int64_t n,
    const double* masses,
    const double* pos_x, const double* pos_y, const double* pos_z,
    double* acc_x, double* acc_y, double* acc_z,
    double G, double softening)
{
    // 归零输出
    std::memset(acc_x, 0, n * sizeof(double));
    std::memset(acc_y, 0, n * sizeof(double));
    std::memset(acc_z, 0, n * sizeof(double));

    const double soft2 = softening * softening;

    // O(N²) 直接求和，对称性利用：F_ij = -F_ji
    for (int64_t i = 0; i < n; ++i) {
        const double xi = pos_x[i];
        const double yi = pos_y[i];
        const double zi = pos_z[i];
        double axi = 0.0, ayi = 0.0, azi = 0.0;

        for (int64_t j = i + 1; j < n; ++j) {
            double dx = xi - pos_x[j];
            double dy = yi - pos_y[j];
            double dz = zi - pos_z[j];

            double dist2 = dx*dx + dy*dy + dz*dz + soft2;
            double dist  = std::sqrt(dist2);
            double dist3 = dist2 * dist;  // dist³

            double fx = G * dx / dist3;
            double fy = G * dy / dist3;
            double fz = G * dz / dist3;

            // a_i += m_j * f_ij
            double mj_times_fx = masses[j] * fx;
            double mj_times_fy = masses[j] * fy;
            double mj_times_fz = masses[j] * fz;

            axi += mj_times_fx;
            ayi += mj_times_fy;
            azi += mj_times_fz;

            // a_j += -m_i * f_ij (牛顿第三定律)
            double mi_times_fx = masses[i] * fx;
            double mi_times_fy = masses[i] * fy;
            double mi_times_fz = masses[i] * fz;

            acc_x[j] -= mi_times_fx;
            acc_y[j] -= mi_times_fy;
            acc_z[j] -= mi_times_fz;
        }

        acc_x[i] += axi;
        acc_y[i] += ayi;
        acc_z[i] += azi;
    }
}

void nbody_test_particle_accel(
    int64_t n_sources,
    const double* source_masses,
    const double* source_x, const double* source_y, const double* source_z,
    double target_x, double target_y, double target_z,
    double* ax, double* ay, double* az,
    double G, double softening)
{
    *ax = 0.0;
    *ay = 0.0;
    *az = 0.0;

    const double soft2 = softening * softening;

    for (int64_t j = 0; j < n_sources; ++j) {
        double dx = target_x - source_x[j];
        double dy = target_y - source_y[j];
        double dz = target_z - source_z[j];

        double dist2 = dx*dx + dy*dy + dz*dz + soft2;
        double dist  = std::sqrt(dist2);
        double dist3 = dist2 * dist;

        double factor = G * source_masses[j] / dist3;

        *ax -= factor * dx;
        *ay -= factor * dy;
        *az -= factor * dz;
    }
}
