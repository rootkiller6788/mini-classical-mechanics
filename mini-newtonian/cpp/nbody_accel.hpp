#pragma once
// nbody_accel.hpp — N体引力加速度热内核 (C++17)
// Julia 通过 ccall 调用 extern "C" 接口

#include <cstddef>
#include <cstdint>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * 计算 N 体系统所有粒子的引力加速度
 *
 * @param n         粒子数
 * @param masses    质量数组 [n] (kg)
 * @param pos_x     位置 x 分量 [n] (m)
 * @param pos_y     位置 y 分量 [n] (m)
 * @param pos_z     位置 z 分量 [n] (m)
 * @param acc_x     输出加速度 x 分量 [n] (m/s²)
 * @param acc_y     输出加速度 y 分量 [n] (m/s²)
 * @param acc_z     输出加速度 z 分量 [n] (m/s²)
 * @param G         万有引力常数
 * @param softening 软化因子（避免奇点，对直接求和法必要）
 */
void nbody_accel_direct(
    int64_t n,
    const double* masses,
    const double* pos_x, const double* pos_y, const double* pos_z,
    double* acc_x, double* acc_y, double* acc_z,
    double G, double softening
);

/**
 * 单粒子受 N-1 个源粒子的引力加速度
 * 用于受限制多体问题（大质量源不移动）
 *
 * @param n_sources 源粒子数
 * @param source_masses 源质量数组
 * @param source_x/y/z 源位置
 * @param target_x/y/z 测试粒子位置
 * @param ax/ay/az    输出加速度
 * @param G          引力常数
 * @param softening  软化因子
 */
void nbody_test_particle_accel(
    int64_t n_sources,
    const double* source_masses,
    const double* source_x, const double* source_y, const double* source_z,
    double target_x, double target_y, double target_z,
    double* ax, double* ay, double* az,
    double G, double softening
);

#ifdef __cplusplus
}
#endif
