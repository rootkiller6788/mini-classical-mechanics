/**
 * vec3.h — 三维欧几里得向量 (L1 Definitions / L3 Mathematical Structures)
 *
 * 所有牛顿力学量的基础算术类型。
 * 定理来源: Euclidean geometry axioms (Euclid's Elements, Book I)
 * 课程对标: MIT 8.012 §1.1-1.3 Vector algebra
 *            Goldstein §1.2 Vectors in Classical Mechanics
 *
 * 设计决策: 栈分配值语义, 零外部依赖 (仅 <math.h>)
 */

#ifndef VEC3_H
#define VEC3_H

#include <math.h>

/* M_PI is POSIX, not C11 standard; define if absent */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================
 * Vec3 — 三维向量 (栈分配, 值语义)
 * ================================================================ */
typedef struct {
    double x, y, z;
} Vec3;

/* 构造函数 */
Vec3 vec3_zero(void);
Vec3 vec3_make(double x, double y, double z);
Vec3 vec3_from_scalar(double s);

/* L1 基本运算 */
Vec3 vec3_add(Vec3 a, Vec3 b);
Vec3 vec3_sub(Vec3 a, Vec3 b);
Vec3 vec3_neg(Vec3 v);
Vec3 vec3_scale(Vec3 v, double s);
Vec3 vec3_div(Vec3 v, double s);

/* L3 内积与外积 */
double vec3_dot(Vec3 a, Vec3 b);
Vec3 vec3_cross(Vec3 a, Vec3 b);
double vec3_scalar_triple(Vec3 a, Vec3 b, Vec3 c);   /* a·(b×c) */
Vec3 vec3_vector_triple(Vec3 a, Vec3 b, Vec3 c);      /* a×(b×c) = b(a·c)-c(a·b) */

/* L3 范数与度量 */
double vec3_norm(Vec3 v);
double vec3_norm2(Vec3 v);
Vec3 vec3_normalize(Vec3 v);
double vec3_distance(Vec3 a, Vec3 b);
double vec3_angle(Vec3 a, Vec3 b);

/* L3 向量分解 */
Vec3 vec3_project(Vec3 a, Vec3 b);
Vec3 vec3_reject(Vec3 a, Vec3 b);
Vec3 vec3_reflect(Vec3 v, Vec3 normal);
Vec3 vec3_lerp(Vec3 a, Vec3 b, double t);

/* L1 坐标系转换 */
typedef struct { double r, theta, phi; } Spherical;
typedef struct { double rho, phi, z; } Cylindrical;

Spherical vec3_to_spherical(Vec3 v);
Vec3    vec3_from_spherical(double r, double theta, double phi);
Cylindrical vec3_to_cylindrical(Vec3 v);
Vec3    vec3_from_cylindrical(double rho, double phi, double z);

/* L1 坐标系基向量 */
Vec3 spherical_radial_unit(double theta, double phi);
Vec3 spherical_polar_unit(double theta, double phi);
Vec3 spherical_azimuthal_unit(double phi);
Vec3 cylindrical_radial_unit(double phi);
Vec3 cylindrical_angular_unit(double phi);

/* L3: Basis decomposition and rotations */
void vec3_decompose_basis(Vec3 v, Vec3 e1, Vec3 e2, Vec3 e3,
                          double *c1, double *c2, double *c3);
Vec3 vec3_rotate_axis_angle(Vec3 v, Vec3 axis, double angle);

#endif /* VEC3_H */
