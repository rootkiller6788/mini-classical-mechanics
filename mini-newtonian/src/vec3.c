/**
 * vec3.c — 三维欧几里得向量实现
 *
 * 每个函数实现一个独立几何/代数知识点。
 * 复杂度: 所有函数 O(1), 仅涉及标量浮点运算。
 */
#include "vec3.h"
#include <math.h>
#include <float.h>

/* ---------- 构造函数 ---------- */

Vec3 vec3_zero(void) {
    Vec3 v = {0.0, 0.0, 0.0};
    return v;
}

Vec3 vec3_make(double x, double y, double z) {
    Vec3 v = {x, y, z};
    return v;
}

Vec3 vec3_from_scalar(double s) {
    Vec3 v = {s, s, s};
    return v;
}

/* ---------- L1 基本运算 ---------- */

Vec3 vec3_add(Vec3 a, Vec3 b) {
    Vec3 r = {a.x + b.x, a.y + b.y, a.z + b.z};
    return r;
}

Vec3 vec3_sub(Vec3 a, Vec3 b) {
    Vec3 r = {a.x - b.x, a.y - b.y, a.z - b.z};
    return r;
}

Vec3 vec3_neg(Vec3 v) {
    Vec3 r = {-v.x, -v.y, -v.z};
    return r;
}

Vec3 vec3_scale(Vec3 v, double s) {
    Vec3 r = {v.x * s, v.y * s, v.z * s};
    return r;
}

Vec3 vec3_div(Vec3 v, double s) {
    Vec3 r = {v.x / s, v.y / s, v.z / s};
    return r;
}

/* ---------- L3 内积与外积 ---------- */

double vec3_dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 vec3_cross(Vec3 a, Vec3 b) {
    Vec3 r = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return r;
}

/**
 * 标量三重积 a·(b×c)
 * 几何意义: 由 a,b,c 张成的平行六面体的有符号体积
 * 知识点: 行列式性质 — 偶置换不变号, 奇置换变号
 */
double vec3_scalar_triple(Vec3 a, Vec3 b, Vec3 c) {
    return vec3_dot(a, vec3_cross(b, c));
}

/**
 * 向量三重积 a×(b×c) = b(a·c) - c(a·b)
 * 知识点: Lagrange 恒等式 / BAC-CAB 法则
 * 定理来源: Gibbs vector calculus (1881), §50
 */
Vec3 vec3_vector_triple(Vec3 a, Vec3 b, Vec3 c) {
    double ac = vec3_dot(a, c);
    double ab = vec3_dot(a, b);
    Vec3 r = {
        b.x * ac - c.x * ab,
        b.y * ac - c.y * ab,
        b.z * ac - c.z * ab
    };
    return r;
}

/* ---------- L3 范数与度量 ---------- */

double vec3_norm(Vec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

double vec3_norm2(Vec3 v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

Vec3 vec3_normalize(Vec3 v) {
    double n = vec3_norm(v);
    if (n < 1e-15) {
        return vec3_zero();
    }
    return vec3_div(v, n);
}

double vec3_distance(Vec3 a, Vec3 b) {
    return vec3_norm(vec3_sub(a, b));
}

/**
 * 向量夹角 (L3: 余弦定理在向量空间的推广)
 * θ = acos( (a·b) / (|a||b|) )
 * 值域 [0, π]
 */
double vec3_angle(Vec3 a, Vec3 b) {
    double na = vec3_norm(a);
    double nb = vec3_norm(b);
    if (na < 1e-15 || nb < 1e-15) return 0.0;
    double cos_theta = vec3_dot(a, b) / (na * nb);
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    return acos(cos_theta);
}

/* ---------- L3 向量分解 ---------- */

/**
 * 向量投影: proj_b(a) = (a·b̂) b̂
 * 知识点: Gram-Schmidt 正交化的基本步骤
 * 定理来源: 希尔伯特空间中的投影定理
 */
Vec3 vec3_project(Vec3 a, Vec3 b) {
    double nb2 = vec3_norm2(b);
    if (nb2 < 1e-30) return vec3_zero();
    return vec3_scale(b, vec3_dot(a, b) / nb2);
}

/**
 * 正交补: rej_b(a) = a - proj_b(a)
 * 满足 a = proj_b(a) + rej_b(a) 且 rej_b(a)·b = 0
 */
Vec3 vec3_reject(Vec3 a, Vec3 b) {
    return vec3_sub(a, vec3_project(a, b));
}

/**
 * 镜面反射: 以 normal 为法向的平面反射
 * v_ref = v - 2(v·n̂)n̂
 * 知识点: Householder 变换的几何原型
 */
Vec3 vec3_reflect(Vec3 v, Vec3 normal) {
    Vec3 n = vec3_normalize(normal);
    double d = 2.0 * vec3_dot(v, n);
    return vec3_sub(v, vec3_scale(n, d));
}

/**
 * 线性插值: lerp(a, b, t) = (1-t)a + tb
 * t=0 → a, t=1 → b
 */
Vec3 vec3_lerp(Vec3 a, Vec3 b, double t) {
    Vec3 r = {
        a.x + t * (b.x - a.x),
        a.y + t * (b.y - a.y),
        a.z + t * (b.z - a.z)
    };
    return r;
}

/* ---------- L1 坐标系转换 ---------- */

/**
 * 笛卡尔 → 球坐标 (r, θ, φ)
 * θ = polar angle from +z (colatitude), φ = azimuthal angle from +x
 * 知识点: 坐标图卡 (chart) 的局部同胚
 */
Spherical vec3_to_spherical(Vec3 v) {
    Spherical s;
    s.r = vec3_norm(v);
    s.theta = (s.r > 0.0) ? acos(v.z / s.r) : 0.0;
    s.phi = atan2(v.y, v.x);
    return s;
}

Vec3 vec3_from_spherical(double r, double theta, double phi) {
    double st = sin(theta);
    Vec3 v = {
        r * st * cos(phi),
        r * st * sin(phi),
        r * cos(theta)
    };
    return v;
}

/**
 * 笛卡尔 → 柱坐标 (ρ, φ, z)
 */
Cylindrical vec3_to_cylindrical(Vec3 v) {
    Cylindrical c;
    c.rho = sqrt(v.x * v.x + v.y * v.y);
    c.phi = atan2(v.y, v.x);
    c.z = v.z;
    return c;
}

Vec3 vec3_from_cylindrical(double rho, double phi, double z) {
    Vec3 v = {rho * cos(phi), rho * sin(phi), z};
    return v;
}

/* ---------- L1 坐标系基向量 ---------- */

Vec3 spherical_radial_unit(double theta, double phi) {
    double st = sin(theta);
    Vec3 v = {st * cos(phi), st * sin(phi), cos(theta)};
    return v;
}

Vec3 spherical_polar_unit(double theta, double phi) {
    double ct = cos(theta);
    Vec3 v = {ct * cos(phi), ct * sin(phi), -sin(theta)};
    return v;
}

Vec3 spherical_azimuthal_unit(double phi) {
    Vec3 v = {-sin(phi), cos(phi), 0.0};
    return v;
}

Vec3 cylindrical_radial_unit(double phi) {
    Vec3 v = {cos(phi), sin(phi), 0.0};
    return v;
}

Vec3 cylindrical_angular_unit(double phi) {
    Vec3 v = {-sin(phi), cos(phi), 0.0};
    return v;
}


/* L3: Gram-Schmidt orthogonalization for three vectors
 * Given an orthonormal basis {e1, e2, e3}, express v in this basis.
 * v = (v.e1)e1 + (v.e2)e2 + (v.e3)e3
 * Knowledge point: basis decomposition, coordinate-free form
 */
void vec3_decompose_basis(Vec3 v, Vec3 e1, Vec3 e2, Vec3 e3,
                          double *c1, double *c2, double *c3) {
    *c1 = vec3_dot(v, e1);
    *c2 = vec3_dot(v, e2);
    *c3 = vec3_dot(v, e3);
}

/* L3: Rotation about an arbitrary axis (Rodrigues rotation formula)
 * v_rot = v*cos(a) + (k x v)*sin(a) + k*(k.v)*(1-cos(a))
 * where k is the unit rotation axis and a is the angle.
 * Theorem: Rodrigues (1840) - SO(3) rotation formula
 */
Vec3 vec3_rotate_axis_angle(Vec3 v, Vec3 axis, double angle) {
    Vec3 k = vec3_normalize(axis);
    double c = cos(angle);
    double s = sin(angle);
    Vec3 term1 = vec3_scale(v, c);
    Vec3 term2 = vec3_scale(vec3_cross(k, v), s);
    Vec3 term3 = vec3_scale(k, vec3_dot(k, v) * (1.0 - c));
    return vec3_add(term1, vec3_add(term2, term3));
}
