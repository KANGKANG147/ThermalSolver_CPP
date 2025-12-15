#pragma once // 防止头文件被重复包含
#include <cmath>
#include <iostream>

const double EPSILON = 1e-6;
const double SIGMA = 5.67e-8;

struct Vec3 {
    double x, y, z;
    Vec3 operator+(const Vec3& v) const { return { x + v.x, y + v.y, z + v.z }; }
    Vec3 operator-(const Vec3& v) const { return { x - v.x, y - v.y, z - v.z }; }
    Vec3 operator*(double s) const { return { x * s, y * s, z * s }; }
    Vec3 operator/(double s) const { return { x / s, y / s, z / s }; }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 normalize(const Vec3& v) { double len = std::sqrt(dot(v, v)); return len > EPSILON ? v / len : v; }
inline Vec3 cross(const Vec3 & a, const Vec3 & b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }

// --- MCRT 新增数学工具 ---

// 0到1随机数
inline double random_double() {
    return rand() / (RAND_MAX + 1.0);
}

// 建立局部坐标系 (Orthonormal Basis)
inline void build_frame(const Vec3& n, Vec3& u, Vec3& v) {
    if (std::abs(n.x) > 0.9) v = { 0, 1, 0 };
    else v = { 1, 0, 0 };
    u = normalize(cross(v, n));
    v = cross(n, u);
}

// 余弦加权半球采样 (Cosine Weighted Hemisphere Sample)
inline Vec3 sample_hemisphere(const Vec3& normal) {
    double r1 = random_double();
    double r2 = random_double();
    double phi = 2.0 * 3.14159 * r1;
    double r = std::sqrt(r2);

    double x = r * std::cos(phi);
    double y = r * std::sin(phi);
    double z = std::sqrt(1.0 - r2);

    Vec3 u, v;
    build_frame(normal, u, v);

    return {
        x * u.x + y * v.x + z * normal.x,
        x * u.y + y * v.y + z * normal.y,
        x * u.z + y * v.z + z * normal.z
    };
}