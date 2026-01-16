#pragma once // 防止头文件被重复包含
#include <cmath>
#include <iostream>
#include <random> // 引入随机数库

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
//计算向量长度
inline double length(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalize(const Vec3& v) { double len = std::sqrt(dot(v, v)); return len > EPSILON ? v / len : v; }
inline Vec3 cross(const Vec3 & a, const Vec3 & b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }

// --- MCRT 新增数学工具 ---

// 0到1随机数
inline double random_double() {
    // static thread_local 保证生成器只初始化一次，且每个线程一份
    static thread_local std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
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

// 用于 QMC 采样的 Halton 序列生成器 (基数为 base)
// index: 采样点的序号 (1, 2, 3...)
inline double halton_sequence(int index, int base) {
    double f = 1.0;
    double r = 0.0;
    while (index > 0) {
        f = f / (double)base;
        r = r + f * (double)(index % base);
        index = index / base;
    }
    return r;
}

// 2. QMC 圆锥采样辅助函数
// 输入:
//   central_dir: 太阳中心方向
//   theta_max_rad: 圆锥半角 (弧度)
//   u, v: 两个 [0,1] 的随机数 (来自 Halton 序列)
static Vec3 sample_cone_deterministic(const Vec3& central_dir, double theta_max_rad, double u, double v) {
    // 1. 构建局部坐标系 (Frame)
    // 假设 MathUtils.h 中已有 build_frame，如果没有，这里内联实现一个
    Vec3 w = normalize(central_dir);
    Vec3 p, q;

    // 简单的构建正交基方法
    if (std::abs(w.x) > 0.9) p = { 0, 1, 0 };
    else p = { 1, 0, 0 };
    q = normalize(cross(p, w)); // q = p x w
    p = cross(w, q);            // p = w x q (注意顺序)

    // 2. 在单位圆锥内生成局部坐标 (Uniform sampling on spherical cap)
    // z 轴朝向 w
    // cos_theta 在 [cos(theta_max), 1] 之间均匀分布
    double cos_theta_max = std::cos(theta_max_rad);
    double z = 1.0 + u * (cos_theta_max - 1.0); // map u(0~1) to [1 ~ cos_max] -> z

    double sin_theta = std::sqrt(std::max(0.0, 1.0 - z * z));
    double phi = 2.0 * 3.14159265358979323846 * v;

    double x = sin_theta * std::cos(phi);
    double y = sin_theta * std::sin(phi);

    // 3. 转换到世界坐标: result = x*p + y*q + z*w
    return {
        x * p.x + y * q.x + z * w.x,
        x * p.y + y * q.y + z * w.y,
        x * p.z + y * q.z + z * w.z
    };
}