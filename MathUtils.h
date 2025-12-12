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
