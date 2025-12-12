#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>
#include "MathUtils.h" // 引用你之前的数学库
#include "CoreTypes.h" // 引用三角形定义

// 轴对齐包围盒(Axis-Aligned Bounding Box)
struct AABB {
    Vec3 min, max;

    AABB();// 初始化为一个无效的盒子

    // 将一个点加入盒子，扩展盒子范围
    void expand(const Vec3& p);

    // 将另一个盒子合并进来
    void expand(const AABB& box);

    Vec3 center() const;

    // 光线与盒子求交检测 (Slab Method 优化版)
    // 如果相交，返回 true，dist 存储到盒子的最近距离
    bool intersect(const Vec3& origin, const Vec3& invDir, double t_max_limit) const;
};

// BVH 树节点
struct BVHNode {
    AABB bbox;
    BVHNode* left = nullptr;
    BVHNode* right = nullptr;
    std::vector<int> node_indices; // 叶子节点存储的是 ThermalNode 在数组中的下标

    bool is_leaf() const { return !node_indices.empty(); }
    ~BVHNode(); // 析构函数负责清理内存
};

// BVH 加速器主类
class BVHAccel {
public:
    BVHAccel();
    ~BVHAccel();

    // 构建 BVH 树
    void build(const std::vector<ThermalNode>& nodes);

    // 阴影遮挡查询：返回 true 表示被遮挡
    // origin: 射线起点, dir: 方向, max_dist: 最大距离, self_id: 发射者ID(防止自遮挡)
    bool intersect_shadow(const Vec3& origin, const Vec3& dir, double max_dist, int self_id) const;

private:
    BVHNode* root = nullptr;
    const std::vector<ThermalNode>* ref_nodes = nullptr; // 持有外部数据的指针

    // 递归构建
    BVHNode* recursive_build(std::vector<int>& indices);

    // 递归查询
    bool recursive_intersect(BVHNode* node, const Vec3& origin, const Vec3& dir, const Vec3& invDir, double max_dist, int self_id) const;
};