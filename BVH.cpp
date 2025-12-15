#include "BVH.h"
#include <cmath>
#include <iostream>
#include <algorithm>

// ==========================================
// AABB 实现
// ==========================================
AABB::AABB() {
    double inf = 1e30;
    min = { inf, inf, inf };
    max = { -inf, -inf, -inf };
}

void AABB::expand(const Vec3& p) {
    min.x = std::fmin(min.x, p.x); min.y = std::fmin(min.y, p.y); min.z = std::fmin(min.z, p.z);
    max.x = std::fmax(max.x, p.x); max.y = std::fmax(max.y, p.y); max.z = std::fmax(max.z, p.z);
}

void AABB::expand(const AABB& box) {
    expand(box.min);
    expand(box.max);
}

Vec3 AABB::center() const {
    return (min + max) * 0.5;
}

bool AABB::intersect(const Vec3& origin, const Vec3& invDir, double t_max_limit) const {
    double tx1 = (min.x - origin.x) * invDir.x;
    double tx2 = (max.x - origin.x) * invDir.x;
    double tmin = std::fmin(tx1, tx2);
    double tmax = std::fmax(tx1, tx2);

    double ty1 = (min.y - origin.y) * invDir.y;
    double ty2 = (max.y - origin.y) * invDir.y;
    tmin = std::fmax(tmin, std::fmin(ty1, ty2));
    tmax = std::fmin(tmax, std::fmax(ty1, ty2));

    double tz1 = (min.z - origin.z) * invDir.z;
    double tz2 = (max.z - origin.z) * invDir.z;
    tmin = std::fmax(tmin, std::fmin(tz1, tz2));
    tmax = std::fmin(tmax, std::fmax(tz1, tz2));

    return tmax >= tmin && tmax > 0 && tmin < t_max_limit;
}

// ==========================================
// BVHNode 实现
// ==========================================
BVHNode::~BVHNode() {
    if (left) delete left;
    if (right) delete right;
}

// ==========================================
// 辅助：射线-三角形求交 (静态函数)
// ==========================================
static bool ray_tri_intersect(const Vec3& orig, const Vec3& dir, const Triangle& tri, double max_dist) {
    Vec3 v0 = tri.v0;
    Vec3 v1 = tri.v1;
    Vec3 v2 = tri.v2;

    Vec3 e1 = v1 - v0;
    Vec3 e2 = v2 - v0;
    Vec3 h = cross(dir, e2);
    double a = dot(e1, h);

    if (a > -EPSILON && a < EPSILON) return false; // 平行

    double f = 1.0 / a;
    Vec3 s = orig - v0;
    double u = f * dot(s, h);
    if (u < 0.0 || u > 1.0) return false;

    Vec3 q = cross(s, e1);
    double v = f * dot(dir, q);
    if (v < 0.0 || u + v > 1.0) return false;

    double t = f * dot(e2, q);
    return (t > 1e-4 && t < max_dist); // 1e-4 是防止自我遮挡的容差
}
static double ray_tri_intersect_t(const Vec3& orig, const Vec3& dir, const Triangle& tri, double max_dist) {
    Vec3 v0 = tri.v0; Vec3 v1 = tri.v1; Vec3 v2 = tri.v2;
    Vec3 e1 = v1 - v0; Vec3 e2 = v2 - v0; Vec3 h = cross(dir, e2);
    double a = dot(e1, h);
    if (a > -1e-8 && a < 1e-8) return -1.0;
    double f = 1.0 / a; Vec3 s = orig - v0; double u = f * dot(s, h);
    if (u < 0.0 || u > 1.0) return -1.0;
    Vec3 q = cross(s, e1); double v = f * dot(dir, q);
    if (v < 0.0 || u + v > 1.0) return -1.0;
    double t = f * dot(e2, q);
    if (t > 1e-4 && t < max_dist) return t;
    return -1.0;
}

// ==========================================
// BVHAccel 实现
// ==========================================
BVHAccel::BVHAccel() {}
BVHAccel::~BVHAccel() { if (root) delete root; }

void BVHAccel::build(const std::vector<ThermalNode>& nodes) {
    if (nodes.empty()) return;
    ref_nodes = &nodes;

    // 创建初始索引列表
    std::vector<int> indices(nodes.size());
    for (int i = 0; i < nodes.size(); ++i) indices[i] = i;

    // 清理旧树
    if (root) delete root;

    std::cout << "[BVH] Building acceleration structure for " << nodes.size() << " objects..." << std::endl;
    root = recursive_build(indices);
}

BVHNode* BVHAccel::recursive_build(std::vector<int>& indices) {
    BVHNode* node = new BVHNode();

    // 1. 计算当前所有物体的总包围盒
    AABB bounds;
    for (int idx : indices) {
        const auto& t_node = (*ref_nodes)[idx];
        for (const auto& tri : t_node.geometry_tris) {
            bounds.expand(tri.v0);
            bounds.expand(tri.v1);
            bounds.expand(tri.v2);
        }
    }
    node->bbox = bounds;

    // 2. 终止条件：物体少于等于2个，或者树太深（这里简单判断数量）
    if (indices.size() <= 2) {
        node->node_indices = indices;
        return node;
    }

    // 3. 寻找最长轴进行分割
    Vec3 extent = bounds.max - bounds.min;
    int axis = 0;
    if (extent.y > extent.x) axis = 1;
    if (extent.z > ((axis == 0) ? extent.x : extent.y)) axis = 2;

    // 4. 根据物体质心排序
    // 质心我们用该 Node 第一个三角形的质心近似，或者计算 Node 的平均质心
    // 这里简单用 node.centroid
    int mid = indices.size() / 2;
    std::nth_element(indices.begin(), indices.begin() + mid, indices.end(),
        [&](int a, int b) {
            const auto& nodeA = (*ref_nodes)[a];
            const auto& nodeB = (*ref_nodes)[b];
            if (axis == 0) return nodeA.centroid.x < nodeB.centroid.x;
            if (axis == 1) return nodeA.centroid.y < nodeB.centroid.y;
            return nodeA.centroid.z < nodeB.centroid.z;
        });

    // 5. 递归构建子节点
    std::vector<int> left_idxs(indices.begin(), indices.begin() + mid);
    std::vector<int> right_idxs(indices.begin() + mid, indices.end());

    node->left = recursive_build(left_idxs);
    node->right = recursive_build(right_idxs);

    return node;
}

bool BVHAccel::intersect_shadow(const Vec3& origin, const Vec3& dir, double max_dist, int self_id) const {
    if (!root) return false;
    // 预计算方向倒数，加速 AABB 检测
    Vec3 invDir = { 1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z };
    return recursive_intersect(root, origin, dir, invDir, max_dist, self_id);
}

bool BVHAccel::recursive_intersect(BVHNode* node, const Vec3& origin, const Vec3& dir, const Vec3& invDir, double max_dist, int self_id) const {
    // 1. 盒子检测
    if (!node->bbox.intersect(origin, invDir, max_dist)) return false;

    // 2. 叶子节点：遍历内部所有三角形
    if (node->is_leaf()) {
        for (int idx : node->node_indices) {
            if (idx == self_id) continue; // 忽略自己

            const auto& target = (*ref_nodes)[idx];
            for (const auto& tri : target.geometry_tris) {
                if (ray_tri_intersect(origin, dir, tri, max_dist)) {
                    return true; // 只要碰到一个，就是被遮挡
                }
            }
        }
        return false;
    }

    // 3. 中间节点：递归
    // 简单的优化：先查哪个盒子近？这里暂且按顺序查
    if (recursive_intersect(node->left, origin, dir, invDir, max_dist, self_id)) return true;
    if (recursive_intersect(node->right, origin, dir, invDir, max_dist, self_id)) return true;

    return false;
}

// 【MCRT 接口实现】
HitInfo BVHAccel::intersect_closest(const Vec3& origin, const Vec3& dir, double max_dist, int self_id) const {
    HitInfo closest;
    closest.has_hit = false;
    closest.t = max_dist;

    if (!root) return closest;
    Vec3 invDir = { 1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z };

    recursive_intersect_closest(root, origin, dir, invDir, max_dist, self_id, closest);
    return closest;
}

void BVHAccel::recursive_intersect_closest(BVHNode* node, const Vec3& origin, const Vec3& dir, const Vec3& invDir, double max_dist, int self_id, HitInfo& closest_hit) const {
    // 如果当前包围盒比已知的最近交点还远，直接跳过
    if (!node->bbox.intersect(origin, invDir, closest_hit.t)) return;

    if (node->is_leaf()) {
        for (int idx : node->node_indices) {
            if (idx == self_id) continue;
            const auto& target = (*ref_nodes)[idx];

            for (const auto& tri : target.geometry_tris) {
                double t = ray_tri_intersect_t(origin, dir, tri, closest_hit.t);

                if (t > 0.0 && t < closest_hit.t) {
                    closest_hit.has_hit = true;
                    closest_hit.t = t;
                    closest_hit.hit_node_index = idx;

                    // 【核心 TAITherm 逻辑：判断正反面】
                    Vec3 tri_normal = normalize(cross(tri.v1 - tri.v0, tri.v2 - tri.v0));
                    // 光线方向 dot 法线 < 0 说明是迎面击中 (Front)
                    if (dot(dir, tri_normal) < 0.0) {
                        closest_hit.hit_front_side = true;
                    }
                    else {
                        closest_hit.hit_front_side = false;
                    }
                }
            }
        }
    }
    else {
        // 中间节点：为了求最近，必须都查一遍
        recursive_intersect_closest(node->left, origin, dir, invDir, max_dist, self_id, closest_hit);
        recursive_intersect_closest(node->right, origin, dir, invDir, max_dist, self_id, closest_hit);
    }
}