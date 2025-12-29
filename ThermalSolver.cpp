#include "ThermalSolver.h"
#include <iostream>
#include <map>
#include <algorithm>
#include <omp.h>

// 辅助结构体
struct VertexKey {
    long long x, y, z;
    bool operator<(const VertexKey& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }
};

struct EdgeKey {
    int v1, v2;
    bool operator<(const EdgeKey& other) const {
        if (v1 != other.v1) return v1 < other.v1;
        return v2 < other.v2;
    }
};

VertexKey to_key(const Vec3& v) {
    const double SCALE = 1000.0;
    return { (long long)std::round(v.x * SCALE), (long long)std::round(v.y * SCALE), (long long)std::round(v.z * SCALE) };
}

// 辅助函数：计算点 p 到线段 (v1-v2) 的最短距离的平方
// 用于计算 Centroid 到 Edge 的垂直距离
double point_segment_distance_sq(const Vec3& p, const Vec3& v1, const Vec3& v2) {
    Vec3 ab = v2 - v1;
    Vec3 ap = p - v1;
    double len_sq = dot(ab, ab);
    if (len_sq < 1e-8) return dot(ap, ap); // v1 == v2 (退化边)

    // 投影参数 t
    double t = dot(ap, ab) / len_sq;

    // 限制 t 在 [0, 1] 范围内 (限制在线段上)
    // 对于热传导，通常我们需要的是“到无限长直线的垂直距离”还是“到线段的距离”？
    // 物理上，如果质心投影在边之外，说明网格扭曲很大。
    // 为了稳健，我们通常允许投影点在线段延伸线上，或者限制在端点。
    // 这里使用限制在线段上的逻辑，保证距离是物理真实的“最近路程”。
    if (t < 0.0) t = 0.0;
    else if (t > 1.0) t = 1.0;

    Vec3 closest = v1 + (ab * t);
    Vec3 d = p - closest;
    return dot(d, d);
}

// ThermalSolver.cpp

void ThermalSolver::ensure_scene_scale() {
    // 如果已经计算过 (大于0)，直接返回，避免重复计算
    if (scene_scale > 1e-9) return;

    if (nodes.empty()) {
        scene_scale = 1.0; // 默认值防崩溃
        return;
    }

    Vec3 min_b = { 1e30, 1e30, 1e30 };
    Vec3 max_b = { -1e30, -1e30, -1e30 };

    for (const auto& node : nodes) {
        if (node.centroid.x < min_b.x) min_b.x = node.centroid.x;
        if (node.centroid.y < min_b.y) min_b.y = node.centroid.y;
        if (node.centroid.z < min_b.z) min_b.z = node.centroid.z;

        if (node.centroid.x > max_b.x) max_b.x = node.centroid.x;
        if (node.centroid.y > max_b.y) max_b.y = node.centroid.y;
        if (node.centroid.z > max_b.z) max_b.z = node.centroid.z;
    }

    Vec3 diag = max_b - min_b;
    // 计算对角线长度
    scene_scale = std::sqrt(dot(diag, diag));

    // 防止单点或极小模型导致 scale 为 0
    if (scene_scale < 1e-6) scene_scale = 1.0;

    std::cout << "[Scene] Computed Scene Scale: " << scene_scale << " m" << std::endl;
}

void ThermalSolver::build_topology() {
    std::cout << "[Topology] Building lateral connections (Advanced Quad Support)..." << std::endl;

    std::vector<Vec3> unique_verts;

    // 1. 顶点焊接 (保持不变)
    std::map<VertexKey, int> vert_map;
    auto get_vert_id = [&](const Vec3& v) -> int {
        VertexKey key = to_key(v);
        if (vert_map.find(key) == vert_map.end()) {
            int new_id = (int)unique_verts.size();
            vert_map[key] = new_id;
            unique_verts.push_back(v);
            return new_id;
        }
        return vert_map[key];
        };

    // 2. 全局边注册 (Edge Registry)
    // 记录每一条边被哪些 Node 引用过
    // Map: EdgeKey -> List of Node Indices
    std::map<EdgeKey, std::vector<int>> edge_registry;

    for (int i = 0; i < nodes.size(); ++i) {
        // ★★★ 修正点1：遍历该 Node 下所有的三角形 ★★★
        for (const auto& tri : nodes[i].geometry_tris) {
            int id0 = get_vert_id(tri.v0);
            int id1 = get_vert_id(tri.v1);
            int id2 = get_vert_id(tri.v2);

            int v[3] = { id0, id1, id2 };

            // 注册该三角形的三条边
            for (int k = 0; k < 3; ++k) {
                int a = v[k];
                int b = v[(k + 1) % 3];
                if (a > b) std::swap(a, b); // 保证从小到大排序，确保 EdgeKey 唯一

                auto& list = edge_registry[{a, b}];
                if (list.empty() || list.back() != i) {
                    list.push_back(i);
                }
                // 注意：这里 push 进去的是 Node 的索引 i，而不是三角形的索引
            }
        }
    }

    // 3. 建立连接并过滤内部边
    int links_count = 0;

    for (auto& kv : edge_registry) {
        const EdgeKey& edge_key = kv.first;
        const std::vector<int>& shared_nodes = kv.second;

        // 情况 A: 这条边只出现过 1 次 -> 它是模型的物理边界（边缘），没有邻居。
        if (shared_nodes.size() != 2) continue;

        // 情况 B: 这条边出现了 2 次（最常见的情况）
        // 这里有两种可能：
        // 1. 两个不同的 Node 共享这条边 -> 是我们要的邻居！
        // 2. 同一个 Node 内部的两个三角形共享这条边（对角线） -> 这是内部边，要忽略！

        // 我们这里只处理 shared_nodes 的前两个。
        // (如果有 >2 个，说明是非流形几何，比如三片叶子连在一条轴上，暂不处理)
        int n1 = shared_nodes[0];
        int n2 = shared_nodes[1];

        // ★★★ 修正点2：过滤内部对角线 ★★★
        if (n1 == n2) {
            continue;
        }

        ThermalNode& node1 = nodes[n1];
        ThermalNode& node2 = nodes[n2];

        // 1. 计算真实的接触边长 (保留你的认可)
        const Vec3& p1 = unique_verts[edge_key.v1];
        const Vec3& p2 = unique_verts[edge_key.v2];
        double L_edge_real = std::sqrt(dot(p1 - p2, p1 - p2));
        if (L_edge_real < 1e-5) continue;

        // --- 2. ★核心修正★：非均匀网格加权 (Centroid Bias) ---
        // 不再使用 dist = distance(c1, c2)
        // 而是使用 dist = dist(c1, edge) + dist(c2, edge)
        // 这自动处理了大小网格不匹配的问题，也隐含处理了角度问题

        double d1_sq = point_segment_distance_sq(node1.centroid, p1, p2);
        double d2_sq = point_segment_distance_sq(node2.centroid, p1, p2);

        double d1 = std::sqrt(d1_sq);
        double d2 = std::sqrt(d2_sq);

        // 有效传热距离 (Total Effective Path Length)
        double dist_effective = d1 + d2;

        // 安全钳位：防止距离过小 (比如质心恰好在边上)
        if (dist_effective < 1e-4) dist_effective = 1e-4;

        // 3. 接触面积
        double t_avg = (node1.thickness + node2.thickness) / 2.0;
        double A_contact = L_edge_real * t_avg;

        // ★★★ 关键升级：异质材料导热系数 (调和平均) ★★★
        // 如果 k1=45(钢), k2=0.04(绝热), 结果约为 0.08 (接近绝热), 这才是对的。
        // 如果用算术平均 (45+0.04)/2 = 22.5, 绝热层就失效了。
        double k1 = node1.k_mat;
        double k2 = node2.k_mat;
        double k_interface = 0.0;

        if (k1 < EPSILON || k2 < EPSILON) {
            k_interface = 0.0; // 只要有一方不导热，整体就不导热
        }
        else {
            // 调和平均公式： 2*k1*k2 / (k1+k2)
            k_interface = (2.0 * k1 * k2) / (k1 + k2);
        }

        double conductance = (k_interface * A_contact) / dist_effective;

        node1.neighbors.push_back({ n2, conductance });
        node2.neighbors.push_back({ n1, conductance });
        links_count++;
    }

    std::cout << "[Topology] Created " << links_count << " lateral links." << std::endl;
    
    // --- 构建 BVH ---
    bvh.build(nodes);
}

// ==========================================
// 2. MCRT 核心实现 把一个 ThermalNode 视为两个独立的辐射源
// ==========================================
void ThermalSolver::calculate_view_factors(int samples) {
    // --- 1. 计算特征码并尝试加载缓存 ---
    size_t checksum = compute_geometry_checksum(samples);
    std::string cache_file = "vf_cache_" + std::to_string(checksum) + ".bin";

    if (load_vf_cache(cache_file, samples)) {
        // 如果加载成功，直接返回，跳过漫长的计算！
        return;
    }

    std::cout << "[MCRT] Calculating View Factors (" << samples << " rays/node)..." << std::endl;

    // 可以在这里打印一下最大线程数
    std::cout << "Max Threads available: " << omp_get_max_threads() << std::endl;

    // 确保场景尺度已计算
    ensure_scene_scale();
    // [优化1] 基于场景尺度的自适应 Bias
    // 取场景大小的万分之一作为偏移量，且不小于 1e-5 米
    const double ADAPTIVE_BIAS = std::max(1e-5, scene_scale * 1e-4);

#pragma omp parallel for 
    for (int i = 0; i < nodes.size(); ++i) {
        ThermalNode& node = nodes[i];

        // 计算节点的【正面】辐射
        // 临时统计 map: key -> pair<目标ID, 目标正反>, value -> 命中次数
        std::map<std::pair<int, bool>, int> hits_front;
        int sky_hits_front = 0;
        Vec3 origin_f = node.centroid + node.normal * ADAPTIVE_BIAS;// 发射点：质心沿法线向外偏移一点点

        for (int s = 0; s < samples; ++s) {
            // 1. 生成基于法线(Normal)的随机射线
            Vec3 dir = sample_hemisphere(node.normal);
            // 查找最近的物体，距离上限设大一点
            HitInfo hit = bvh.intersect_closest(origin_f, dir, 1e20, i);
            // 击中了别的物体 (target)
            if (hit.has_hit) hits_front[{hit.hit_node_index, hit.hit_front_side}]++;
            // 击中天空
            else sky_hits_front++;
        }

        // 归一化并存入 rad_links_front
        node.rad_links_front.clear();
        for (auto const& entry : hits_front) {
            auto const& map_key = entry.first;
            int count = entry.second;
            double vf = (double)count / samples;
            if (vf > 0.001) node.rad_links_front.push_back({ map_key.first, map_key.second, vf });
        }
        node.vf_sky_front = (double)sky_hits_front / samples;

        // 计算节点的【背面】辐射
        if (node.bc_back.type != CONV_INSULATED) {
            std::map<std::pair<int, bool>, int> hits_back;
            int sky_hits_back = 0;
            Vec3 origin_b = node.centroid - node.normal * ADAPTIVE_BIAS;// 向内偏移
            Vec3 normal_b = node.normal * -1.0;// 反向法线

            for (int s = 0; s < samples; ++s) {
                // 1. 生成基于反向法线(-Normal)的随机射线
                Vec3 dir = sample_hemisphere(normal_b);
                HitInfo hit = bvh.intersect_closest(origin_b, dir, 1e20, i);
                if (hit.has_hit) hits_back[{hit.hit_node_index, hit.hit_front_side}]++;
                else sky_hits_back++;
            }

            node.rad_links_back.clear();
            for (auto const& entry : hits_back) {
                auto const& map_key = entry.first;
                int count = entry.second;
                double vf = (double)count / samples;
                if (vf > 0.001) node.rad_links_back.push_back({ map_key.first, map_key.second, vf });
            }
            node.vf_sky_back = (double)sky_hits_back / samples;
        }
        else {
            node.vf_sky_back = 0.0;
        }
    }
    std::cout << "[MCRT] Done." << std::endl;

    // --- 2. 计算完成后保存缓存 ---
    save_vf_cache(cache_file, samples);
}

void ThermalSolver::update_shadows(const Vec3& sun_dir) {
    if (nodes.empty()) return;

    // --- 2. 设定自适应 Bias(防自遮挡) ---
    // 经验值：取场景尺度的万分之一
    // 例如：100米的大楼 -> Bias = 1cm
    //       0.1米的手机 -> Bias = 0.01mm
    // 同时设置一个极小值 (如 1e-6) 防止 Bias 变为 0
    ensure_scene_scale();
    double adaptive_bias = std::max(1e-6, scene_scale * 1e-4);

    Vec3 ray_dir = normalize(sun_dir); 
    const double MAX_DIST = 1.0e20; // 太阳视为无限远
    
    // 开启 OpenMP 并行计算，BVH 是只读的，所以并行是安全的
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nodes.size(); ++i) {
        ThermalNode& receiver = nodes[i];
        Vec3 face_normal = receiver.normal;
        if (dot(receiver.normal, ray_dir) < 0.0) face_normal = receiver.normal * -1.0;
        Vec3 origin = receiver.centroid + (face_normal * adaptive_bias);
        
        // --- BVH 加速查询 ---
        // 参数：起点, 方向, 最大距离, 排除的ID(自己)
        bool blocked = bvh.intersect_shadow(origin, ray_dir, MAX_DIST, i);
        receiver.shadow_factor = blocked ? 0.0 : 1.0;
    }
}

// ==========================================
// 辅助：辐射线性化系数计算
// ==========================================
double ThermalSolver::calc_h_rad(double T_surf_K, double T_env_K, double epsilon) {
    // 防止除零错误
    if (std::abs(T_surf_K - T_env_K) < 0.1) {
        return 4.0 * SIGMA * epsilon * std::pow(T_surf_K, 3);
    }
    return SIGMA * epsilon * (std::pow(T_surf_K, 2) + std::pow(T_env_K, 2)) * (T_surf_K + T_env_K);
}

// ==========================================
// 对流系数 
// ==========================================
// 返回 pair: {h_value, T_fluid}
std::pair<double, double> ThermalSolver::get_convection_params(const ThermalNode& node, bool is_front, const WeatherData& w) {
    const ConvectionBC& bc = is_front ? node.bc_front : node.bc_back;
    double T_surf = is_front ? node.T_front : node.T_back;

    if (bc.type == CONV_INSULATED) return { 0.0, 0.0 };

    if (bc.type == CONV_FIXED_H_T) {
        double T_f = (std::abs(bc.fixed_fluid_T) < 0.001) ? w.air_temp : bc.fixed_fluid_T;
        return { bc.fixed_h, T_f };
    }

    // CONV_WIND
    double T_fluid = w.air_temp;
    double h_total = bc.wind_coeff_A + bc.wind_coeff_B * w.wind;
    return { h_total, T_fluid };
}

// 计算几何特征码 (简单的哈希组合)
size_t ThermalSolver::compute_geometry_checksum(int samples) {
    size_t seed = 0;

    // 1. 混入节点数量
    seed ^= std::hash<size_t>{}(nodes.size()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    // 2. 混入采样参数 (如果采样数变了，缓存必须失效)
    seed ^= std::hash<int>{}(samples)+0x9e3779b9 + (seed << 6) + (seed >> 2);

    // 3. 混入几何特征 (采样部分节点的质心)
    // 为了速度，只采样 10% 的节点，或者步长为 10
    size_t step = std::max((size_t)1, nodes.size() / 100);
    for (size_t i = 0; i < nodes.size(); i += step) {
        const auto& c = nodes[i].centroid;
        // 简单的异或哈希
        auto h1 = std::hash<double>{}(c.x);
        auto h2 = std::hash<double>{}(c.y);
        auto h3 = std::hash<double>{}(c.z);
        seed ^= h1 ^ (h2 << 1) ^ (h3 << 2);
    }
    return seed;
}

void ThermalSolver::save_vf_cache(const std::string& filename, int samples) {
    std::cout << "[Cache] Saving View Factors to " << filename << "..." << std::endl;
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "[Cache] Failed to open file for writing!" << std::endl;
        return;
    }

    // 1. 写入头信息 (魔数 + 版本 + 节点数)
    const char* MAGIC = "VF_BIN";
    int version = 1;
    size_t num_nodes = nodes.size();

    out.write(MAGIC, 6);
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));
    out.write(reinterpret_cast<const char*>(&num_nodes), sizeof(num_nodes));
    out.write(reinterpret_cast<const char*>(&samples), sizeof(samples));

    // 2. 写入每个节点的数据
    for (const auto& node : nodes) {
        // --- Front ---
        out.write(reinterpret_cast<const char*>(&node.vf_sky_front), sizeof(double));
        size_t count_f = node.rad_links_front.size();
        out.write(reinterpret_cast<const char*>(&count_f), sizeof(size_t));

        if (count_f > 0) {
            // 直接写入 vector 的内存块 (Pod Type 才可以这样写，RadLink 是 Pod 结构)
            // 注意：如果 RadLink 包含 std::vector 等动态对象则不能这样写，
            // 但你的 RadLink {int, bool, float/double} 是纯数据，可以直接写。
            out.write(reinterpret_cast<const char*>(node.rad_links_front.data()),
                count_f * sizeof(RadLink));
        }

        // --- Back ---
        out.write(reinterpret_cast<const char*>(&node.vf_sky_back), sizeof(double));
        size_t count_b = node.rad_links_back.size();
        out.write(reinterpret_cast<const char*>(&count_b), sizeof(size_t));

        if (count_b > 0) {
            out.write(reinterpret_cast<const char*>(node.rad_links_back.data()),
                count_b * sizeof(RadLink));
        }
    }

    //  获取文件大小
    // out.tellp() 会返回当前写入指针的位置，也就是文件的总字节数
    long long file_size_bytes = out.tellp();

    out.close();
    std::cout << "[Cache] Save complete (" << file_size_bytes / 1024 / 1024 << " MB)." << std::endl;
}

bool ThermalSolver::load_vf_cache(const std::string& filename, int samples) {
    std::ifstream in(filename, std::ios::binary);

    //  检查文件流状态，如果没打开(即文件不存在)，直接返回 false
    if (!in.is_open()) {
        return false;
    }

    std::cout << "[Cache] Found cache file: " << filename << ", loading..." << std::endl;

    // 1. 校验头信息
    char magic[6];
    in.read(magic, 6);
    if (std::string(magic, 6) != "VF_BIN") {
        std::cout << "[Cache] Invalid file format." << std::endl;
        return false;
    }

    int ver;
    in.read(reinterpret_cast<char*>(&ver), sizeof(ver));

    size_t cached_nodes;
    in.read(reinterpret_cast<char*>(&cached_nodes), sizeof(cached_nodes));
    if (cached_nodes != nodes.size()) {
        std::cout << "[Cache] Node count mismatch (" << cached_nodes << " vs " << nodes.size() << "). Ignoring." << std::endl;
        return false;
    }

    int cached_samples;
    in.read(reinterpret_cast<char*>(&cached_samples), sizeof(cached_samples));
    if (cached_samples != samples) {
        std::cout << "[Cache] Sample count mismatch. Ignoring." << std::endl;
        return false;
    }

    // 2. 读取数据
    for (auto& node : nodes) {
        // --- Front ---
        in.read(reinterpret_cast<char*>(&node.vf_sky_front), sizeof(double));
        size_t count_f;
        in.read(reinterpret_cast<char*>(&count_f), sizeof(size_t));

        node.rad_links_front.resize(count_f);
        if (count_f > 0) {
            in.read(reinterpret_cast<char*>(node.rad_links_front.data()),
                count_f * sizeof(RadLink));
        }

        // --- Back ---
        in.read(reinterpret_cast<char*>(&node.vf_sky_back), sizeof(double));
        size_t count_b;
        in.read(reinterpret_cast<char*>(&count_b), sizeof(size_t));

        node.rad_links_back.resize(count_b);
        if (count_b > 0) {
            in.read(reinterpret_cast<char*>(node.rad_links_back.data()),
                count_b * sizeof(RadLink));
        }
    }

    std::cout << "[Cache] Loaded successfully!" << std::endl;
    return true;
}

void ThermalSolver::solve_radiosity_system(double env_temp_K) {
    const double SIGMA = 5.67e-8;
    double E_env = SIGMA * std::pow(env_temp_K, 4.0);

    // 收敛参数
    const double CONVERGENCE_TOL = 1e-6; // 辐射度收敛阈值 (W/m2)
    const int MAX_RAD_ITERS = 100;       // 最大迭代次数防止死循环

    // ---------------------------------------------------------
    // 1. [优化] 预计算发射项 E (避免在迭代中重复 pow 计算)
    //    同时完成 J 的初始化
    // ---------------------------------------------------------
    // 使用 std::vector 临时存储 E，避免修改 ThermalNode 结构体
    std::vector<double> E_emit_f(nodes.size());
    std::vector<double> E_emit_b(nodes.size());

#pragma omp parallel for schedule(static)
    for (int i = 0; i < nodes.size(); ++i) {
        ThermalNode& node = nodes[i];

        // Front
        double T_f_K = node.T_front_next + 273.15;
        E_emit_f[i] = node.ir_emissivity * SIGMA * std::pow(T_f_K, 4.0);
        node.J_front = E_emit_f[i] + (1.0 - node.ir_emissivity) * node.vf_sky_front * E_env;

        // Back
        double T_b_K = node.T_back_next + 273.15;
        E_emit_b[i] = node.ir_emissivity * SIGMA * std::pow(T_b_K, 4.0);
        node.J_back = E_emit_b[i] + (1.0 - node.ir_emissivity) * node.vf_sky_back * E_env;
    }

    // 2. 迭代求解 (Gauss-Seidel)
    int iter = 0;
    double final_max_diff = 0.0;
    for (; iter < MAX_RAD_ITERS; ++iter) {
        double max_diff = 0.0; // 记录本轮最大变化
        for (int i = 0; i < nodes.size(); ++i) {
            ThermalNode& node = nodes[i];
            double rho = 1.0 - node.ir_emissivity;

            // --- 更新 Front J ---
            double H_inc_front = 0.0;
            for (const auto& link : node.rad_links_front) {
                // 对方如果是正面，取 J_front；如果是背面，取 J_back
                double J_target = link.target_is_front ? nodes[link.target_node_idx].J_front
                    : nodes[link.target_node_idx].J_back;
                H_inc_front += link.view_factor * J_target;
            }
            H_inc_front += node.vf_sky_front * E_env;

            double J_new_f = E_emit_f[i] + rho * H_inc_front;
            double diff_f = std::abs(J_new_f - node.J_front);
            if (diff_f > max_diff) max_diff = diff_f;
            node.J_front = J_new_f;// 更新

            // --- 更新 Back J ---
            if (node.bc_back.type != CONV_INSULATED) {
                double H_inc_back = 0.0;
                for (const auto& link : node.rad_links_back) {
                    double J_target = link.target_is_front ? nodes[link.target_node_idx].J_front
                        : nodes[link.target_node_idx].J_back;
                    H_inc_back += link.view_factor * J_target;
                }
                H_inc_back += node.vf_sky_back * E_env;

                double J_new_b = E_emit_b[i] + rho * H_inc_back;
                double diff_b = std::abs(J_new_b - node.J_back);
                if (diff_b > max_diff) max_diff = diff_b;
                node.J_back = J_new_b;
            }
        }
        // 如果变化很小，提前退出
        final_max_diff = max_diff;
        if (max_diff < CONVERGENCE_TOL) {
            iter++; // 
            break;
        }
    }

    // [建议] 日志输出：未收敛警告
    if (iter >= MAX_RAD_ITERS) {
        std::cerr << "[WARNING] Radiosity solver NOT converged after " << MAX_RAD_ITERS
            << " iters. Max Resid: " << final_max_diff << std::endl;
    }

    // 3. 计算净热流 (使用高精度辐照度法)
    // 公式: Q_net = Area * epsilon * (H_inc - sigma * T^4)
    // 这避免了 epsilon -> 1 时的除零分母问题
#pragma omp parallel for schedule(static)
    for (int i = 0; i < nodes.size(); ++i) {
        ThermalNode& node = nodes[i];
        double eps = node.ir_emissivity;
        // --- Front Heat Flux ---
        // 必须最后重新计算一次 H_inc，因为对于高发射率物体，J 里不包含 H 的信息
        double H_final_f = node.vf_sky_front * E_env;
        for (const auto& link : node.rad_links_front) {
            double J_target = link.target_is_front ? nodes[link.target_node_idx].J_front
                : nodes[link.target_node_idx].J_back;
            H_final_f += link.view_factor * J_target;
        }

        // Front Q_rad
        // Q = Absorbed - Emitted
        node.Q_rad_front = node.area * (eps * H_final_f - E_emit_f[i]);

        if (node.bc_back.type != CONV_INSULATED) {
            double H_final_b = node.vf_sky_back * E_env;
            for (const auto& link : node.rad_links_back) {
                double J_target = link.target_is_front ? nodes[link.target_node_idx].J_front
                    : nodes[link.target_node_idx].J_back;
                H_final_b += link.view_factor * J_target;
            }

            // Back Q_rad
            node.Q_rad_back = node.area * (eps * H_final_b - E_emit_b[i]);
		} else {
			node.Q_rad_back = 0.0;
		}
    }
}

void ThermalSolver::solve_step(double dt, double hour, const Vec3& sun_dir, WeatherSystem& weather, bool is_steady_init) {
    int N = nodes.size();
    int DOFs = 2 * N;
    WeatherData w = weather.get_weather(hour);

    // ==========================================
    // 1. 收敛控制参数 (TAItherm 标准)
    // ==========================================
    double TOL_RESID = 1e-2;       // 温度容差 (Tolerance)
    double TOL_SLOPE = 1e-5;       // 容差斜率 (Tolerance Slope)
    int max_iters = is_steady_init ? 100 : 50; // 初始最大迭代数
    int current_iter = 0;

    // [修正2] 定义松弛因子 (0.0 < omega <= 1.0)
    // 0.6 意味着新值取 60%，旧值保留 40%。这能有效抑制 T^4 的震荡。
    double relaxation = is_steady_init ? 0.5 : 0.8;

    // 状态记录变量
    double max_resid_curr = 0.0;    // 当前步最大温差 (Delta T)
    double max_resid_prev = 0.0;    // 上一步最大温差
    double slope = 0.0;             // 容差斜率
    bool converged = false;

    // 如果是稳态初始化，我们将 dt 设为一个巨大的数，甚至可以忽略 mass 项
    double eff_dt = is_steady_init ? 1e4 : dt;

    // 1. 预计算环境参数
    bool use_lwir = (w.lwir > 10.0);
    // 计算天空等效辐射温度 (Kelvin)
    double T_sky_K = 0.0;
    if (use_lwir) {
        // 如果有长波辐射实测值: T_sky = (LWIR / sigma)^0.25
        T_sky_K = std::pow(w.lwir / SIGMA, 0.25);
    }
    else {
        // 简易估算: 空气温度 * 0.9 (转为K)
        T_sky_K = (w.air_temp + 273.15) * 0.9;
    }
   
    

    // 初始化待求温度 (初猜值)
    for (int i = 0; i < N; ++i) {
        // 使用当前时刻的值作为非线性迭代的起点
        nodes[i].T_front_next = nodes[i].T_front;
        nodes[i].T_back_next = nodes[i].T_back;
    }

    // ==========================================
    // 2. 非线性迭代循环 (Outer Loop)
    // ==========================================
    while (current_iter < max_iters) {
        // =============================================================
    // [STEP 0] 求解长波辐射网络 (Radiosity)
    // 这一步计算了包含多重反射和天空辐射的净热流 Q_rad_front/back
    // =============================================================
        solve_radiosity_system(T_sky_K);
        // --- A. 组装线性系统 Ax = b ---
        // 注意：必须在循环内组装，因为 h_rad (辐射线性化系数) 随温度变化
        MatrixBuilder mb(DOFs);
        Vector b(DOFs, 0.0);
        Vector x(DOFs, 0.0);

        // 归一化太阳向量，避免重复计算
        Vec3 sun_vec = normalize(sun_dir);

        for (int i = 0; i < N; ++i) {
            ThermalNode& node = nodes[i];
            int idx_F = 2 * i; // Front 自由度索引
            int idx_B = 2 * i + 1; // Back 自由度索引

            // 初猜值：使用上一轮迭代算出来的 _next 温度
            x[idx_F] = node.T_front_next;
            x[idx_B] = node.T_back_next;

            if (node.group_type == TYPE_ASSIGNED) {
                // 对角线设为1，右端项设为目标温度 -> T = T_target
                mb.add(idx_F, idx_F, 1.0); b[idx_F] = node.T_front;
                mb.add(idx_B, idx_B, 1.0); b[idx_B] = node.T_back;
                continue;
            }

            // --- 构建 Front 方程 ---
            // (C/dt + SumG + h_rad) * T_new - Sum(K*T_neigh) = Q_old + Q_src + h_rad*T_old

            // C: 热惯性系数 (Mass * Cp / dt)
            double C = node.mass_node / eff_dt;
            double diag_F = C;// 对角线基础值 (惯性)
            double rhs_F = C * node.T_front;// 右端项基础值 (历史温度惯性)

            // 辐射线性化
            double Tk_f = node.T_front_next + 273.15;
            double h_rad_f = 4.0 * node.ir_emissivity * SIGMA * node.area * std::pow(Tk_f, 3.0);
            diag_F += h_rad_f;
            rhs_F += node.Q_rad_front + h_rad_f * node.T_front_next;

            // 对流
            auto conv_f = get_convection_params(node, true, w);
            diag_F += conv_f.first * node.area;// hA 加到左边
            rhs_F += conv_f.first * node.area * conv_f.second;// hA * T_fluid 加到右边

            // 内部热源
            rhs_F += 0.5 * node.Q_gen_total;

            // 导热 (法向)
            double K = node.conductance;
            diag_F += K;
            mb.add(idx_F, idx_B, -K);

            // 横向导热 (Front) - 承担 50%
            for (auto& link : node.neighbors) {
                double k_lat_part = link.conductance * 0.5; // 分摊 50%
                diag_F += k_lat_part;
                mb.add(idx_F, 2 * link.neighbor_idx, -k_lat_part);
            }

            // ---------------------------------------------------------
            // 构建 Back 方程 (idx_B)
            // 方程: (C + SumG) * T_B_new - Sum(K_neigh * T_neigh_new) = Q_sources
            // ---------------------------------------------------------
            double diag_B = C;              // 对角线基础值 (惯性)
            double rhs_B = C * node.T_back; // 右端项基础值 (历史温度惯性)

            // 法向导热 (连接 Front)
            diag_B += K; // 增加自身稳定性
            mb.add(idx_B, idx_F, -K); // 耦合项: -K * T_Front

            // 对流 (Back)
            auto conv_b = get_convection_params(node, false, w);
            diag_B += conv_b.first * node.area;
            rhs_B += conv_b.first * node.area * conv_b.second;

            // 内部热源 (Back Share 50%)
            rhs_B += 0.5 * node.Q_gen_total;

            // ---------------------------------------------------------
            // ★★★ 核心修改：双面太阳光照计算 ★★★
            // ---------------------------------------------------------
            double cos_theta = dot(node.normal, sun_vec);

            // 计算总吸收功率 (注意这里用 abs，不再强制为 0)
            // 公式: Solar * Area * Absorp * Shadow * |cos(theta)|
            double solar_power = w.solar * node.area * node.solar_absorp * node.shadow_factor * std::abs(cos_theta);

            node.Q_solar_absorbed = solar_power; // 记录总吸收用于调试

            if (cos_theta > 0.0) {
                // 太阳在正面 (Front Facing) -> 加热 Front 节点
                rhs_F += solar_power;
            }
            else {
                // 太阳在背面 (Back Facing) -> 加热 Back 节点
                rhs_B += solar_power;
            }

            // 辐射 (Back) - 仅当非绝热时计算
            if (node.bc_back.type != CONV_INSULATED) {
                // 基础 MCRT 辐射热流
                rhs_B += node.Q_rad_back;

                // 辐射线性化 (Back)
                double Tk_b = node.T_back_next + 273.15;
                double h_rad_b = 4.0 * node.ir_emissivity * SIGMA * node.area * std::pow(Tk_b, 3.0);

                diag_B += h_rad_b;               // 加到左边
                rhs_B += h_rad_b * node.T_back_next;  // 加到右边
            }

            // 横向导热 (Back) - 承担剩余 50% 
            // 允许热量在背面节点之间横向流动
            for (auto& link : node.neighbors) {
                double k_lat_part = link.conductance * 0.5; // 分摊 50%
                diag_B += k_lat_part;

                // 注意：邻居的 Back 节点索引是 [2 * neighbor_idx + 1]
                mb.add(idx_B, 2 * link.neighbor_idx + 1, -k_lat_part);
            }

            // 写入 Front 矩阵行
            mb.add(idx_F, idx_F, diag_F);
            b[idx_F] = rhs_F;

            // 写入 Back 矩阵行
            mb.add(idx_B, idx_B, diag_B);
            b[idx_B] = rhs_B;
        }

        // --- B. AMG 求解 ---
        // 2. 组装矩阵
        SparseMatrix A = mb.build();

        // 3. AMG Setup (如果是第一次，或者拓扑改变了，或者是线性化导致矩阵变了)
        // 对于非线性热问题，矩阵A随温度变化 (h_rad)，所以每步都要 Setup (或者每几步)
        amg_solver.setup(A);

        // 4. AMG Solve
        // 初始猜测 x 已经是 T_front_next
        amg_solver.solve(b, x);

        // --- C. 计算残差和斜率 ---
        max_resid_curr = 0.0;
        for (int i = 0; i < N; ++i) {
            // Front
            double val_F = x[2 * i];
            // [安全钳位]
            if (std::isnan(val_F)) val_F = nodes[i].T_front_next;
            if (val_F < -270.0) val_F = -270.0;
            if (val_F > 5000.0) val_F = 5000.0;

            // [应用松弛] New = Old + omega * (Calc - Old)
            double T_F_new = nodes[i].T_front_next + relaxation * (val_F - nodes[i].T_front_next);

            // Back
            double val_B = x[2 * i + 1];
            if (std::isnan(val_B)) val_B = nodes[i].T_back_next;
            if (val_B < -270.0) val_B = -270.0;
            if (val_B > 5000.0) val_B = 5000.0;

            double T_B_new = nodes[i].T_back_next + relaxation * (val_B - nodes[i].T_back_next);
            double diff_F = std::abs(T_F_new - nodes[i].T_front_next);
            double diff_B = std::abs(T_B_new - nodes[i].T_back_next);
            max_resid_curr = std::max(max_resid_curr, std::max(diff_F, diff_B));

            // 更新温度用于下一次迭代
            nodes[i].T_front_next = T_F_new;
            nodes[i].T_back_next = T_B_new;
        }

        if (current_iter > 0) {
            slope = std::abs(max_resid_curr - max_resid_prev);
        }
        max_resid_prev = max_resid_curr;
        current_iter++;

        // 调试打印 (每 10 步看一次)
        if (is_steady_init && current_iter % 10 == 0) {
            std::cout << " Iter " << current_iter
                << " Resid: " << max_resid_curr
                << " Slope: " << slope << std::endl;
        }

        // --- D. 检查收敛 ---
        // 逻辑：只有当 (残差 < 阈值) 且 (斜率 < 阈值) 时才算收敛
        if (current_iter > 2 && max_resid_curr < TOL_RESID && slope < TOL_SLOPE) {
            converged = true;
            break;
        }

        // ==========================================
        // 3. 交互式收敛警告 (Convergence Warning)
        // ==========================================
        if (current_iter >= max_iters && !converged) {
            std::cout << "\n\n[WARNING] Convergence Warning at Hour " << hour << std::endl;
            std::cout << "  Iterations reached maximum: " << max_iters << std::endl;
            std::cout << "  Current Residual: " << max_resid_curr << " (Target: " << TOL_RESID << ")" << std::endl;
            std::cout << "  Current Slope:    " << slope << " (Target: " << TOL_SLOPE << ")" << std::endl;
            std::cout << "------------------------------------------------" << std::endl;
            std::cout << "Select action:\n";
            std::cout << "  [1] Increase Iterations (+50) and Continue\n";
            std::cout << "  [2] Accept Current Solution and Proceed\n";
            std::cout << "Input choice (1 or 2): ";

            int choice;
            if (std::cin >> choice) {
                if (choice == 1) {
                    max_iters += 50;
                    std::cout << ">>> Extending iterations to " << max_iters << "..." << std::endl;
                }
                else {
                    std::cout << ">>> Accepting current solution." << std::endl;
                    converged = true; // 强制设为收敛以跳出
                    break;
                }
            }
            else {
                // 防止输入错误死循环
                std::cin.clear();
                std::cin.ignore(10000, '\n');
                std::cout << "Invalid input. Accepting solution." << std::endl;
                break;
            }
        }
    }

    // 5. 更新回物理节点
    for (int i = 0; i < N; ++i) {
        // 更新显示值
        nodes[i].T_front = nodes[i].T_front_next;
        nodes[i].T_back = nodes[i].T_back_next;
    }
}

