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
    std::cout << "[MCRT] Calculating View Factors (" << samples << " rays/node)..." << std::endl;

    #pragma omp parallel for 
    for (int i = 0; i < nodes.size(); ++i) {
        ThermalNode& node = nodes[i];

        // 计算节点的【正面】辐射
        // 临时统计 map: key -> pair<目标ID, 目标正反>, value -> 命中次数
        std::map<std::pair<int, bool>, int> hits_front;
        int sky_hits_front = 0;
        Vec3 origin_f = node.centroid + node.normal * 0.01;// 发射点：质心沿法线向外偏移一点点

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
            Vec3 origin_b = node.centroid - node.normal * 0.01;// 向内偏移
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
}

void ThermalSolver::update_shadows(const Vec3& sun_dir) {
    Vec3 ray_dir = normalize(sun_dir); 
    const double BIAS = 0.05; // 偏移量，防止自己遮挡自己
    const double MAX_DIST = 1.0e20; // 太阳视为无限远
    
    // 开启 OpenMP 并行计算，BVH 是只读的，所以并行是安全的
    #pragma omp parallel for 
    for (int i = 0; i < nodes.size(); ++i) {
        ThermalNode& receiver = nodes[i];
        Vec3 face_normal = receiver.normal;
        if (dot(receiver.normal, ray_dir) < 0.0) face_normal = receiver.normal * -1.0;
        Vec3 origin = receiver.centroid + (face_normal * BIAS);
        
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

void ThermalSolver::solve_radiosity_system(double env_temp_K) {
    const double SIGMA = 5.67e-8;
    double E_env = SIGMA * std::pow(env_temp_K, 4.0);

    // 1. 初始化 J (基于当前温度)
    for (auto& node : nodes) {
        double T_f_K = node.T_front + 273.15;
        double T_b_K = node.T_back + 273.15;

        double E_f = node.ir_emissivity * SIGMA * std::pow(T_f_K, 4.0);
        double E_b = node.ir_emissivity * SIGMA * std::pow(T_b_K, 4.0);

        // 初始猜测：只有自身发射
        node.J_front = E_f + (1.0 - node.ir_emissivity) * E_env * node.vf_sky_front;
        node.J_back = E_b + (1.0 - node.ir_emissivity) * E_env * node.vf_sky_back;
    }

    // 2. 迭代求解 (Gauss-Seidel)
    int iterations = 10;
    for (int k = 0; k < iterations; ++k) {
        for (auto& node : nodes) {
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

            double T_f_K = node.T_front + 273.15;
            double E_f = node.ir_emissivity * SIGMA * std::pow(T_f_K, 4.0);
            node.J_front = E_f + rho * H_inc_front;

            // --- 更新 Back J ---
            double H_inc_back = 0.0;
            for (const auto& link : node.rad_links_back) {
                double J_target = link.target_is_front ? nodes[link.target_node_idx].J_front
                    : nodes[link.target_node_idx].J_back;
                H_inc_back += link.view_factor * J_target;
            }
            H_inc_back += node.vf_sky_back * E_env;

            double T_b_K = node.T_back + 273.15;
            double E_b = node.ir_emissivity * SIGMA * std::pow(T_b_K, 4.0);
            node.J_back = E_b + rho * H_inc_back;
        }
    }

    // 3. 计算净热流 Q_rad (Watts)
    for (auto& node : nodes) {
        double eps = node.ir_emissivity;
        if (eps > 0.999) eps = 0.999;

        // Front Q_rad
        double T_f_K = node.T_front + 273.15;
        double E_f = SIGMA * std::pow(T_f_K, 4.0);
        // Q_net_gain = - Q_net_loss
        node.Q_rad_front = -(node.area * eps * (E_f - node.J_front) / (1.0 - eps));

        // Back Q_rad
        double T_b_K = node.T_back + 273.15;
        double E_b = SIGMA * std::pow(T_b_K, 4.0);
        node.Q_rad_back = -(node.area * eps * (E_b - node.J_back) / (1.0 - eps));
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

            // 太阳辐射 
            double n_dot_s = dot(node.normal, normalize(sun_dir));
            rhs_F += w.solar * node.area * node.solar_absorp * node.shadow_factor * std::max(0.0, n_dot_s);

            // 内部热源
            rhs_F += 0.5 * node.Q_gen_total;

            // 导热 (法向)
            double K = node.conductance;
            diag_F += K;
            mb.add(idx_F, idx_B, -K);

            // 导热 (横向)
            for (auto& link : node.neighbors) {
                diag_F += link.conductance;
                mb.add(idx_F, 2 * link.neighbor_idx, -link.conductance);
            }

            // 写入 Front 矩阵行
            mb.add(idx_F, idx_F, diag_F);
            b[idx_F] = rhs_F;

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

            // 横向导热 (Back)
            // 目前模型假设横向导热主要通过 Front 节点连接，Back 节点通常无横向连接
            // 如果您的模型支持多层横向导热，可以在此添加循环

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

