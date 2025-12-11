#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <algorithm>
#include <iomanip>
#include <regex>
#include <omp.h>

// ==========================================
// 1. 基础数学工具
// ==========================================
const double EPSILON = 1e-6;
const double SIGMA = 5.67e-8;

struct Vec3 {
    double x, y, z;
    Vec3 operator+(const Vec3& v) const { return { x + v.x, y + v.y, z + v.z }; }
    Vec3 operator-(const Vec3& v) const { return { x - v.x, y - v.y, z - v.z }; }
    Vec3 operator*(double s) const { return { x * s, y * s, z * s }; }
    Vec3 operator/(double s) const { return { x / s, y / s, z / s }; }
};
double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
Vec3 normalize(const Vec3& v) { double len = std::sqrt(dot(v, v)); return len > EPSILON ? v / len : v; }
Vec3 cross(const Vec3& a, const Vec3& b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }

// ==========================================
// 2. 天气系统 (重型清洗版 - 解决格式错乱)
// ==========================================
struct WeatherData {
    double time_hour;
    double air_temp;
    double solar;
    double wind;
    double lwir;
};

class WeatherSystem {
public:
    std::vector<WeatherData> data_points;

    bool load_weather(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cout << "[Error] Cannot open weather file: " << filename << std::endl;
            return false;
        }

        // 1. 读取整个文件到字符串
        std::stringstream buffer;
        buffer << file.rdbuf();
        std::string content = buffer.str();

        std::cout << "[Weather] File read. Size: " << content.size() << " chars." << std::endl;

        // 2. 找到数据起始点 (TIME 之后)
        // 这样可以跳过文件开头的 "04" 这种无关数字
        size_t data_start = content.find("TIME");
        if (data_start != std::string::npos) {
            content = content.substr(data_start);
        }

        // 3. 正则清洗 (这是最关键的一步)
        try {
            // A. 去除 这种标签
            std::regex re_source("\\[.*?\\]");
            content = std::regex_replace(content, re_source, " ");

            // B. 去除所有字母 (TIME, AIRT 等标题)
            std::regex re_letters("[A-Za-z]");
            content = std::regex_replace(content, re_letters, " ");

            // C. 将冒号等标点替换为空格 (防止 09:20 读错)
            std::regex re_punct("[:]");
            content = std::regex_replace(content, re_punct, " ");
        }
        catch (const std::regex_error& e) {
            std::cout << "[Error] Regex failed: " << e.what() << std::endl;
            return false;
        }

        // 4. 流式读取纯数字
        std::stringstream ss(content);
        std::vector<double> all_numbers;
        double val;
        while (ss >> val) {
            all_numbers.push_back(val);
        }

        std::cout << "[Weather] Extracted " << all_numbers.size() << " numbers." << std::endl;

        // 5. 按列重组 (你的文件是 8 列数据)
        // 0:TIME, 1:AIRT, 2:SOLAR, 3:WIND, 4:HUMID, 5:CLOUD, 6:LWIR, 7:WINDIR(or Rain)
        int stride = 8;

        if (all_numbers.size() < stride) {
            std::cout << "[Error] Not enough data found!" << std::endl;
            return false;
        }

        double last_raw_time = -1.0;
        double day_offset = 0.0;

        for (size_t i = 0; i + stride <= all_numbers.size(); i += stride) {
            double raw_time = all_numbers[i];     // HHMM
            double airt = all_numbers[i + 1];
            double solar = all_numbers[i + 2];
            double wind = all_numbers[i + 3];
            double lwir = all_numbers[i + 6];
            // 格式解析 HHMM -> Hour
            int hh = (int)(raw_time / 100);
            int mm = (int)(raw_time) % 100;

            // 基本校验
            if (hh > 24 || mm >= 60) continue; // 跳过非法时间

            // 处理跨天 (2355 -> 0000)
            if (last_raw_time >= 0 && raw_time < last_raw_time) {
                // 如果时间突然变小 (且不是一点点波动)，说明过了一天
                // 比如 2355 -> 0000
                day_offset += 24.0;
            }
            last_raw_time = raw_time;

            double final_hour = hh + mm / 60.0 + day_offset;

            data_points.push_back({ final_hour, airt, solar, wind, lwir });
        }

        std::cout << "[Weather] Successfully parsed " << data_points.size() << " time steps." << std::endl;
        if (!data_points.empty()) {
            std::cout << "          Range: " << data_points.front().time_hour << "h -> "
                << data_points.back().time_hour << "h" << std::endl;
        }

        return true;
    }

    WeatherData get_weather(double query_hour) {
        if (data_points.empty()) return { query_hour, 20.0, 0.0, 1.0 };

        if (query_hour <= data_points.front().time_hour) return data_points.front();
        if (query_hour >= data_points.back().time_hour) return data_points.back();

        // 线性查找插值
        for (size_t i = 0; i < data_points.size() - 1; ++i) {
            if (query_hour >= data_points[i].time_hour && query_hour < data_points[i + 1].time_hour) {
                const auto& p1 = data_points[i];
                const auto& p2 = data_points[i + 1];

                double ratio = (query_hour - p1.time_hour) / (p2.time_hour - p1.time_hour);

                WeatherData res;
                res.time_hour = query_hour;
                res.air_temp = p1.air_temp + ratio * (p2.air_temp - p1.air_temp);
                res.solar = p1.solar + ratio * (p2.solar - p1.solar);
                res.wind = p1.wind + ratio * (p2.wind - p1.wind);
                res.lwir = p1.lwir + ratio * (p2.lwir - p1.lwir);
                return res;
            }
        }
        return data_points.back();
    }
};

WeatherSystem g_weather;

// ==========================================
// 3. 核心数据结构升级 (支持 Assigned 和 Insulated)
// ==========================================

// 3.1 材质定义
struct Material {
    double k;           // 导热率
    double rho;         // 密度
    double Cp;          // 比热
    double absorptivity; // 太阳吸收率
    double emissivity;   // 红外发射率
};

// 3.2 对流/边界类型枚举
enum ConvectionType {
    CONV_WIND,          // 随天气风速变化
    CONV_FIXED_H_T,     // 固定 h 和 T
    CONV_INSULATED      // ★★★ 新增：绝热 (Q=0) ★★★
};

struct ConvectionBC {
    ConvectionType type;
    double fixed_h;
    double fixed_fluid_T;
    double wind_coeff_A;
    double wind_coeff_B;
};

// 3.3 Group 类型枚举
enum GroupType {
    TYPE_CALCULATED,    // 计算温度 (正常物理对象)
    TYPE_ASSIGNED       // ★★★ 新增：指定温度 (恒温源) ★★★
};

// 3.4 部件属性配置
struct PartProperty {
    std::string material_name;
    double thickness;
    double initial_temp; // 如果是 Assigned 类型，这个值就是全过程的固定温度

    GroupType group_type; // Assigned or Calculated

    ConvectionBC front_bc;
    ConvectionBC back_bc;
};

// 全局配置库
std::map<std::string, Material> MAT_LIB;
std::map<std::string, PartProperty> PROJECT_CONFIG;
// ==========================================
// 3.5 配置加载系统 (新增：对接 Python)
// ==========================================
struct GlobalConfig {
    std::string obj_file = "chuan.tai";
    std::string weather_file = "weather.txt";
    double start_time = 6.0;
    double end_time = 16.0;
    double dt = 0.5;
};

class ConfigSystem {
public:
    GlobalConfig global;

    // 辅助：解析 BC 参数行，例如 "FIXED 5.0 20.0"
    ConvectionBC parse_bc(std::stringstream& ss) {
        std::string type_str;
        ss >> type_str;

        ConvectionBC bc;
        // 默认值
        bc.type = CONV_WIND;
        bc.wind_coeff_A = 5.7; bc.wind_coeff_B = 3.8;
        bc.fixed_h = 5.0; bc.fixed_fluid_T = 0.0;

        if (type_str == "WIND") {
            bc.type = CONV_WIND;
            // 可选：读取自定义风速系数 "WIND 5.7 3.8"
            if (!ss.eof()) ss >> bc.wind_coeff_A >> bc.wind_coeff_B;
        }
        else if (type_str == "FIXED") {
            bc.type = CONV_FIXED_H_T;
            ss >> bc.fixed_h >> bc.fixed_fluid_T;
        }
        else if (type_str == "INSULATED") {
            bc.type = CONV_INSULATED;
        }
        return bc;
    }

    bool load_config(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cout << "[Config] Warning: " << filename << " not found. Using defaults." << std::endl;
            return false; // 返回 false 表示需要使用硬编码默认值
        }

        std::cout << "[Config] Reading " << filename << "..." << std::endl;

        std::string line;
        PartProperty current_prop; // 临时存当前 Group 属性
        std::string current_name = "";
        bool inside_group = false;

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue; // 跳过注释空行
            std::stringstream ss(line);
            std::string key;
            ss >> key;

            // --- 全局参数 ---
            if (key == "OBJ_FILE") ss >> global.obj_file;
            else if (key == "WEATHER_FILE") ss >> global.weather_file;
            else if (key == "TIME_PARAMS") ss >> global.start_time >> global.end_time >> global.dt;

            // --- Group 解析 ---
            else if (key == "BEGIN_GROUP") {
                inside_group = true;
                // 重置临时变量为默认值
                current_prop = { "Steel", 0.01, 20.0, TYPE_CALCULATED,
                                 {CONV_WIND, 0,0,10,3}, {CONV_FIXED_H_T, 5,0,0,0} };
                current_name = "Default";
            }
            else if (key == "END_GROUP") {
                if (inside_group && !current_name.empty()) {
                    PROJECT_CONFIG[current_name] = current_prop;
                    // std::cout << "  Configured Group: " << current_name << std::endl;
                }
                inside_group = false;
            }
            else if (inside_group) {
                if (key == "NAME") {
                    // 读取剩余整行作为名字（防止名字里有空格）
                    std::getline(ss, current_name);
                    // 去掉开头的空格
                    size_t first = current_name.find_first_not_of(" \t");
                    if (first != std::string::npos) current_name = current_name.substr(first);
                }
                else if (key == "TYPE") {
                    std::string val; ss >> val;
                    current_prop.group_type = (val == "ASSIGNED") ? TYPE_ASSIGNED : TYPE_CALCULATED;
                }
                else if (key == "MAT") ss >> current_prop.material_name;
                else if (key == "THICK") ss >> current_prop.thickness;
                else if (key == "INIT_TEMP") ss >> current_prop.initial_temp;
                else if (key == "BC_FRONT") current_prop.front_bc = parse_bc(ss);
                else if (key == "BC_BACK")  current_prop.back_bc = parse_bc(ss);
            }
        }
        return true;
    }
};
ConfigSystem g_config;
// ==========================================
// 4. 项目配置中心 (Setup Project)
// ==========================================
void setup_project_properties() {
    // 1. 初始化材质库 (这部分通常是固定的，或者也可以从文件读)
    MAT_LIB["Steel"] = { 45.0, 7850.0, 460.0, 0.7, 0.9 };
    MAT_LIB["Iron"] = { 40.0, 7800.0, 500.0, 0.8, 0.85 };
    MAT_LIB["Glass"] = { 1.0,  2500.0, 840.0, 0.1, 0.9 };
    MAT_LIB["Insulation"] = { 0.04, 50.0, 1000.0, 0.5, 0.9 };

    // 2. 尝试加载 Python 生成的 config.txt
    if (!g_config.load_config("config.txt")) {
        // 如果没找到文件，可以在这里保留一些硬编码的默认值用于测试
        std::cout << "[Config] Using internal hardcoded defaults." << std::endl;

        // ... (把你原来的 PROJECT_CONFIG["Hull"] = ... 代码放在这里作为保底)
        ConvectionBC bc_weather = { CONV_WIND, 0, 0, 10.0, 3.0 };
        ConvectionBC bc_cabin = { CONV_FIXED_H_T, 5.0, 0.0, 0, 0 };
        PROJECT_CONFIG["Default"] = { "Steel", 0.01, 20.0, TYPE_CALCULATED, bc_weather, bc_cabin };
    }
}

PartProperty get_part_property(const std::string& group_name) {
    if (PROJECT_CONFIG.count(group_name)) return PROJECT_CONFIG[group_name];
    for (const auto& kv : PROJECT_CONFIG) {
        if (group_name.find(kv.first) != std::string::npos && kv.first != "Default") return kv.second;
    }
    return PROJECT_CONFIG["Default"];
}

// ==========================================
// 5. 计算节点
// ==========================================
struct Triangle { Vec3 v0, v1, v2; };

// 唯一的顶点 ID (用于焊接)
struct VertexKey {
    long long x, y, z;
    bool operator<(const VertexKey& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return z < other.z;
    }
};

// 简单的边结构 (由两个顶点ID组成)
struct EdgeKey {
    int v1, v2;
    bool operator<(const EdgeKey& other) const {
        if (v1 != other.v1) return v1 < other.v1;
        return v2 < other.v2;
    }
};

// 将坐标量化为整数 Key (容差 1mm)
VertexKey to_key(const Vec3& v) {
    const double SCALE = 1000.0; // 精度 1mm
    return { (long long)std::round(v.x * SCALE), (long long)std::round(v.y * SCALE), (long long)std::round(v.z * SCALE) };
}

struct EdgeLink {
    int neighbor_idx;
    double conductance;
};

struct ThermalNode {
    double T_front, T_back;
    double T_front_next, T_back_next;

    double mass_node, conductance, area;
    double solar_absorp, ir_emissivity;

    ConvectionBC bc_front, bc_back;
    GroupType group_type; // 记录类型：ASSIGNED 或 CALCULATED

    Vec3 centroid, normal;
    double shadow_factor;
    std::string part_name;
    std::vector<Triangle> geometry_tris;

    // 横向邻居列表 
    std::vector<EdgeLink> neighbors;
    double k_mat; // 需要存一下材料的 k 值，用于计算横向导热
    double thickness; // 存一下厚度
};
std::vector<ThermalNode> g_nodes;

// ==========================================
// 6. 几何与加载
// ==========================================
bool ray_intersects_triangle(const Vec3& ray_origin, const Vec3& ray_dir, const Triangle& tri) {
    Vec3 edge1 = tri.v1 - tri.v0; Vec3 edge2 = tri.v2 - tri.v0; Vec3 h = cross(ray_dir, edge2);
    double a = dot(edge1, h); if (a > -EPSILON && a < EPSILON) return false;
    double f = 1.0 / a; Vec3 s = ray_origin - tri.v0; double u = f * dot(s, h);
    if (u < 0.0 || u > 1.0) return false; Vec3 q = cross(s, edge1); double v = f * dot(ray_dir, q);
    if (v < 0.0 || u + v > 1.0) return false; double t = f * dot(edge2, q); return t > EPSILON;
}

void update_shadows(const Vec3& sun_dir) {
    Vec3 ray_dir = normalize(sun_dir); const double BIAS = 0.05;
#pragma omp parallel for 
    for (int i = 0; i < g_nodes.size(); ++i) {
        ThermalNode& receiver = g_nodes[i];
        Vec3 face_normal = receiver.normal;
        if (dot(receiver.normal, ray_dir) < 0.0) face_normal = receiver.normal * -1.0;
        Vec3 origin = receiver.centroid + (face_normal * BIAS);
        bool blocked = false;
        for (int j = 0; j < g_nodes.size(); ++j) {
            if (i == j) continue;
            for (const auto& tri : g_nodes[j].geometry_tris) {
                if (ray_intersects_triangle(origin, ray_dir, tri)) { blocked = true; break; }
            } if (blocked) break;
        } receiver.shadow_factor = blocked ? 0.0 : 1.0;
    }
}

bool init_model(const std::string& filename) {
    setup_project_properties(); // 加载配置
    std::ifstream file(filename); if (!file.is_open()) return false;
    std::vector<Vec3> temp_verts; std::string line, group = "Default";
    std::cout << "Loading model..." << std::endl;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; std::stringstream ss(line); std::string type; ss >> type;
        if (type == "v") { double x, y, z; ss >> x >> y >> z; temp_verts.push_back({ x, y, z }); }
        else if (type == "g") { std::string temp; std::getline(ss, temp); size_t first = temp.find_first_not_of(' '); group = (first != std::string::npos) ? temp.substr(first) : "Default"; }
        else if (type == "f") {
            std::vector<int> idxs; int idx; while (ss >> idx) idxs.push_back(idx - 1); if (idxs.size() < 3) continue;
            std::vector<Triangle> tris; Vec3 v0 = temp_verts[idxs[0]], v1 = temp_verts[idxs[1]], v2 = temp_verts[idxs[2]];
            tris.push_back({ v0, v1, v2 }); double area = 0.5 * std::sqrt(dot(cross(v1 - v0, v2 - v0), cross(v1 - v0, v2 - v0)));
            if (idxs.size() == 4) { Vec3 v3 = temp_verts[idxs[3]]; tris.push_back({ v0, v2, v3 }); area += 0.5 * std::sqrt(dot(cross(v2 - v0, v3 - v0), cross(v2 - v0, v3 - v0))); }
            if (area < 1e-6) continue;

            PartProperty prop = get_part_property(group);
            Material mat = MAT_LIB[prop.material_name];

            ThermalNode node;
            node.part_name = group;
            node.area = area;
            node.geometry_tris = tris;
            Vec3 center_sum = { 0,0,0 }; for (auto& t : tris) center_sum = center_sum + t.v0 + t.v1 + t.v2;
            node.centroid = center_sum / (double)(tris.size() * 3);
            node.normal = normalize(cross(tris[0].v1 - tris[0].v0, tris[0].v2 - tris[0].v0));

            double total_mass = area * prop.thickness * mat.rho;
            node.mass_node = (total_mass * mat.Cp) / 2.0;
            node.conductance = (mat.k * area) / prop.thickness;
            node.solar_absorp = mat.absorptivity;
            node.ir_emissivity = mat.emissivity;
            node.bc_front = prop.front_bc;
            node.bc_back = prop.back_bc;
            node.group_type = prop.group_type; // 记录类型
            node.k_mat = mat.k;
            node.thickness = prop.thickness;

            // 初始温度：如果是 Assigned，这个值将作为永久固定值
            node.T_front = node.T_back = prop.initial_temp;
            node.shadow_factor = 1.0;
            g_nodes.push_back(node);
        }
    }
    std::cout << "Nodes: " << g_nodes.size() << std::endl;
    return true;
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

// ==========================================
// 全局辅助：存储焊接后的唯一顶点
// ==========================================
std::vector<Vec3> g_unique_verts; // 新增：存坐标

void build_topology() {
    std::cout << "[Topology] Building lateral connections (Advanced Quad Support)..." << std::endl;

    // 全局辅助清理
    g_unique_verts.clear();
    for (auto& n : g_nodes) n.neighbors.clear();

    // 1. 顶点焊接 (保持不变)
    std::map<VertexKey, int> vert_map;
    auto get_vert_id = [&](const Vec3& v) -> int {
        VertexKey key = to_key(v);
        if (vert_map.find(key) == vert_map.end()) {
            int new_id = (int)g_unique_verts.size();
            vert_map[key] = new_id;
            g_unique_verts.push_back(v);
            return new_id;
        }
        return vert_map[key];
    };

    // 2. 全局边注册 (Edge Registry)
    // 记录每一条边被哪些 Node 引用过
    // Map: EdgeKey -> List of Node Indices
    std::map<EdgeKey, std::vector<int>> edge_registry;

    for (int i = 0; i < g_nodes.size(); ++i) {
        // ★★★ 修正点1：遍历该 Node 下所有的三角形 ★★★
        for (const auto& tri : g_nodes[i].geometry_tris) {
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

        ThermalNode& node1 = g_nodes[n1];
        ThermalNode& node2 = g_nodes[n2];

        // 1. 计算真实的接触边长 (保留你的认可)
        const Vec3& p1 = g_unique_verts[edge_key.v1];
        const Vec3& p2 = g_unique_verts[edge_key.v2];
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
}

// ==========================================
// 辅助：辐射线性化系数计算
// ==========================================
double calc_h_rad(double T_surf_K, double T_env_K, double epsilon) {
    // 防止除零错误
    if (std::abs(T_surf_K - T_env_K) < 0.1) {
        return 4.0 * SIGMA * epsilon * std::pow(T_surf_K, 3);
    }
    return SIGMA * epsilon * (std::pow(T_surf_K, 2) + std::pow(T_env_K, 2)) * (T_surf_K + T_env_K);
}

// ==========================================
// 辅助：高级对流系数 (TAITherm 风格)
// ==========================================
// 返回 pair: {h_value, T_fluid}
std::pair<double, double> get_convection_params(const ThermalNode& node, bool is_front, const WeatherData& w) {
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

// ==========================================
// 核心：隐式求解步 (Gauss-Seidel 迭代)
// ==========================================
void solve_step(double dt, double hour, const Vec3& sun_dir, bool is_steady_init = false) {
    WeatherData w = g_weather.get_weather(hour);

    // 如果是稳态初始化，我们将 dt 设为一个巨大的数，甚至可以忽略 mass 项
    double eff_dt = is_steady_init ? 1e6 : dt;

    // 1. 预计算环境参数
    bool use_lwir = (w.lwir > 10.0);
    double T_sky_K = use_lwir ? 0.0 : (w.air_temp + 273.15) * 0.9;
    // 注意：如果有实测 LWIR，sky temperature 并不直接用。这里简化处理。

    // 2. 初始化 Guess 温度 (用上一时刻的温度作为初猜值)
    for (auto& node : g_nodes) {
        node.T_front_next = node.T_front;
        node.T_back_next = node.T_back;
    }

    // 3. 高斯-赛德尔迭代 (Gauss-Seidel Loop)
    // 增加迭代次数可以提高精度，一般 5-10 次足够
    int iterations = is_steady_init ? 50 : 5;

    for (int iter = 0; iter < iterations; ++iter) {

        // 并行计算时，Gauss-Seidel 严格来说需要红黑排序，但对于热传导，
        // 直接并行通常也能收敛，或者使用 OpenMP 的原子操作(虽然慢)。
        // 为了简单展示逻辑，这里暂不开启 omp，确保逻辑清晰。
        // #pragma omp parallel for 
        for (int i = 0; i < g_nodes.size(); ++i) {
            ThermalNode& node = g_nodes[i];

            // --- 处理 Assigned 类型 (恒温) ---
            if (node.group_type == TYPE_ASSIGNED) {
                // 温度保持不变，不需要计算平衡方程
                continue;
            }

            // 获取当前迭代的最新温度估计值
            double Tf_curr = node.T_front_next;
            double Tb_curr = node.T_back_next;

            // ==========================
            // 节点 1: Front Surface
            // ==========================
            double sum_G = 0.0; // 总热导 (Conductance Sum)
            double sum_Q = 0.0; // 总热源 (Source Sum)

            // A. 热容项 (惯性)
            double C_term = node.mass_node / eff_dt;
            sum_G += C_term;
            sum_Q += C_term * node.T_front; // 上一时刻的温度

            // B. 太阳辐射 (作为固定热源 Q)
            double n_dot_s = dot(node.normal, normalize(sun_dir));
            // 只有当法线朝向太阳时才计算 (n_dot_s > 0)
            double cos_theta = std::max(0.0, n_dot_s);
            double Q_sun = w.solar * node.area * node.solar_absorp * node.shadow_factor * cos_theta;
            sum_Q += Q_sun;

            // C. 对流 (线性化: h * (T_fluid - T_surf) -> h*T_fluid - h*T_surf)
            auto conv_f = get_convection_params(node, true, w);
            double G_conv_f = conv_f.first * node.area;
            sum_G += G_conv_f;
            sum_Q += G_conv_f * conv_f.second; // h * T_fluid

            // D. 环境辐射 (线性化)
            double G_rad_f = 0.0;
            double T_rad_env = 0.0;
            if (use_lwir) {
                // Q = eps * A * (LWIR - sigma*T^4) -> 很难完全线性化，
                // 这里采用简化策略：把它当成针对天空温度的辐射
                // 此时 T_sky^4 对应 w.lwir / sigma
                double T_sky_equiv = std::pow(w.lwir / SIGMA, 0.25);
                G_rad_f = calc_h_rad(Tf_curr + 273.15, T_sky_equiv, node.ir_emissivity) * node.area;
                T_rad_env = T_sky_equiv - 273.15;
            }
            else {
                G_rad_f = calc_h_rad(Tf_curr + 273.15, T_sky_K, node.ir_emissivity) * node.area;
                T_rad_env = T_sky_K - 273.15;
            }
            sum_G += G_rad_f;
            sum_Q += G_rad_f * T_rad_env;

            // E. 内部导热 (连接 Back 面)
            // Q_cond = K * (Tb_curr - Tf_curr) -> +K*Tb_curr - K*Tf_curr
            sum_G += node.conductance;
            sum_Q += node.conductance * Tb_curr; // 这里的 Tb_curr 是最新迭代值

            // F. 横向导热 (Lateral Conduction)
            // 逻辑：Q_lat = sum( K_ij * (T_neighbor - T_self) )
            // 在隐式方程中：
            // sum_Q += K_ij * T_neighbor
            // sum_G += K_ij

            for (const auto& link : node.neighbors) {
                int n_idx = link.neighbor_idx;
                double K_lat = link.conductance;

                // 获取邻居的最新温度 (使用 Front 温度，因为对于薄壳，横向导热通常看作均温或 Front 层主导)
                // 也可以更精细：Front 传 Front，Back 传 Back。这里简化为 Front 层横向导热。
                double T_neighbor = g_nodes[n_idx].T_front_next;

                sum_G += K_lat;
                sum_Q += K_lat * T_neighbor;
            }

            // >>> 求解 Front 新温度 <<<
            // sum_G * T_new = sum_Q
            node.T_front_next = sum_Q / sum_G;


            // ==========================
            // 节点 2: Back Surface
            // ==========================
            sum_G = 0.0;
            sum_Q = 0.0;

            // A. 热容
            sum_G += C_term;
            sum_Q += C_term * node.T_back;

            // B. 内部发热 (如有)
            double Q_gen = 0.0;
            sum_Q += Q_gen;

            // C. 对流
            auto conv_b = get_convection_params(node, false, w);
            double G_conv_b = conv_b.first * node.area;
            sum_G += G_conv_b;
            sum_Q += G_conv_b * conv_b.second;

            // D. 内部导热 (连接 Front 面)
            // 注意：这里要用刚刚更新过的 T_front_next，这叫 Gauss-Seidel 此时更新
            sum_G += node.conductance;
            sum_Q += node.conductance * node.T_front_next;

            // >>> 求解 Back 新温度 <<<
            node.T_back_next = sum_Q / sum_G;
        }
    }

    // 4. 更新状态
    for (auto& node : g_nodes) {
        node.T_front = node.T_front_next;
        node.T_back = node.T_back_next;
    }
}

// ==========================================
// 7. 可视化导出 (VTK 格式 - 用于 Paraview 查看云图)
// ==========================================
void export_vtk(const std::string& filename, double current_time) {
    std::ofstream file(filename);
    if (!file.is_open()) return;

    // 1. 统计总的顶点数和三角形数
    int total_points = 0;
    int total_cells = 0;
    for (const auto& node : g_nodes) {
        total_points += node.geometry_tris.size() * 3;
        total_cells += node.geometry_tris.size();
    }

    // 2. 写入 VTK 头部信息 (Legacy ASCII Format)
    file << "# vtk DataFile Version 3.0\n";
    file << "Thermal Simulation Result Time=" << current_time << "\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    // 3. 写入几何点坐标 (POINTS)
    file << "POINTS " << total_points << " float\n";
    for (const auto& node : g_nodes) {
        for (const auto& tri : node.geometry_tris) {
            // 写入三角形的三个顶点 v0, v1, v2
            file << tri.v0.x << " " << tri.v0.y << " " << tri.v0.z << "\n";
            file << tri.v1.x << " " << tri.v1.y << " " << tri.v1.z << "\n";
            file << tri.v2.x << " " << tri.v2.y << " " << tri.v2.z << "\n";
        }
    }

    // 4. 写入拓扑结构 (POLYGONS)
    file << "POLYGONS " << total_cells << " " << total_cells * 4 << "\n";
    int point_idx = 0;
    for (const auto& node : g_nodes) {
        for (size_t i = 0; i < node.geometry_tris.size(); ++i) {
            // 格式: 3(表示三角形) id0 id1 id2
            file << "3 " << point_idx << " " << point_idx + 1 << " " << point_idx + 2 << "\n";
            point_idx += 3;
        }
    }

    // 5. 写入温度数据 (CELL_DATA) - 这里相当于给每个三角形上色
    file << "CELL_DATA " << total_cells << "\n";
    file << "SCALARS Temperature float 1\n";
    file << "LOOKUP_TABLE default\n";

    for (const auto& node : g_nodes) {
        for (size_t i = 0; i < node.geometry_tris.size(); ++i) {
            // 只要是属于这个 Node 的三角形，都由该 Node 的温度决定颜色
            // 这里输出的是正面温度 (T_front)
            file << node.T_front << "\n";
        }
    }

    file.close();
    // std::cout << "Exported: " << filename << std::endl;
}
// ==========================================
// 主程序
// ==========================================
int main() {
    // 1. 加载配置 (材质和Group设置)
    setup_project_properties();

    // 2. 加载模型 (文件名来自配置)
    if (!init_model(g_config.global.obj_file)) return -1;

    build_topology();

    // 3. 加载天气 (文件名来自配置)
    if (!g_weather.load_weather(g_config.global.weather_file)) return -1;

    // 输出过滤
    std::vector<int> surf_indices;
    for (int i = 0; i < g_nodes.size(); ++i) {
        if (g_nodes[i].part_name.find("Deck") != std::string::npos ||
            g_nodes[i].part_name.find("Hull") != std::string::npos ||
            g_nodes[i].part_name.find("Stack") != std::string::npos ||
            g_nodes[i].part_name.find("Engine") != std::string::npos ||
            g_nodes[i].part_name.find("Room") != std::string::npos){ // 把 Engine 也输出来看看 Assigned 效果
            surf_indices.push_back(i);
        }
    }

    std::ofstream out_csv("results_advanced.csv");
    out_csv << "Time(h)";
    for (int idx : surf_indices) out_csv << "," << g_nodes[idx].part_name << "_N" << idx;
    out_csv << "\n";

    // 4. 使用配置的时间参数
    double start_h = g_config.global.start_time;
    double end_h = g_config.global.end_time;
    double dt = g_config.global.dt;

    std::cout << "\n[1/3] Initializing..." << std::endl;
    double start_sun_ang = 3.14159 * (start_h - 6.0) / 12.0;
    Vec3 start_sun = { std::cos(start_sun_ang), 0.5 * std::cos(start_sun_ang), std::sin(start_sun_ang) };
    if (g_weather.get_weather(start_h).solar <= 1.0) start_sun = { 0,0,-1 };
    update_shadows(start_sun);

    // 新逻辑：利用隐式法的大时间步特性，直接跳跃到稳态
// 只需要执行 1 次，或者为了保险起见执行 5-10 次即可
// 因为单次调用内部已经包含了 50 次 Gauss-Seidel 迭代
    for (int k = 0; k < 10; k++) {
        solve_step(dt, start_h, start_sun, true);

        // 可选：打印一下进度，让你知道它活着
        if (k % 2 == 0) std::cout << "." << std::flush;
    }
    std::cout << " -> Done." << std::endl;

    std::cout << "\n[2/3] Simulating..." << std::endl;
    int steps = (int)((end_h - start_h) * 3600.0 / dt);
    double last_output = -9999;

    // 创建一个文件夹存放结果是个好习惯，这里直接存根目录
    int frame_count = 0;

    for (int i = 0; i <= steps; i++) {
        double current_time = i * dt;
        double hour = start_h + (current_time / 3600.0);

        double sun_ang = 3.14159 * (hour - 6.0) / 12.0;
        Vec3 sun_dir = { std::cos(sun_ang), 0.5 * std::cos(sun_ang), std::sin(sun_ang) };
        if (g_weather.get_weather(hour).solar <= 1.0) sun_dir = { 0,0,-1 };

        if (i % (int)(900 / dt) == 0 && hour <= 19.0) update_shadows(sun_dir);

        solve_step(dt, hour, sun_dir, false);

        if (current_time - last_output >= 1200.0 - 0.1 || i == 0) {
            out_csv << std::fixed << std::setprecision(2) << hour;
            for (int idx : surf_indices) out_csv << "," << g_nodes[idx].T_front;
            out_csv << "\n";

            // VTK 输出
            std::string vtk_name = "sim_" + std::to_string(frame_count++) + ".vtk";
            export_vtk(vtk_name, hour);

            last_output = current_time;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << (double)i / steps * 100.0 << "%" << std::flush;
        }
    }

    std::cout << "\n[3/3] Done." << std::endl;
    return 0;
}