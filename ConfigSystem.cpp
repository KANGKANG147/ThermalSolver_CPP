#include "ConfigSystem.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

void ConfigSystem::init_defaults() {
    // 1. 初始化材质库 (这部分通常是固定的，或者也可以从文件读)
    //  材质定义
    /*struct Material {
        double k;           // 导热率
        double rho;         // 密度
        double Cp;          // 比热
        double absorptivity; // 太阳吸收率
        double emissivity;   // 红外发射率
    };*/
    mat_lib["Steel"] = { 52.019, 7768.98, 460.967, 0.74, 0.74 };
    mat_lib["Iron"] = { 40.0, 7800.0, 500.0, 0.8, 0.85 };
    mat_lib["Glass"] = { 0.81, 2800.0, 800.0, 0.1, 0.9 };
    mat_lib["Insulation"] = { 0.04, 50.0, 1000.0, 0.05, 0.95 };

    // 设置硬编码默认配置作为保底
    ConvectionBC bc_weather = { CONV_WIND, 0, 0, 10.0, 3.0 };
    ConvectionBC bc_cabin = { CONV_FIXED_H_T, 5.0, 0.0, 0, 0 };
    project_config["Default"] = { "Steel", 0.01, 20.0, 0.0, TYPE_CALCULATED, bc_weather, bc_cabin };
}

bool ConfigSystem::load_config(const std::string& filename) {
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
        if (key == "OBJ_FILE") ss >> this->settings.obj_file;
        else if (key == "WEATHER_FILE") ss >> this->settings.weather_file;
        else if (key == "SIM_START") {
            ss >> settings.start_date_time.year
                >> settings.start_date_time.month
                >> settings.start_date_time.day
                >> settings.start_date_time.hour;
        }
        else if (key == "SIM_END") {
            ss >> settings.end_date_time.year
                >> settings.end_date_time.month
                >> settings.end_date_time.day
                >> settings.end_date_time.hour;
        }
        else if (key == "DT") {
            ss >> settings.dt;
        }

        // --- Group 解析 ---
        else if (key == "BEGIN_GROUP") {
            inside_group = true;
            // 重置临时变量为默认值
            current_prop = { "Steel", 0.01, 20.0, 0.0, TYPE_CALCULATED,
                             {CONV_WIND, 0,0,10,3}, {CONV_FIXED_H_T, 5,0,0,0} };
            current_name = "Default";
        }
        else if (key == "END_GROUP") {
            if (inside_group && !current_name.empty()) {
                this->project_config[current_name] = current_prop;
                // std::cout << "  Configured Group: " << current_name << std::endl;
            }
            inside_group = false;
        }

        // 解析地理位置
        else if (key == "GEO_LOCATION") {
            ss >> this->settings.latitude
                >> this->settings.longitude
                >> this->settings.time_zone;
        }
        // 解析日期
        else if (key == "SIM_DATE") {
            ss >> this->settings.year
                >> this->settings.month
                >> this->settings.day;
        }
        // 解析模型朝向
        else if (key == "MODEL_HEADING") {
            ss >> this->settings.north_angle;
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
            else if (key == "HEAT_GEN_VOL") ss >> current_prop.volumetric_heat_gen;
            else if (key == "BC_FRONT") current_prop.front_bc = parse_bc(ss);
            else if (key == "BC_BACK")  current_prop.back_bc = parse_bc(ss);
        }
    }
    return true;
}

// 注意：我把 init_model 的逻辑移到了这里，并稍微改造了一下参数
bool ConfigSystem::load_obj_model(const std::string& filename, std::vector<ThermalNode>& out_nodes) {
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
            Material mat = this->mat_lib[prop.material_name];

            ThermalNode node;
            node.part_name = group;
            node.area = area;
            node.geometry_tris = tris;
            Vec3 center_sum = { 0,0,0 }; 
            for (auto& t : tris) center_sum = center_sum + t.v0 + t.v1 + t.v2;
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

            // 计算该节点的总发热功率 (W) = (W/m^3) * Area * Thickness 
            node.Q_gen_total = prop.volumetric_heat_gen * node.area * node.thickness;

            // 初始温度：如果是 Assigned，这个值将作为永久固定值
            node.T_front = node.T_back = prop.initial_temp;
            node.shadow_factor = 1.0;
            out_nodes.push_back(node);
        }
    }
    std::cout << "Nodes: " << out_nodes.size() << std::endl;
    return true;
}

void ConfigSystem::export_vtk(const std::string& filename, double current_time, const std::vector<ThermalNode>& nodes) {
    std::ofstream file(filename);
    if (!file.is_open()) return;

    // 1. 统计总的顶点数和三角形数
    int total_points = 0;
    int total_cells = 0;
    for (const auto& node : nodes) {
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
    for (const auto& node : nodes) {
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
    for (const auto& node : nodes) {
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

    for (const auto& node : nodes) {
        for (size_t i = 0; i < node.geometry_tris.size(); ++i) {
            // 只要是属于这个 Node 的三角形，都由该 Node 的温度决定颜色
            // 这里输出的是正面温度 (T_front)
            file << node.T_front << "\n";
        }
    }

    file.close();
    // std::cout << "Exported: " << filename << std::endl;
}

void ConfigSystem::export_results_tai_format(const std::string& filename, const std::vector<ThermalNode>& nodes) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing." << std::endl;
        return;
    }

    // 1. 设置高精度输出 
    out << std::fixed << std::setprecision(14);

    std::string current_group = "";

    // 2. 设置高精度输出
    out << std::fixed << std::setprecision(14);

    // 3. 遍历所有节点（面片）输出温度
    for (const auto& node : nodes) {
        // 3. 检测组名是否发生变化
        // 如果当前节点的 part_name 与上一个不同，说明进入了新 Group
        if (node.part_name != current_group) {
            current_group = node.part_name;

            // 如果 OBJ 中未定义组名，代码默认为 "Default"，此处可直接输出
            // 格式示例: g Body 或 g Hull
            out << "g " << current_group << "\n";
        }

        // 4. 输出温度值 (格式: f <temp>)
        out << "f " << node.T_front << "\n";
    }

    out.close();
    std::cout << "[Export] Results exported to " << filename << " (tai format)" << std::endl;
}

// 私有辅助函数
PartProperty ConfigSystem::get_part_property(const std::string& group_name) {
    if (project_config.count(group_name)) return project_config[group_name];
    for (const auto& kv : project_config) {
        if (group_name.find(kv.first) != std::string::npos && kv.first != "Default") return kv.second;
    }
    return project_config["Default"];
}
ConvectionBC ConfigSystem::parse_bc(std::stringstream& ss) {
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