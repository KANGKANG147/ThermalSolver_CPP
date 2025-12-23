#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "MathUtils.h"
#include "ConfigSystem.h"
#include "WeatherSystem.h"
#include "ThermalSolver.h"

// ==========================================
// 主程序
// ==========================================
int main() {
    // 1. 初始化系统
    ConfigSystem config;
    WeatherSystem weather;
    ThermalSolver solver;

    // 2. 加载配置和数据
    config.init_defaults();
    config.load_config("config.txt"); // 如果找不到，会用代码里的默认值

    // 注意：这里把 config 和 solver 连起来了，加载的模型直接放进 solver.nodes
    if (!config.load_obj_model(config.settings.obj_file, solver.nodes)) {
        std::cerr << "Failed to load OBJ model: " << config.settings.obj_file << std::endl;
        return -1;
    }

    // 3. 加载天气 (文件名来自配置)
    if (!weather.load_weather(config.settings.weather_file)) {
        std::cerr << "Failed to load weather file: " << config.settings.weather_file << std::endl;
        return -1;
    }

    // 3. 建立物理拓扑
    solver.build_topology();

    // 计算角系数 (预计算)
    // 建议采样数 1000 以上，越多越准，但越慢
    solver.calculate_view_factors(1000);
    // 输出过滤
  /*  std::vector<int> surf_indices;
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
    */
    // 4. 准备模拟参数
    double start_h = config.settings.start_time;
    double end_h = config.settings.end_time;
    double dt = config.settings.dt;

    // 5. 稳态初始化 (预热)
    std::cout << "\n[1/3] Pre-heating (Steady State Init)..." << std::endl;
    double start_sun_ang = 3.14159 * (start_h - 6.0) / 12.0;
    Vec3 start_sun = { std::cos(start_sun_ang), 0.5 * std::cos(start_sun_ang), std::sin(start_sun_ang) };
    if (weather.get_weather(start_h).solar <= 1.0) start_sun = { 0,0,-1 };

    solver.update_shadows(start_sun);

    // 预热 10 步，每步内部迭代 50 次
    for (int k = 0; k < 10; k++) {
        solver.solve_step(dt, start_h, start_sun, weather, true);
    }
    std::cout << " Done." << std::endl;

    // 6. 瞬态模拟循环
    std::cout << "\n[2/3] Simulating Transient..." << std::endl;

    std::ofstream out_csv("results.csv");
    out_csv << "Time(h),Node_0_Temp_Front,Node_0_Temp_Back\n";

    int steps = (int)((end_h - start_h) * 3600.0 / dt);
    double last_output_time = -9999;
    int frame_count = 0;

    for (int i = 0; i <= steps; i++) {
        double current_sec = i * dt;
        double hour = start_h + (current_sec / 3600.0);

        // 计算太阳位置
        double sun_ang = 3.14159 * (hour - 6.0) / 12.0;
        Vec3 sun_dir = { std::cos(sun_ang), 0.5 * std::cos(sun_ang), std::sin(sun_ang) };
        if (weather.get_weather(hour).solar <= 1.0) sun_dir = { 0,0,-1 };

        // 更新阴影 (策略优化：如果步长很大，建议每步都更)
        // 如果 dt > 900s，则每步都更；否则每15分钟更
        if (dt >= 900.0 || (i % (int)(900 / dt) == 0)) {
            if (hour <= 19.0) solver.update_shadows(sun_dir);
        }

        // 求解一步
        solver.solve_step(dt, hour, sun_dir, weather, false);

        // 输出结果 (每20分钟)
        if (current_sec - last_output_time >= 1200.0 - 0.1 || i == 0) {
            // CSV 记录第0个节点的温度
            if (!solver.nodes.empty()) {
                out_csv << std::fixed << std::setprecision(2) << hour << "," << solver.nodes[0].T_front << "\n";
            }

            // VTK
            std::string vtk_name = "sim_" + std::to_string(frame_count++) + ".vtk";
            config.export_vtk(vtk_name, hour, solver.nodes);

            last_output_time = current_sec;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << (double)i / steps * 100.0 << "%" << std::flush;
        }
    }

    std::cout << "\n[3/3] Simulation Complete." << std::endl;
    out_csv.close();
    return 0;
}