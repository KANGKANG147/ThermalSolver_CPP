#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "MathUtils.h"
#include "ConfigSystem.h"
#include "WeatherSystem.h"
#include "ThermalSolver.h"
#include "SolarSystem.h"

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
   
    // 5. 稳态初始化 (预热)
    std::cout << "\n[1/3] Pre-heating (Steady State Init)..." << std::endl;
    DateTime curr_time = config.settings.start_date_time; // 这里的 curr_time 是动态的

    SolarVector start_sun = SolarSystem::calculate_model_sun(
        curr_time.year, curr_time.month, curr_time.day,
        curr_time.hour,
        config.settings.time_zone,
        config.settings.latitude, config.settings.longitude,
        config.settings.north_angle
    );

    // 获取天气数据 
    // 逻辑：稳态初始化直接取仿真开始时刻对应的天气
    // 假设 Weather File 的 Day 1 = Simulation Start Day
    // Weather Query Hour = (0天 * 24) + current_hour
    double init_weather_query = curr_time.hour;
    WeatherData w_init = weather.get_weather(init_weather_query);

    Vec3 init_sun_dir = start_sun.dir;
    // 如果 SPA 说是晚上，或者天气文件说辐射极低
    if (start_sun.is_night || w_init.solar <= 1.0) {
        init_sun_dir = { 0, 0, -1 };
    }

    // 更新一次初始阴影
    if (!start_sun.is_night) solver.update_shadows(init_sun_dir);

    //稳态初始化
    solver.solve_step(config.settings.dt, init_weather_query, init_sun_dir, weather, true);
    std::cout << " Done." << std::endl;

    // 6. 瞬态模拟循环
    std::cout << "\n[2/3] Simulating Transient..." << std::endl;

    std::ofstream out_csv("results.csv");
    out_csv << "Time(h),Node_0_Temp_Front,Node_0_Temp_Back\n";

    // 计算总步数 (基于秒数差)
    double total_seconds = get_duration_seconds(config.settings.start_date_time, config.settings.end_date_time);
    int steps = (int)(total_seconds / config.settings.dt);

    // 如果 start > end，steps 会是负数，提示错误
    if (steps <= 0) {
        std::cerr << "Error: Start time is later than End time!" << std::endl;
        return -1;
    }

    double last_output_time = -9999;
    int frame_count = 0;

    // 记录开始时的 Day of Year，用于计算经过的天数
    int start_doy = config.settings.start_date_time.get_day_of_year();
    int start_year = config.settings.start_date_time.year;

    for (int i = 0; i <= steps; i++) {
        double elapsed_sec = i * config.settings.dt;
        // ==========================================
        // ★★★ 核心修改：使用 NREL SPA 计算太阳 ★★★
        // ==========================================
        SolarVector sun = SolarSystem::calculate_model_sun(
            config.settings.year,
            config.settings.month,
            config.settings.day,
            curr_time.hour,
            config.settings.time_zone,
            config.settings.latitude,
            config.settings.longitude,
            config.settings.north_angle
        );

        // 2. 计算天气查询时间 (Weather Query Time)
        // 算法： (当前年的DOY - 开始年的DOY) * 24.0 + 当前小时
        // 这将生成一个连续增加的时间：6.0 -> 24.0 -> 30.0 (第二天6点)
        // 正好对应 WeatherSystem 解析出来的连续数据。

        int current_doy = curr_time.get_day_of_year();
        // 加上年份带来的天数差 (防止跨年出错)
        int day_diff = (current_doy + (curr_time.year - start_year) * 365) - start_doy;

        double weather_query_hour = (day_diff * 24.0) + curr_time.hour;

        // 3. 查天气
        // WeatherSystem 内部存储的是类似 5.9, 23.9, 24.0, 48.0 的数据
        // 所以我们传 30.0 进去，它就能找到第二天的数据，而不是第一天的。
        WeatherData w_curr = weather.get_weather(weather_query_hour);

        // 最终太阳向量
        Vec3 sun_dir = sun.dir;

        // 如果是晚上，或者天气文件说没太阳，或者 SPA 算出在地平线以下
        if (sun.is_night || w_curr.solar <= 1.0) {
            sun_dir = { 0, 0, -1 }; // 强制设为无效方向
        }

        // 更新阴影 (仅白天且每隔一段时间)
        if (!sun.is_night && (config.settings.dt >= 900.0 || (i % (int)(900 / config.settings.dt) == 0))) {
            solver.update_shadows(sun_dir);
        }

        // 求解一步
        solver.solve_step(config.settings.dt, weather_query_hour, sun_dir, weather, false);

        // 输出结果 (每20分钟)
        if (elapsed_sec - last_output_time >= 1200.0 - 0.1 || i == 0) {
            // CSV 记录第0个节点的温度
            if (!solver.nodes.empty()) {
                out_csv << std::fixed << std::setprecision(2) 
                    << elapsed_sec / 3600.0 << ","
                    << curr_time.year << "," << curr_time.month << "," << curr_time.day << "," << curr_time.hour << ","
                    << solver.nodes[0].T_front << "\n";
            }

            // VTK
            std::string vtk_name = "sim_" + std::to_string(frame_count++) + ".vtk";
            config.export_vtk(vtk_name, weather_query_hour, solver.nodes);

            last_output_time = elapsed_sec;
            std::cout << "\rProgress: " << (int)((double)i / steps * 100.0) << "% "
                << curr_time.year << "-" << curr_time.month << "-" << curr_time.day << " "
                << std::fixed << std::setprecision(1) << curr_time.hour << "h" << std::flush;
        }
        // 时间推进
        advance_time(curr_time, config.settings.dt);
    }

    std::cout << "\n[3/3] Simulation Complete." << std::endl;
    out_csv.close();
    return 0;
}