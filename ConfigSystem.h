#pragma once
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include "CoreTypes.h"

// 背景类型枚举
enum BackgroundType {
    BG_SEA = 0,
    BG_GROUND = 1
};

struct GlobalSettings {
    std::string obj_file = "chuan.tai";
    std::string weather_file = "weather.txt";
    DateTime start_date_time;
    DateTime end_date_time;
    double dt = 300.0;

    // [新增] 世界坐标系参数
    double latitude = 39.90;    // 默认北京
    double longitude = 116.40;
    double time_zone = 8.0;     // UTC+8
    double north_angle = 0.0;   // 模型Y轴就是正北

    // --- 背景控制参数 ---
    bool enable_background = true;      // 总开关：是否考虑背景
    BackgroundType background_type = BG_SEA; // 背景类型：海面 或 地面

    double water_temp = 15.0;           // 海水温度 (C)
    double sea_albedo = 0.1;            // 海面反照率 (0.0 - 1.0)

    double ground_temp = 20.0;          // [新增] 地面温度 (C)
    double ground_albedo = 0.2;         // [新增] 地面反照率 (0.0 - 1.0)

    // [新增] 模拟日期
    int year = 2024;
    int month = 7;
    int day = 21;
};

class ConfigSystem {
public:
    GlobalSettings settings;
    std::map<std::string, Material> mat_lib;
    std::map<std::string, PartProperty> project_config;

    void init_defaults();
    bool load_config(const std::string& filename);
    bool load_obj_model(const std::string& filename, std::vector<ThermalNode>& out_nodes);
    // 创建配置中定义的流体节点
    void create_fluid_nodes(std::vector<ThermalNode>& out_nodes);
    void export_vtk(const std::string& filename, double current_time, const std::vector<ThermalNode>& nodes);
    void export_results_tai_format(const std::string& filename, const std::vector<ThermalNode>& nodes);
private:
    ConvectionBC parse_bc(std::stringstream& ss);
    PartProperty get_part_property(const std::string& group_name);
};