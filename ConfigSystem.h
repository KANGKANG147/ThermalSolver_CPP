#pragma once
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include "CoreTypes.h"

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
    void export_vtk(const std::string& filename, double current_time, const std::vector<ThermalNode>& nodes);

private:
    ConvectionBC parse_bc(std::stringstream& ss);
    PartProperty get_part_property(const std::string& group_name);
};