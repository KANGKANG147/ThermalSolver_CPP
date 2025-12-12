#pragma once
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include "CoreTypes.h"

struct GlobalSettings {
    std::string obj_file = "chuan.tai";
    std::string weather_file = "weather.txt";
    double start_time = 6.0;
    double end_time = 16.0;
    double dt = 0.5;
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