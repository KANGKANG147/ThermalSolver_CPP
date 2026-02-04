#pragma once
#include <vector>
#include "CoreTypes.h"
#include "WeatherSystem.h" 
#include "BVH.h"
#include "LinearAlgebra.h"
#include "EigenSolverAdapter.h"
#include "SolarRadiation.h"
#include "ConfigSystem.h"
#include <fstream>
#include <filesystem>

class ThermalSolver {
public:
    std::vector<ThermalNode> nodes; // 所有的热节点
    BVHAccel bvh;

    EigenSolverAdapter eigen_solver; // 求解器实例

    // --- 背景环境参数 ---
    bool enable_background = true;      // 开关
    BackgroundType bg_type = BG_SEA;    // 类型
    double sea_temp_K = 288.15;         // 海水温度初始化为 288.15K (15C)
    double ground_temp_K = 293.15;      // 地面温度
    double sea_albedo = 0.1;            // 海面反照率
    double ground_albedo = 0.2;         // 地面反照率

    // 设置背景参数的接口
    void set_background_params(bool enable, int type, double water_C, double ground_C, double s_albedo, double g_albedo) {
        enable_background = enable;
        bg_type = (BackgroundType)type;
        sea_temp_K = water_C + 273.15;
        ground_temp_K = ground_C + 273.15;
        sea_albedo = s_albedo;
        ground_albedo = g_albedo;
    }

    // 场景特征尺度 (用于计算自适应 Bias)
    // 初始化为 0.0，表示还没计算过
    double scene_scale = 0.0;

    // 辅助函数：计算场景包围盒和尺度
    void ensure_scene_scale();

    // 构建拓扑连接（顶点焊接、横向导热）
    void build_topology();

    // 解析耦合关系：将 PartName 映射到 node index
    void resolve_couplings();

    void enforce_reciprocity();

    // 更新阴影遮挡
    void update_shadows(const Vec3& sun_dir);

    // MCRT 核心入口：计算角系数
    void calculate_view_factors(int samples);

    // 求解辐射度矩阵 (同时处理 Front 和 Back)
    void solve_radiosity_system(double sky_temp_K);

    // 核心计算步
    void solve_step(double dt, double hour, const Vec3& sun_dir, 
                    double zenith_deg, int day_of_year,
                    WeatherSystem& weather, bool is_steady_init);

private:
    // 内部辅助函数
    std::pair<double, double> get_convection_params(const ThermalNode& node, bool is_front, const WeatherData& w);
    double calc_h_rad(double T_surf_K, double T_env_K, double epsilon);

    // 缓存系统辅助函数
    size_t compute_geometry_checksum(int samples);
    bool load_vf_cache(const std::string& filename, int samples);
    void save_vf_cache(const std::string& filename, int samples);
};
