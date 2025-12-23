#pragma once
#include <vector>
#include "CoreTypes.h"
#include "WeatherSystem.h" 
#include "BVH.h"
#include "LinearAlgebra.h"
#include "AMGSolver.h"

class ThermalSolver {
public:
    std::vector<ThermalNode> nodes; // 所有的热节点
    BVHAccel bvh;

    // AMG 求解器实例
    AMGSolver amg_solver;
    bool amg_initialized = false;

    // 构建拓扑连接（顶点焊接、横向导热）
    void build_topology();

    // 更新阴影遮挡
    void update_shadows(const Vec3& sun_dir);

    // MCRT 核心入口：计算角系数
    void calculate_view_factors(int samples);

    // 求解辐射度矩阵 (同时处理 Front 和 Back)
    void solve_radiosity_system(double env_temp_K);

    // 核心计算步
    void solve_step(double dt, double hour, const Vec3& sun_dir, WeatherSystem& weather, bool is_steady_init);

private:
    // 内部辅助函数
    std::pair<double, double> get_convection_params(const ThermalNode& node, bool is_front, const WeatherData& w);
    double calc_h_rad(double T_surf_K, double T_env_K, double epsilon);
};
