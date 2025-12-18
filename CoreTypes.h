#pragma once
#include <string>
#include <vector>
#include "MathUtils.h"

//  材质定义
struct Material {
    double k;           // 导热率
    double rho;         // 密度
    double Cp;          // 比热
    double absorptivity; // 太阳吸收率
    double emissivity;   // 红外发射率
};

//  对流/边界类型枚举
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

//  Group 类型枚举
enum GroupType {
    TYPE_CALCULATED,    // 计算温度 (正常物理对象)
    TYPE_ASSIGNED       // ★★★ 新增：指定温度 (恒温源) ★★★
};

//  部件属性配置
struct PartProperty {
    std::string material_name;
    double thickness;
    double initial_temp; // 如果是 Assigned 类型，这个值就是全过程的固定温度

    // 新增：体积热源强度 (W/m^3)
    double volumetric_heat_gen = 0.0;

    GroupType group_type; // Assigned or Calculated

    ConvectionBC front_bc;
    ConvectionBC back_bc;
};

// --- 网格节点 ---
struct Triangle { Vec3 v0, v1, v2; };
struct EdgeLink {
    int neighbor_idx;
    double conductance;
};

// 辐射链接：记录“我”看到了“谁”的“哪一面”
struct RadLink {
    int target_node_idx;     // 目标 Node ID
    bool target_is_front;    // true=看到的是目标的正面, false=看到的是背面
    double view_factor;      // 角系数
};

struct ThermalNode {
    double T_front, T_back;
    double T_front_next, T_back_next;

    double mass_node, conductance, area;
    double solar_absorp, ir_emissivity;

    // 新增：该节点的总内部热源功率 (Watts)
    // 注意：求解时 Front 和 Back 各分一半
    double Q_gen_total = 0.0;

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

    // --- MCRT ---
    // 正面 (Front) 的辐射数据
    std::vector<RadLink> rad_links_front;
    double vf_sky_front; // 正面看到天空的系数

    // 背面 (Back) 的辐射数据
    std::vector<RadLink> rad_links_back;
    double vf_sky_back;  // 背面看到天空的系数

    // 辐射度 J (W/m^2): 离开表面的总能量
    double J_front = 0.0;
    double J_back = 0.0;

    // 净辐射热流 Q_rad (Watts): 最终算出来用于加热节点的能量
    double Q_rad_front = 0.0;
    double Q_rad_back = 0.0;
};