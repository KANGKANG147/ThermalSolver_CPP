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