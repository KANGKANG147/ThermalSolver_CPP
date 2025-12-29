#pragma once
#include <cmath>
#include <iostream>
#include "MathUtils.h"

// 因为 spa.h 是 C 语言写的，在 C++ 里引用需要加 extern "C"
extern "C" {
    #include "spa.h"
}

struct SolarVector {
    Vec3 dir;       // 模型坐标系下的单位光照向量
    double azimuth; // 地理方位角 (度)
    double zenith;  // 天顶角 (度)
    bool is_night;  // 是否是晚上
};

class SolarSystem {
public:
    // 计算太阳向量
    // 参数：
    // date: 年月日
    // time_hour: 当地时间 (例如 14.5 表示 14:30)
    // time_zone: 时区 (中国是 +8)
    // lat/lon: 纬度/经度
    // north_angle: 模型Y轴相对于正北的角度 (模型朝向)
    static SolarVector calculate_model_sun(int year, int month, int day,
        double time_hour, double time_zone,
        double lat, double lon,
        double north_angle)
    {
        spa_data spa; // 声明 NREL SPA 的数据结构

        // 1. 初始化基础参数
        spa.year = year;
        spa.month = month;
        spa.day = day;
        spa.timezone = time_zone;

        // 把小数时间 (14.5) 拆解为 时:分:秒
        spa.hour = (int)time_hour;
        double rem_min = (time_hour - spa.hour) * 60.0;
        spa.minute = (int)rem_min;
        spa.second = (int)((rem_min - spa.minute) * 60.0);

        // 地理参数
        spa.longitude = lon;
        spa.latitude = lat;
        spa.elevation = 0; // 海拔 (米)，影响很小，暂设0

        // 2. 气象修正参数 (用于计算大气折射)
        // Delta T: 地球自转与原子时的差值，2024年约 69秒，精度要求不高可用默认
        spa.delta_t = 67.0;
        spa.pressure = 1013.25; // 标准大气压 (mbar)
        spa.temperature = 15.0; // 平均气温 (C)

        // slope 和 azm_rotation 用于计算坡面辐射，我们只求空间向量，设为0即可
        spa.slope = 0;
        spa.azm_rotation = 0;
        spa.atmos_refract = 0.5667; // 标准大气折射率
        spa.function = SPA_ZA;      // 我们只需要 Zenith 和 Azimuth，不需要算辐射值

        // 3. 调用 C 语言库函数
        int res = spa_calculate(&spa);

        if (res != 0) {
            std::cerr << "[Solar Error] SPA calculation failed. Code: " << res << std::endl;
            return { {0,0,-1}, 0, 0, true }; // 报错返回无效值
        }

        // 4. 构造输出
        SolarVector out;
        out.azimuth = spa.azimuth;
        out.zenith = spa.zenith;
        std::cout << "azimuth" << spa.azimuth<< "zenith" << spa.zenith << std::endl;
        // 判断是否是夜晚
        // Zenith = 0 (头顶), 90 (地平线). 大于90即太阳下山
        // 考虑一点大气折射余量 (90.5度)
        if (spa.zenith > 90.5) {
            out.is_night = true;
            out.dir = { 0, 0, -1 }; // 指向地下
            return out;
        }
        out.is_night = false;

        // 5. 【核心】坐标变换：从地平坐标系 -> 模型坐标系

        // 将角度转弧度
        double deg2rad = 3.14159265359 / 180.0;
        double az_rad = spa.azimuth * deg2rad;
        double ze_rad = spa.zenith * deg2rad;
        double heading_rad = north_angle * deg2rad;

        // 计算相对方位角 (Relative Azimuth)
        // 太阳相对于模型正前方的夹角
        double alpha = az_rad - heading_rad;

        // 球坐标转笛卡尔坐标 (Z轴向上, Y轴向前, X轴向右)
        // S_z = cos(Zenith)
        // S_horizontal = sin(Zenith)
        // S_x = S_horizontal * sin(alpha)
        // S_y = S_horizontal * cos(alpha)

        out.dir.z = std::cos(ze_rad);
        double horizontal_comp = std::sin(ze_rad);
        out.dir.x = horizontal_comp * std::sin(alpha);
        out.dir.y = horizontal_comp * std::cos(alpha);

        // 归一化 (防误差)
        out.dir = normalize(out.dir);

        return out;
    }
};