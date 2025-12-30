#include "SolarRadiation.h"
#include <iostream> // 用于可能的调试输出

// 常量定义
static const double PI = 3.14159265358979323846;
static const double SOLAR_CONSTANT = 1367.0; // 太阳常数 W/m^2

double SolarRadiation::calc_extraterrestrial_irradiance(int doy) {
    // 地球轨道偏心率修正因子 (Eccentricity correction factor)
    double B = (doy - 1) * 360.0 / 365.0; // 角度 B (度)
    double B_rad = B * PI / 180.0;

    // I0 = Gsc * (1 + 0.033 * cos(360n/365)) 
    // 这里使用更精确的 Spencer (1971) 级数展开形式
    double E0 = 1.000110 + 0.034221 * std::cos(B_rad) + 0.001280 * std::sin(B_rad)
        + 0.000719 * std::cos(2.0 * B_rad) + 0.000077 * std::sin(2.0 * B_rad);

    return SOLAR_CONSTANT * E0;
}

SolarComponents SolarRadiation::split_GHI_Erbs(double GHI, double zenith_deg, int doy) {
    SolarComponents result = { 0.0, 0.0, GHI, 0.0 };

    // 1. 物理边界过滤：如果没有光，或者太阳快下山了(Zenith > 87度)，直接当作无辐射
    // 避免除零和由于余弦效应造成的巨大误差
    if (GHI < 1.0 || zenith_deg > 89.5) {
        return result; // 全零
    }

    double zenith_rad = zenith_deg * PI / 180.0;
    double cos_theta_z = std::cos(zenith_rad);

    // 钳位 cos_theta_z 防止除零 (对应 Zenith > 88.8度)
    if (cos_theta_z < 0.02) cos_theta_z = 0.02;

    // 2. 计算大气层外辐射 I0
    double I0 = calc_extraterrestrial_irradiance(doy);

    // 3. 计算晴空指数 kt (Clearness Index)
    // kt = GHI / (I0 * cos(theta_z))
    // 代表有多少比例的太阳光穿透了大气层
    double kt = GHI / (I0 * cos_theta_z);

    // 物理约束：kt 理论上在 0 到 1 之间 (极少数晴朗干燥高海拔可能略超 1)
    if (kt < 0.0) kt = 0.0;
    if (kt > 1.0) kt = 1.0;
    result.kt = kt;

    // 4. Erbs 相关式计算漫射比 kd (Diffuse Fraction = DHI / GHI)
    double kd = 0.0;

    if (kt <= 0.22) {
        // 阴天：几乎全是漫射
        kd = 1.0 - 0.09 * kt;
    }
    else if (kt <= 0.80) {
        // 多云到晴：混合区域 (Erbs 核心公式)
        kd = 0.9511 - 0.1604 * kt + 4.388 * std::pow(kt, 2)
            - 16.638 * std::pow(kt, 3) + 12.336 * std::pow(kt, 4);
    }
    else {
        // 非常晴朗：漫射比例很低
        kd = 0.165;
    }

    // 5. 分离计算
    // DHI (水平漫射)
    result.DHI = GHI * kd;

    // DNI (法向直射)
    // 几何关系：GHI = DNI * cos(z) + DHI
    // -> DNI = (GHI - DHI) / cos(z)
    result.DNI = (GHI - result.DHI) / cos_theta_z;

    // 最后的物理安全检查
    if (result.DNI < 0.0) result.DNI = 0.0;
    // DNI 不太可能超过大气层外辐射太多
    if (result.DNI > I0) result.DNI = I0;

    return result;
}