#pragma once
#include <cmath>
#include <algorithm>

// 辐射分量结构体
struct SolarComponents {
    double DNI; // 法向直射辐射 (Direct Normal Irradiance) [W/m^2]
    double DHI; // 水平漫射辐射 (Diffuse Horizontal Irradiance) [W/m^2]
    double GHI; // 水平总辐射 (Global Horizontal Irradiance) [W/m^2]
    double kt;  // 晴空指数 (Clearness Index) [0-1]
};

class SolarRadiation {
public:
    /**
     * @brief 使用 Erbs 模型将水平总辐射(GHI)分离为直射(DNI)和漫射(DHI)
     * * @param GHI        水平总辐射 (来自气象文件) [W/m^2]
     * @param zenith_deg 太阳天顶角 [度]
     * @param doy        当前是一年中的第几天 (Day of Year) [1-366]
     * @return SolarComponents 分离后的分量
     */
    static SolarComponents split_GHI_Erbs(double GHI, double zenith_deg, int doy);

private:
    // 计算大气层外太阳辐射 (Extraterrestrial Radiation)
    static double calc_extraterrestrial_irradiance(int doy);
};
