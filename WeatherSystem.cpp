#include "WeatherSystem.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>

bool WeatherSystem::load_weather(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "[Error] Cannot open weather file: " << filename << std::endl;
        return false;
    }

    // 1. 读取整个文件到字符串
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();

    std::cout << "[Weather] File read. Size: " << content.size() << " chars." << std::endl;

    // 2. 找到数据起始点 (TIME 之后)
    // 这样可以跳过文件开头的 "04" 这种无关数字
    size_t data_start = content.find("TIME");
    if (data_start != std::string::npos) {
        content = content.substr(data_start);
    }

    // 3. 正则清洗 (这是最关键的一步)
    try {
        // A. 去除 这种标签
        std::regex re_source("\\[.*?\\]");
        content = std::regex_replace(content, re_source, " ");

        // B. 去除所有字母 (TIME, AIRT 等标题)
        std::regex re_letters("[A-Za-z]");
        content = std::regex_replace(content, re_letters, " ");

        // C. 将冒号等标点替换为空格 (防止 09:20 读错)
        std::regex re_punct("[:]");
        content = std::regex_replace(content, re_punct, " ");
    }
    catch (const std::regex_error& e) {
        std::cout << "[Error] Regex failed: " << e.what() << std::endl;
        return false;
    }

    // 4. 流式读取纯数字
    std::stringstream ss(content);
    std::vector<double> all_numbers;
    double val;
    while (ss >> val) {
        all_numbers.push_back(val);
    }

    std::cout << "[Weather] Extracted " << all_numbers.size() << " numbers." << std::endl;

    // 5. 按列重组 (你的文件是 8 列数据)
    // 0:TIME, 1:AIRT, 2:SOLAR, 3:WIND, 4:HUMID, 5:CLOUD, 6:LWIR, 7:WINDIR(or Rain)
    int stride = 8;

    if (all_numbers.size() < stride) {
        std::cout << "[Error] Not enough data found!" << std::endl;
        return false;
    }

    double last_raw_time = -1.0;
    double day_offset = 0.0;

    for (size_t i = 0; i + stride <= all_numbers.size(); i += stride) {
        double raw_time = all_numbers[i];     // HHMM
        double airt = all_numbers[i + 1];
        double solar = all_numbers[i + 2];
        double wind = all_numbers[i + 3];
        double lwir = all_numbers[i + 6];
        // 格式解析 HHMM -> Hour
        int hh = (int)(raw_time / 100);
        int mm = (int)(raw_time) % 100;

        // 基本校验
        if (hh > 24 || mm >= 60) continue; // 跳过非法时间

        // 处理跨天 (2355 -> 0000)
        if (last_raw_time >= 0 && raw_time < last_raw_time) {
            // 如果时间突然变小 (且不是一点点波动)，说明过了一天
            // 比如 2355 -> 0000
            day_offset += 24.0;
        }
        last_raw_time = raw_time;

        double final_hour = hh + mm / 60.0 + day_offset;

        data_points.push_back({ final_hour, airt, solar, wind, lwir });
    }

    std::cout << "[Weather] Successfully parsed " << data_points.size() << " time steps." << std::endl;
    if (!data_points.empty()) {
        std::cout << "          Range: " << data_points.front().time_hour << "h -> "
            << data_points.back().time_hour << "h" << std::endl;
    }

    return true;
}

WeatherData WeatherSystem::get_weather(double query_hour) {
    if (data_points.empty()) return { query_hour, 20.0, 0.0, 1.0 };

    if (query_hour <= data_points.front().time_hour) return data_points.front();
    if (query_hour >= data_points.back().time_hour) return data_points.back();

    // 线性查找插值
    for (size_t i = 0; i < data_points.size() - 1; ++i) {
        if (query_hour >= data_points[i].time_hour && query_hour < data_points[i + 1].time_hour) {
            const auto& p1 = data_points[i];
            const auto& p2 = data_points[i + 1];

            double ratio = (query_hour - p1.time_hour) / (p2.time_hour - p1.time_hour);

            WeatherData res;
            res.time_hour = query_hour;
            res.air_temp = p1.air_temp + ratio * (p2.air_temp - p1.air_temp);
            res.solar = p1.solar + ratio * (p2.solar - p1.solar);
            res.wind = p1.wind + ratio * (p2.wind - p1.wind);
            res.lwir = p1.lwir + ratio * (p2.lwir - p1.lwir);
            return res;
        }
    }
    return data_points.back();
}