#pragma once
#include <vector>
#include <string>
#include "CoreTypes.h"

struct WeatherData {
    double time_hour;
    double air_temp;
    double solar;
    double wind;
    double cloud;// Cloud Cover [0-10]
    double lwir;
};

class WeatherSystem {
public:
    bool load_weather(const std::string& filename);
    WeatherData get_weather(double query_hour);

private:
    std::vector<WeatherData> data_points;
};