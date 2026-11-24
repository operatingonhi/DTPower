#ifndef CONFIG_LOADER_H
#define CONFIG_LOADER_H

#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>

struct EngineGeometry {
    double bore;
    double stroke;
    double connecting_rod_length;
    double compression_ratio;
};

struct OperatingConditions {
    double engine_speed;
    double intake_pressure;
    double intake_temperature;
    double exhaust_pressure;
    std::string intake_composition;
    double wall_temperature;
};

struct SparkParameters {
    double timing_btdc;
    double energy;
    double duration;
    bool enabled;
};

struct ValveParameters {
    double open_angle;
    double close_angle;
    double max_lift;
    double diameter;
    double seat_diameter;
    double discharge_coefficient;
};

struct ValveConfig {
    ValveParameters intake;
    ValveParameters exhaust;
};

struct EngineParameters {
    EngineGeometry geometry;
    OperatingConditions operating_conditions;
    SparkParameters spark;
    ValveConfig valves;
    std::string mechanism_file;
};

struct TimeSteppingParameters {
    double crank_angle_step;
    double min_angle_step;
    double max_angle_step;
    bool adaptive_stepping;
};

struct HeatTransferParameters {
    double woschni_c1;
    double woschni_c2;
    double woschni_c3;
};

struct FrictionParameters {
    double cf_c1;
    double cf_c2;
    double cf_c3;
    double cf_c4;
};

struct SimulationControl {
    double output_interval;
    double start_angle;
    double end_angle;
    bool equilibrium_expansion = true;
};

struct ToleranceParameters {
    double temperature_min;
    double temperature_max;
    double pressure_min;
    double pressure_max;
};

struct SimulationParameters {
    TimeSteppingParameters time_stepping;
    HeatTransferParameters heat_transfer;
    FrictionParameters friction;
    SimulationControl simulation;
    ToleranceParameters tolerances;
};

class ConfigLoader {
public:
    static EngineParameters loadEngineParameters(const std::string& filename);
    static SimulationParameters loadSimulationParameters(const std::string& filename);
    
private:
    static ValveParameters loadValveParameters(const YAML::Node& node);
};

#endif
