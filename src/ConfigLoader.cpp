#include "ConfigLoader.h"
#include <stdexcept>
#include <cmath>

namespace {
double finite_value(const YAML::Node& section, const char* key) {
    const double value = section[key].as<double>();
    if (!std::isfinite(value)) {
        throw std::runtime_error(std::string(key) + " must be finite");
    }
    return value;
}

double positive_value(const YAML::Node& section, const char* key) {
    const double value = finite_value(section, key);
    if (value <= 0.0) {
        throw std::runtime_error(std::string(key) + " must be positive");
    }
    return value;
}

double nonnegative_value(const YAML::Node& section, const char* key) {
    const double value = finite_value(section, key);
    if (value < 0.0) {
        throw std::runtime_error(std::string(key) + " must be nonnegative");
    }
    return value;
}
}

EngineParameters ConfigLoader::loadEngineParameters(const std::string& filename) {
    try {
        YAML::Node config = YAML::LoadFile(filename);
        EngineParameters params;
        
        // Load geometry
        auto geometry = config["geometry"];
        params.geometry.bore = positive_value(geometry, "bore");
        params.geometry.stroke = positive_value(geometry, "stroke");
        params.geometry.connecting_rod_length = positive_value(geometry, "connecting_rod_length");
        params.geometry.compression_ratio = positive_value(geometry, "compression_ratio");
        
        if (params.geometry.compression_ratio <= 1.0) {
            throw std::runtime_error("compression_ratio must be greater than one");
        }
        if (params.geometry.connecting_rod_length <= 0.5 * params.geometry.stroke) {
            throw std::runtime_error("connecting_rod_length must exceed half the stroke");
        }

        // Load operating conditions
        auto operating = config["operating_conditions"];
        params.operating_conditions.engine_speed = positive_value(operating, "engine_speed");
        params.operating_conditions.intake_pressure = positive_value(operating, "intake_pressure");
        params.operating_conditions.intake_temperature = positive_value(operating, "intake_temperature");
        params.operating_conditions.exhaust_pressure = positive_value(operating, "exhaust_pressure");
        params.operating_conditions.intake_composition = operating["intake_composition"].as<std::string>();
        params.operating_conditions.wall_temperature = positive_value(operating, "wall_temperature");
        
        // Load spark parameters
        auto spark = config["spark"];
        params.spark.timing_btdc = finite_value(spark, "timing_btdc");
        params.spark.energy = nonnegative_value(spark, "energy");
        params.spark.duration = nonnegative_value(spark, "duration");
        params.spark.enabled = spark["enabled"].as<bool>();
        
        if (params.spark.enabled && params.spark.energy > 0.0 && params.spark.duration == 0.0) {
            throw std::runtime_error("spark duration must be positive for an active spark");
        }

        // Load valve parameters
        auto valves = config["valves"];
        params.valves.intake = loadValveParameters(valves["intake"]);
        params.valves.exhaust = loadValveParameters(valves["exhaust"]);
        
        // Load mechanism file
        params.mechanism_file = config["mechanism_file"].as<std::string>();
        
        return params;
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to load engine parameters: " + std::string(e.what()));
    }
}

SimulationParameters ConfigLoader::loadSimulationParameters(const std::string& filename) {
    try {
        YAML::Node config = YAML::LoadFile(filename);
        SimulationParameters params;
        
        // Load time stepping parameters
        auto time_stepping = config["time_stepping"];
        params.time_stepping.crank_angle_step = positive_value(time_stepping, "crank_angle_step");
        params.time_stepping.min_angle_step = positive_value(time_stepping, "min_angle_step");
        params.time_stepping.max_angle_step = positive_value(time_stepping, "max_angle_step");
        params.time_stepping.adaptive_stepping = time_stepping["adaptive_stepping"].as<bool>();
        
        if (params.time_stepping.min_angle_step > params.time_stepping.max_angle_step) {
            throw std::runtime_error("min_angle_step must not exceed max_angle_step");
        }

        // Load heat transfer parameters
        auto heat_transfer = config["heat_transfer"];
        params.heat_transfer.woschni_c1 = nonnegative_value(heat_transfer, "woschni_c1");
        params.heat_transfer.woschni_c2 = nonnegative_value(heat_transfer, "woschni_c2");
        params.heat_transfer.woschni_c3 = nonnegative_value(heat_transfer, "woschni_c3");
        
        // Load friction parameters
        auto friction = config["friction"];
        params.friction.cf_c1 = finite_value(friction, "cf_c1");
        params.friction.cf_c2 = finite_value(friction, "cf_c2");
        params.friction.cf_c3 = finite_value(friction, "cf_c3");
        params.friction.cf_c4 = finite_value(friction, "cf_c4");
        
        // Load simulation control
        auto simulation = config["simulation"];
        params.simulation.output_interval = positive_value(simulation, "output_interval");
        params.simulation.start_angle = finite_value(simulation, "start_angle");
        params.simulation.end_angle = finite_value(simulation, "end_angle");
        params.simulation.equilibrium_expansion =
            simulation["equilibrium_expansion"].as<bool>(true);
        
        if (!(params.simulation.end_angle > params.simulation.start_angle) ||
            !std::isfinite(params.simulation.end_angle - params.simulation.start_angle)) {
            throw std::runtime_error("end_angle must define a finite positive range after start_angle");
        }

        // Load tolerances
        auto tolerances = config["tolerances"];
        params.tolerances.temperature_min = positive_value(tolerances, "temperature_min");
        params.tolerances.temperature_max = positive_value(tolerances, "temperature_max");
        params.tolerances.pressure_min = positive_value(tolerances, "pressure_min");
        params.tolerances.pressure_max = positive_value(tolerances, "pressure_max");
        
        return params;
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to load simulation parameters: " + std::string(e.what()));
    }
}

ValveParameters ConfigLoader::loadValveParameters(const YAML::Node& node) {
    ValveParameters params;
    params.open_angle = finite_value(node, "open_angle");
    params.close_angle = finite_value(node, "close_angle");
    params.max_lift = nonnegative_value(node, "max_lift");
    params.diameter = positive_value(node, "diameter");
    params.seat_diameter = positive_value(node, "seat_diameter");
    params.discharge_coefficient = nonnegative_value(node, "discharge_coefficient");
    const double duration = params.close_angle - params.open_angle;
    if (duration < 0.0 || duration > 720.0 || (duration == 0.0 && params.max_lift > 0.0)) {
        throw std::runtime_error("valve close_angle must follow open_angle by up to 720 degrees");
    }
    return params;
}
