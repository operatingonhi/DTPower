#include "ConfigLoader.h"
#include <stdexcept>
#include <iostream>

EngineParameters ConfigLoader::loadEngineParameters(const std::string& filename) {
    try {
        YAML::Node config = YAML::LoadFile(filename);
        EngineParameters params;
        
        // Load geometry
        auto geometry = config["geometry"];
        params.geometry.bore = geometry["bore"].as<double>();
        params.geometry.stroke = geometry["stroke"].as<double>();
        params.geometry.connecting_rod_length = geometry["connecting_rod_length"].as<double>();
        params.geometry.compression_ratio = geometry["compression_ratio"].as<double>();
        
        // Load operating conditions
        auto operating = config["operating_conditions"];
        params.operating_conditions.engine_speed = operating["engine_speed"].as<double>();
        params.operating_conditions.intake_pressure = operating["intake_pressure"].as<double>();
        params.operating_conditions.intake_temperature = operating["intake_temperature"].as<double>();
        params.operating_conditions.exhaust_pressure = operating["exhaust_pressure"].as<double>();
        params.operating_conditions.intake_composition = operating["intake_composition"].as<std::string>();
        params.operating_conditions.wall_temperature = operating["wall_temperature"].as<double>();
        
        // Load spark parameters
        auto spark = config["spark"];
        params.spark.timing_btdc = spark["timing_btdc"].as<double>();
        params.spark.energy = spark["energy"].as<double>();
        params.spark.duration = spark["duration"].as<double>();
        params.spark.enabled = spark["enabled"].as<bool>();
        
        // Load valve parameters
        auto valves = config["valves"];
        params.valves.intake = loadValveParameters(valves["intake"]);
        params.valves.exhaust = loadValveParameters(valves["exhaust"]);
        
        // Load mechanism file
        params.mechanism_file = config["mechanism_file"].as<std::string>();
        
        return params;
        
    } catch (const YAML::Exception& e) {
        throw std::runtime_error("Failed to load engine parameters: " + std::string(e.what()));
    }
}

SimulationParameters ConfigLoader::loadSimulationParameters(const std::string& filename) {
    try {
        YAML::Node config = YAML::LoadFile(filename);
        SimulationParameters params;
        
        // Load time stepping parameters
        auto time_stepping = config["time_stepping"];
        params.time_stepping.crank_angle_step = time_stepping["crank_angle_step"].as<double>();
        params.time_stepping.min_angle_step = time_stepping["min_angle_step"].as<double>();
        params.time_stepping.max_angle_step = time_stepping["max_angle_step"].as<double>();
        params.time_stepping.adaptive_stepping = time_stepping["adaptive_stepping"].as<bool>();
        
        // Load heat transfer parameters
        auto heat_transfer = config["heat_transfer"];
        params.heat_transfer.woschni_c1 = heat_transfer["woschni_c1"].as<double>();
        params.heat_transfer.woschni_c2 = heat_transfer["woschni_c2"].as<double>();
        params.heat_transfer.woschni_c3 = heat_transfer["woschni_c3"].as<double>();
        
        // Load friction parameters
        auto friction = config["friction"];
        params.friction.cf_c1 = friction["cf_c1"].as<double>();
        params.friction.cf_c2 = friction["cf_c2"].as<double>();
        params.friction.cf_c3 = friction["cf_c3"].as<double>();
        params.friction.cf_c4 = friction["cf_c4"].as<double>();
        
        // Load simulation control
        auto simulation = config["simulation"];
        params.simulation.output_interval = simulation["output_interval"].as<double>();
        params.simulation.start_angle = simulation["start_angle"].as<double>();
        params.simulation.end_angle = simulation["end_angle"].as<double>();
        
        // Load tolerances
        auto tolerances = config["tolerances"];
        params.tolerances.temperature_min = tolerances["temperature_min"].as<double>();
        params.tolerances.temperature_max = tolerances["temperature_max"].as<double>();
        params.tolerances.pressure_min = tolerances["pressure_min"].as<double>();
        params.tolerances.pressure_max = tolerances["pressure_max"].as<double>();
        
        return params;
        
    } catch (const YAML::Exception& e) {
        throw std::runtime_error("Failed to load simulation parameters: " + std::string(e.what()));
    }
}

ValveParameters ConfigLoader::loadValveParameters(const YAML::Node& node) {
    ValveParameters params;
    params.open_angle = node["open_angle"].as<double>();
    params.close_angle = node["close_angle"].as<double>();
    params.max_lift = node["max_lift"].as<double>();
    params.diameter = node["diameter"].as<double>();
    params.seat_diameter = node["seat_diameter"].as<double>();
    params.discharge_coefficient = node["discharge_coefficient"].as<double>();
    return params;
}
