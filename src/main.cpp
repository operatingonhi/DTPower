#include "Engine.h"
#include "ConfigLoader.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
std::vector<double> valve_lifts(const std::vector<double>& angles, const ValveParameters& valve) {
    std::vector<double> lifts(angles.size(), 0.0);
    if (valve.max_lift == 0.0) return lifts;

    const double duration = valve.close_angle - valve.open_angle;
    for (std::size_t i = 0; i < angles.size(); ++i) {
        double elapsed = std::fmod(angles[i] - valve.open_angle, 720.0);
        if (elapsed < 0.0) elapsed += 720.0;
        if (elapsed < duration) {
            lifts[i] = valve.max_lift * (std::min(elapsed, duration - elapsed) / (0.5 * duration));
        }
    }
    return lifts;
}
}

int main() {
    try {
        std::cout << "Loading configuration files..." << std::endl;
        
        // Load configuration parameters
        EngineParameters engine_params = ConfigLoader::loadEngineParameters("engine_parameters.yaml");
        SimulationParameters sim_params = ConfigLoader::loadSimulationParameters("simulation_parameters.yaml");
        
        std::cout << "Initializing engine simulation..." << std::endl;
        
        // Create engine with loaded parameters
        Engine engine(engine_params.geometry.bore, 
                     engine_params.geometry.stroke, 
                     engine_params.geometry.connecting_rod_length,
                     engine_params.geometry.compression_ratio, 
                     engine_params.mechanism_file,
                     sim_params,
                     engine_params);
        
        // Set operating conditions
        engine.set_engine_speed(engine_params.operating_conditions.engine_speed);
        engine.set_boundary_conditions(engine_params.operating_conditions.intake_pressure,
                                     engine_params.operating_conditions.intake_temperature,
                                     engine_params.operating_conditions.exhaust_pressure);
        engine.set_intake_composition(engine_params.operating_conditions.intake_composition);
        
        // Configure spark ignition
        engine.set_spark_timing(engine_params.spark.timing_btdc);
        engine.set_spark_energy(engine_params.spark.energy);
        engine.enable_spark(engine_params.spark.enabled);
        
        std::vector<double> valve_angles;
        for (int angle = 0; angle <= 720; ++angle) {
            valve_angles.push_back(angle);
        }
        engine.set_valve_profiles(valve_angles, valve_lifts(valve_angles, engine_params.valves.intake),
                                  valve_angles, valve_lifts(valve_angles, engine_params.valves.exhaust));

        std::cout << "Running engine cycle simulation..." << std::endl;
        
        // Run simulation
        engine.run_cycle();
        
        // Output performance summary
        engine.output_performance_summary("engine_performance.csv");
        
        // Display performance metrics
        std::cout << "\nEngine Performance Results:" << std::endl;
        std::cout << "IMEP: " << engine.get_indicated_mean_effective_pressure() / 1.0e5 << " bar" << std::endl;
        std::cout << "GMEP: " << engine.get_gross_mean_effective_pressure() / 1.0e5 << " bar" << std::endl;
        std::cout << "PMEP: " << engine.get_pumping_mean_effective_pressure() / 1.0e5 << " bar" << std::endl;
        std::cout << "BMEP: " << engine.get_brake_mean_effective_pressure() / 1.0e5 << " bar" << std::endl;
        std::cout << "Indicated Power: " << engine.get_indicated_power() / 1000.0 << " kW" << std::endl;
        std::cout << "Brake Power: " << engine.get_brake_power() / 1000.0 << " kW" << std::endl;
        std::cout << "Cycle Efficiency: " << engine.get_cycle_efficiency() * 100.0 << " %" << std::endl;
        
        // Print detailed energy ledger
        engine.print_energy_ledger();
        
        std::cout << "\nSimulation completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
