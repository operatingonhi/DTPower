#include "Engine.h"
#include "ConfigLoader.h"
#include <iostream>
#include <vector>

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
        
        // Define valve profiles using loaded parameters
        std::vector<double> intake_angles, intake_lifts;
        std::vector<double> exhaust_angles, exhaust_lifts;
        
        // Intake valve profile
        for (double angle = 0; angle <= 720; angle += 1.0) {
            intake_angles.push_back(angle);
            double intake_start = engine_params.valves.intake.open_angle;
            double intake_end = engine_params.valves.intake.close_angle;
            double max_lift = engine_params.valves.intake.max_lift;
            
            // Make valve timing periodic by checking multiple cycles
            bool valve_open = false;
            double lift_angle = 0.0;
            double duration = intake_end - intake_start;
            
            // Check current cycle (0-720°)
            if (angle >= intake_start && angle <= intake_end) {
                valve_open = true;
                lift_angle = angle - intake_start;
            }
            // Check previous cycle (-720 to 0°)
            else if (angle >= (intake_start - 720) && angle <= (intake_end - 720)) {
                valve_open = true;
                lift_angle = angle - (intake_start - 720);
            }
            // Check next cycle (720-1440°)
            else if (angle >= (intake_start + 720) && angle <= (intake_end + 720)) {
                valve_open = true;
                lift_angle = angle - (intake_start + 720);
            }
            
            if (valve_open) {
                if (lift_angle <= duration/2.0) {
                    intake_lifts.push_back(max_lift * lift_angle / (duration/2.0));
                } else {
                    intake_lifts.push_back(max_lift * (duration - lift_angle) / (duration/2.0));
                }
            } else {
                intake_lifts.push_back(0.0);
            }
        }
        
        // Exhaust valve profile
        for (double angle = 0; angle <= 720; angle += 1.0) {
            exhaust_angles.push_back(angle);
            double exhaust_start = engine_params.valves.exhaust.open_angle;
            double exhaust_end = engine_params.valves.exhaust.close_angle;
            double max_lift = engine_params.valves.exhaust.max_lift;
            
            // Make valve timing periodic by checking multiple cycles
            bool valve_open = false;
            double lift_angle = 0.0;
            double duration = exhaust_end - exhaust_start;
            
            // Check current cycle (0-720°)
            if (angle >= exhaust_start && angle <= exhaust_end) {
                valve_open = true;
                lift_angle = angle - exhaust_start;
            }
            // Check previous cycle (-720 to 0°)
            else if (angle >= (exhaust_start - 720) && angle <= (exhaust_end - 720)) {
                valve_open = true;
                lift_angle = angle - (exhaust_start - 720);
            }
            // Check next cycle (720-1440°)
            else if (angle >= (exhaust_start + 720) && angle <= (exhaust_end + 720)) {
                valve_open = true;
                lift_angle = angle - (exhaust_start + 720);
            }
            
            if (valve_open) {
                if (lift_angle <= duration/2.0) {
                    exhaust_lifts.push_back(max_lift * lift_angle / (duration/2.0));
                } else {
                    exhaust_lifts.push_back(max_lift * (duration - lift_angle) / (duration/2.0));
                }
            } else {
                exhaust_lifts.push_back(0.0);
            }
        }
        
        engine.set_valve_profiles(intake_angles, intake_lifts, exhaust_angles, exhaust_lifts);
        
        std::cout << "Running engine cycle simulation..." << std::endl;
        
        // Run simulation
        engine.run_cycle();
        
        // Output results
        engine.output_results("engine_results.csv");
        
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
