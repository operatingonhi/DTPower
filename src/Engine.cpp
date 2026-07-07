#include "Engine.h"
#include "Profiler.h"
#include "DebugLogger.h"
#include <cmath>
#include <iostream>
#include <chrono>


Engine::Engine(double cylinder_bore, double cylinder_stroke, double connecting_rod_length,
               double compression_ratio, const std::string& mechanism_file,
               const SimulationParameters& sim_params, const EngineParameters& engine_params)
    : bore(cylinder_bore), stroke(cylinder_stroke), con_rod_length(connecting_rod_length),
      compression_ratio(compression_ratio), engine_speed(1000.0), crank_angle(0.0),
      crank_angle_step(sim_params.time_stepping.crank_angle_step), time_step(0.0), 
      intake_pressure(101325.0), intake_temperature(300.0), exhaust_pressure(101325.0),
      intake_composition("O2:0.21, N2:0.79"), mechanism_file(mechanism_file),
      fuel_energy_rate(0.0), lower_heating_value(0.0),
      sim_params(sim_params), engine_params(engine_params),
      total_lhv_energy(0.0), indicated_work(0.0), pumping_work(0.0), 
      friction_work(0.0), heat_transfer_loss(0.0), exhaust_energy_loss(0.0),
      reference_enthalpy(0.0), T3(0.0), valves_flowing(false) {
    
    // Calculate clearance volume from compression ratio
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    clearance_volume = displacement / (compression_ratio - 1.0);
    
    // Initialize components with parameters from YAML files
    double T_wall = engine_params.operating_conditions.wall_temperature;
    cylinder = std::make_unique<Cylinder>(bore, stroke, clearance_volume, mechanism_file, T_wall);
    
    // Update heat transfer model with Woschni correlation parameters
    cylinder->update_heat_transfer_parameters(sim_params.heat_transfer.woschni_c1,
                                             sim_params.heat_transfer.woschni_c2,
                                             sim_params.heat_transfer.woschni_c3);
    
    intake_valve = std::make_unique<ValveFlow>(engine_params.valves.intake.diameter,
                                              engine_params.valves.intake.seat_diameter,
                                              engine_params.valves.intake.discharge_coefficient);
    exhaust_valve = std::make_unique<ValveFlow>(engine_params.valves.exhaust.diameter,
                                               engine_params.valves.exhaust.seat_diameter,
                                               engine_params.valves.exhaust.discharge_coefficient);
    
    friction = std::make_unique<FrictionModel>(bore, stroke, con_rod_length,
                                              sim_params.friction.cf_c1,
                                              sim_params.friction.cf_c2,
                                              sim_params.friction.cf_c3,
                                              sim_params.friction.cf_c4);
    
    crankshaft = std::make_unique<Crankshaft>(stroke, con_rod_length, bore);
    
    spark = std::make_unique<SparkModel>(engine_params.spark.timing_btdc,
                                        engine_params.spark.duration,
                                        engine_params.spark.energy);
    
    // Initialize persistent intake gas object
    intake_gas = std::make_unique<CanteraInterface>(mechanism_file);
    
    // Enable combustion logging
    cylinder->enable_combustion_logging("combustion_diagnostics.csv");
    
    calculate_time_step();
}

void Engine::set_engine_speed(double rpm) {
    engine_speed = rpm;
    calculate_time_step();
}

void Engine::set_boundary_conditions(double P_intake, double T_intake, double P_exhaust) {
    intake_pressure = P_intake;
    intake_temperature = T_intake;
    exhaust_pressure = P_exhaust;
    
    // Update the persistent intake gas object with new conditions
    if (intake_gas) {
        intake_gas->set_state_TPX(intake_temperature, intake_pressure, intake_composition);
    }
}

void Engine::set_intake_composition(const std::string& composition) {
    intake_composition = composition;
    
    // Update the persistent intake gas object with new composition
    if (intake_gas) {
        intake_gas->set_state_TPX(intake_temperature, intake_pressure, intake_composition);
    }
}

void Engine::set_valve_profiles(const std::vector<double>& intake_angles,
                               const std::vector<double>& intake_lifts,
                               const std::vector<double>& exhaust_angles,
                               const std::vector<double>& exhaust_lifts) {
    intake_valve->set_lift_profile(intake_angles, intake_lifts);
    exhaust_valve->set_lift_profile(exhaust_angles, exhaust_lifts);
}

void Engine::calculate_time_step() {
    // Convert RPM to rad/s and calculate time step from crank angle step
    double omega = engine_speed * 2.0 * M_PI / 60.0;  // [rad/s]
    time_step = crank_angle_step * M_PI / 180.0 / omega;  // [s]
}

void Engine::update_cylinder_volume() {
    double volume = crankshaft->calculate_cylinder_volume(crank_angle, clearance_volume);
    cylinder->set_volume(volume);
}

void Engine::calculate_valve_flows() {
    double intake_lift = intake_valve->get_lift(crank_angle);
    double exhaust_lift = exhaust_valve->get_lift(crank_angle);
    
    // Track if any mass flow is occurring
    bool mass_flowing = false;
    
    // Specific gas constant for air [J/(kg·K)]
    const double R_air = 287.0;
    
    // Calculate intake mass flow
    if (intake_lift > 0.0) {
        double intake_density = intake_pressure / (R_air * intake_temperature);
        double mass_flow = intake_valve->calculate_mass_flow(
            intake_lift, intake_pressure, cylinder->get_pressure(),
            intake_temperature, intake_density);
        
        if (mass_flow > 0.0) {
            // Flow into cylinder
            intake_gas->set_state_TPX(intake_temperature, intake_pressure, intake_composition);
            cylinder->add_mass_flow(mass_flow, time_step, intake_composition, intake_temperature);
            mass_flowing = true;
        } else if (mass_flow < 0.0) {
            // Backflow out of cylinder
            cylinder->remove_mass_flow(-mass_flow, time_step);
            mass_flowing = true;
        }
    }
    
    // Calculate exhaust mass flow
    if (exhaust_lift > 0.0 && cylinder->get_pressure() > exhaust_pressure) {
        double mass_flow = exhaust_valve->calculate_mass_flow(
            exhaust_lift, cylinder->get_pressure(), exhaust_pressure,
            cylinder->get_temperature(), cylinder->get_density());
        
        if (mass_flow > 0.0) {
            // Calculate exhaust energy loss before removing mass
            double exhaust_enthalpy = cylinder->get_enthalpy();  // [J/kg]
            double enthalpy_difference = exhaust_enthalpy - reference_enthalpy;  // [J/kg]
            double energy_loss_rate = mass_flow * enthalpy_difference;  // [W]
            
            // Accumulate exhaust energy loss (only positive losses)
            if (enthalpy_difference > 0.0) {
                exhaust_energy_loss += energy_loss_rate * time_step;  // [J]
            }
            
            cylinder->remove_mass_flow(mass_flow, time_step);
            mass_flowing = true;
        }
    }
    
    // Store valve flow status for use in update_thermodynamics
    valves_flowing = mass_flowing;
}

void Engine::store_results() {
    results.crank_angles.push_back(crank_angle);
    results.cylinder_pressures.push_back(cylinder->get_pressure());
    results.cylinder_temperatures.push_back(cylinder->get_temperature());
    results.cylinder_volumes.push_back(cylinder->get_volume());
    results.cylinder_masses.push_back(cylinder->get_mass());
    results.heat_release_rates.push_back(cylinder->get_heat_release_rate());
    results.heat_transfer_rates.push_back(cylinder->get_heat_transfer_rate());
    
    // Calculate and store Woschni heat transfer variables
    double T_gas = cylinder->get_temperature();
    double mean_piston_speed = crankshaft->calculate_mean_piston_speed(engine_speed);
    double surface_area = cylinder->get_heat_transfer()->calculate_surface_area(
        cylinder->get_volume(), bore, stroke, clearance_volume);
    double h_conv = cylinder->get_heat_transfer()->calculate_heat_transfer_coefficient(
        cylinder->get_pressure(), T_gas, mean_piston_speed, 0.0);
    
    results.gas_temperatures.push_back(T_gas);
    results.surface_areas.push_back(surface_area);
    results.heat_transfer_coefficients.push_back(h_conv);
    
    // Store friction power
    double piston_velocity = crankshaft->calculate_piston_velocity(crank_angle, engine_speed);
    double friction_power = friction->calculate_instantaneous_friction_power(
        engine_speed, cylinder->get_pressure(), crank_angle, piston_velocity);
    results.friction_powers.push_back(friction_power);
    
    results.gamma_values.push_back(cylinder->get_gamma());
    
    // Store species mass fractions
    std::vector<double> mass_fractions = cylinder->get_gas()->mass_fractions();
    results.species_mass_fractions.push_back(mass_fractions);
    
    // Store pumping work data (pressure difference during intake/exhaust strokes)
    double pumping_pressure_diff = 0.0;
    // Intake stroke: 0-180°, Exhaust stroke: 540-720°
    if ((crank_angle >= 0.0 && crank_angle <= 180.0) || 
        (crank_angle >= 540.0 && crank_angle <= 720.0)) {
        pumping_pressure_diff = cylinder->get_pressure() - intake_pressure;
    }
    results.pumping_works.push_back(pumping_pressure_diff);
    
    // Store current exhaust energy loss rate
    results.exhaust_energy_losses.push_back(exhaust_energy_loss);
    
    // Store valve flows (simplified)
    results.intake_mass_flows.push_back(0.0);
    results.exhaust_mass_flows.push_back(0.0);
}

void Engine::initialize_cycle() {
    crank_angle = sim_params.simulation.start_angle;
    results = EngineResults();  // Clear previous results
    
    // Reset energy tracking variables
    exhaust_energy_loss = 0.0;
    heat_transfer_loss = 0.0;
    T3 = 0.0;  // Reset T3 for new cycle
    valves_flowing = false;  // Reset valve flow status
    
    // Initialize cylinder state at starting angle
    cylinder->get_gas()->set_state_TPX(intake_temperature, intake_pressure, intake_composition);
    
    // Update cylinder volume for the starting position
    update_cylinder_volume();
    
    // Update cylinder mass based on the initialized state
    double initial_density = cylinder->get_gas()->density();
    double initial_volume = cylinder->get_volume();
    cylinder->set_mass(initial_density * initial_volume);
    
    // Initialize previous state for combustion zone tracking
    // This is done in update_thermodynamics, but we need initial values
    // The first call to update_thermodynamics will set these properly
    
    // Get species names for CSV header (only need to do this once)
    if (results.species_names.empty()) {
        results.species_names = cylinder->get_gas()->species_names();
    }
    
    // Start friction integration
    friction->start_friction_integration(crank_angle);
    
    // Calculate Lower Heating Value (LHV) for efficiency calculation
    try {
        auto lhv_gas = std::make_unique<CanteraInterface>(mechanism_file);
        
        // Set to intake conditions for LHV calculation
        lhv_gas->set_state_TPX(intake_temperature, intake_pressure, intake_composition);
        
        // Store initial enthalpy
        double h_initial = lhv_gas->enthalpy();
        
        // Perform equilibrium combustion at constant temperature and pressure
        lhv_gas->get_solution()->thermo()->equilibrate("TP");
        
        // Calculate enthalpy change (LHV per unit mass of mixture)
        double h_final = lhv_gas->enthalpy();
        lower_heating_value = h_initial - h_final; // J/kg of mixture
        
        // Calculate reference enthalpy at intake conditions for exhaust energy loss
        lhv_gas->set_state_TP(intake_temperature, intake_pressure);
        reference_enthalpy = lhv_gas->enthalpy();  // [J/kg]
        
    } catch (const std::exception& e) {
        std::cerr << "Warning: LHV calculation failed: " << e.what() << std::endl;
        lower_heating_value = 0.0;
        fuel_energy_rate = 0.0;
    }
}

void Engine::run_cycle() {
    PROFILE_FUNCTION();
    
    // Run motored cycle first to get motored pressure profile
    std::vector<double> motored_pressures;
    std::vector<double> motored_angles;
    
    std::cout << "\n=== Running motored cycle to establish baseline pressure ===" << std::endl;
    initialize_cycle();  // Initialize for motored run
    run_cycle_internal(true, motored_pressures, motored_angles);
    
    // Now run the fired cycle with motored pressure data
    std::cout << "\n=== Running fired cycle with chemistry ===" << std::endl;
    
    // Reset for fired cycle
    initialize_cycle();
    run_cycle_internal(false, motored_pressures, motored_angles);
}

void Engine::run_cycle_internal(bool is_motored, std::vector<double>& motored_pressures, std::vector<double>& motored_angles) {
    PROFILE_FUNCTION();
    
    // Use simulation parameters instead of hardcoded values
    double target_crank_angle = sim_params.simulation.end_angle;
    double adaptive_angle_step = crank_angle_step;
    double min_angle_step = sim_params.time_stepping.min_angle_step;
    double max_angle_step = sim_params.time_stepping.max_angle_step;
    bool adaptive_stepping = sim_params.time_stepping.adaptive_stepping;
    
    // Track if we've set reference values at IVC
    bool reference_values_set = false;
    const double IVC_angle = 180.0;  // Intake valve closure occurs at ~180° ATDC
    const double IVC_tolerance = 1.0;  // Capture within 1 degree of IVC
    
    // Progress tracking variables
    double start_angle = sim_params.simulation.start_angle;
    double total_angle_range = target_crank_angle - start_angle;
    int last_progress_percent = -1;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    const int progress_bar_width = 50;
    std::string cycle_type = is_motored ? "motored" : "fired";
    std::cout << "Starting " << cycle_type << " cycle simulation..." << std::endl;
    std::cout << "Progress: [";
    for (int i = 0; i < progress_bar_width; ++i) std::cout << " ";
    std::cout << "] 0%" << std::flush;
    
    // For fired cycle, we need motored pressure data for heat transfer calculations
    // For motored cycle, we're building the pressure profile
    
    for (crank_angle = start_angle; crank_angle < target_crank_angle; ) {
        PROFILE_SCOPE("main_loop_iteration");
        
        // Interpolate motored pressure for current crank angle (only used in fired cycle)
        double motored_pressure = 0.0;
        if (!is_motored && !motored_pressures.empty()) {
            motored_pressure = interpolate_motored_pressure(crank_angle, motored_angles, motored_pressures);
        }
        
        // Calculate progress percentage
        double progress = (crank_angle - start_angle) / total_angle_range;
        int progress_percent = static_cast<int>(progress * 100);
        
        // Update progress bar every 2%
        if (progress_percent > last_progress_percent && progress_percent % 2 == 0) {
            last_progress_percent = progress_percent;
            
            // Clear the line and redraw progress bar
            std::cout << "\rProgress: [";
            int filled = static_cast<int>(progress * 50);
            for (int i = 0; i < 50; ++i) {
                if (i < filled) std::cout << "=";
                else if (i == filled) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << progress_percent << "%" << std::flush;
            
            // Show current crank angle and phase
            std::string phase;
            double angle_mod = std::fmod(crank_angle, 720.0);
            if (angle_mod < 0) angle_mod += 720.0;
            
            if (angle_mod < 180) phase = "Intake";
            else if (angle_mod < 360) phase = "Compression";
            else if (angle_mod < 540) phase = "Expansion";
            else phase = "Exhaust";
            
            std::cout << " (" << static_cast<int>(crank_angle) << "_ - " << phase << ")" << std::flush;
        }
        
        // Check if we're in a reactive phase requiring smaller time steps
        const double reactive_temperature_threshold = 800.0;  // [K]
        double spark_power = 0.0;
        bool reactive_phase = false;
        
        if (!is_motored) {
            spark_power = spark->get_spark_power(crank_angle);
            reactive_phase = (cylinder->get_temperature() > reactive_temperature_threshold) || 
                            (spark_power > 0.0) ||
                            (spark->is_spark_active(crank_angle));
        }
        
        // Adjust angle step based on conditions
        if (adaptive_stepping) {
            const double step_reduction_factor = 0.5;
            const double step_increase_factor = 1.1;
            
            if (reactive_phase) {
                adaptive_angle_step = std::max(min_angle_step, adaptive_angle_step * step_reduction_factor);
            } else {
                adaptive_angle_step = std::min(max_angle_step, adaptive_angle_step * step_increase_factor);
            }
        }
        
        // Ensure we don't overshoot target angle
        adaptive_angle_step = std::min(adaptive_angle_step, target_crank_angle - crank_angle);
        
        // Calculate corresponding time step
        double omega = engine_speed * 2.0 * M_PI / 60.0;  // [rad/s]
        double adaptive_time_step = adaptive_angle_step * M_PI / 180.0 / omega;  // [s]
        time_step = adaptive_time_step;
        
        // Update cylinder volume due to piston motion
        {
            PROFILE_SCOPE("update_cylinder_volume");
            update_cylinder_volume();
        }
        
        // Calculate and apply mass flows before thermodynamics
        {
            PROFILE_SCOPE("calculate_valve_flows");
            calculate_valve_flows();
        }
        
        // Handle spark ignition (only for fired cycle)
        if (!is_motored && spark_power > 0.0) {
            PROFILE_SCOPE("spark_ignition");
            // Calculate fuel energy input rate (chemical energy entering cylinder)
            const double strokes_per_cycle = 2.0;  // 4-stroke engine
            double cycles_per_second = engine_speed / (60.0 * strokes_per_cycle);
            double mass_flow_rate = cylinder->get_mass() * cycles_per_second;
            fuel_energy_rate = mass_flow_rate * lower_heating_value;  // [W]
            cylinder->initiate_spark_ignition(spark_power * adaptive_time_step);
        }
        
        // Capture peak temperature during combustion (only for fired cycle)
        if (!is_motored && cylinder->get_temperature() > T3) {
            T3 = cylinder->get_temperature();
        }
        
        // Update thermodynamics
        {
            PROFILE_SCOPE("update_thermodynamics");
            // Only advance chemistry for fired cycle
            if (!is_motored && reactive_temperature_threshold < cylinder->get_temperature()){
                cylinder->advance_chemistry(adaptive_time_step);
            }
            double mean_piston_speed = crankshaft->calculate_mean_piston_speed(engine_speed);
            // For motored run, use 0.0 for motored_pressure (no combustion term)
            // For fired run, use interpolated motored_pressure
            cylinder->heat_loss(adaptive_time_step, mean_piston_speed, motored_pressure, valves_flowing);
        }
        
        // Calculate and accumulate friction work (only for fired cycle)
        if (!is_motored) {
            PROFILE_SCOPE("friction_calculation");
            double piston_velocity = crankshaft->calculate_piston_velocity(crank_angle, engine_speed);
            double friction_power = friction->calculate_instantaneous_friction_power(
                engine_speed, cylinder->get_pressure(), crank_angle, piston_velocity);
            friction->accumulate_friction_work(friction_power, adaptive_time_step);
        }
        
        // Track energy flows for ledger (only for fired cycle)
        if (!is_motored) {
            PROFILE_SCOPE("energy_tracking");
            double heat_loss = cylinder->get_heat_transfer_rate() * adaptive_time_step;
            heat_transfer_loss += std::abs(heat_loss);  // Always positive for losses
        }
        
        crank_angle += adaptive_angle_step;
        
        // Capture reference values at IVC during motored run
        if (is_motored && !reference_values_set && 
            crank_angle >= (IVC_angle - IVC_tolerance) && 
            crank_angle <= (IVC_angle + IVC_tolerance)) {
            
            // Get state at IVC
            double T_ivc = cylinder->get_temperature();  // [K]
            double P_ivc = cylinder->get_pressure();     // [Pa]
            double V_ivc = cylinder->get_volume();       // [m³]
            
            // Set reference values for Woschni heat transfer
            cylinder->set_heat_transfer_reference_values(T_ivc, P_ivc, V_ivc);
            reference_values_set = true;
            
            std::cout << "\nReference values set at IVC (CA=" << crank_angle << "°):" << std::endl;
            std::cout << "  T_ivc = " << T_ivc << " K" << std::endl;
            std::cout << "  P_ivc = " << P_ivc / 1.0e3 << " kPa" << std::endl;
            std::cout << "  V_ivc = " << V_ivc * 1.0e6 << " cm³" << std::endl;
        }
        
        // Store results at specified intervals (starting from 0°)
        if (crank_angle >= 0.0 && 
            (results.crank_angles.empty() || 
             crank_angle - results.crank_angles.back() >= sim_params.simulation.output_interval)) {
            
            if (is_motored) {
                // Store motored pressure profile for use in fired cycle
                motored_pressures.push_back(cylinder->get_pressure());
                motored_angles.push_back(crank_angle);
            }
            
            store_results();
        }
    }
    
    // Complete progress bar and show timing
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "\rProgress: [";
    for (int i = 0; i < progress_bar_width; ++i) std::cout << "=";
    std::cout << "] 100% - Complete!" << std::endl;
    std::cout << "Simulation completed in " << duration.count() << " ms" << std::endl;
    
    // Calculate energy ledger and output results based on cycle type
    if (is_motored) {
        // Output motored results to separate file
        output_results("motored_results.csv");
    } else {
        // Fired cycle: calculate energy ledger and output main results
        calculate_energy_ledger();
        output_results("engine_results.csv");
        
        #ifdef PROFILING_ENABLED
        Profiler::getInstance().printReport();
        #endif
        
        // Disable combustion logging at end of cycle
        cylinder->disable_combustion_logging();
    }
}

double Engine::interpolate_motored_pressure(double angle, const std::vector<double>& angles, const std::vector<double>& pressures) {
    if (angles.empty() || pressures.empty()) {
        return 0.0;  // No motored data available
    }
    
    // Handle angles outside the motored data range by using modulo arithmetic
    // This allows the simulation to run beyond 720 degrees
    double angle_range = angles.back() - angles.front();
    double normalized_angle = angle;
    
    // Normalize angle to be within the motored data range
    while (normalized_angle > angles.back()) {
        normalized_angle -= 720.0;  // Subtract one complete cycle
    }
    while (normalized_angle < angles.front()) {
        normalized_angle += 720.0;  // Add one complete cycle
    }
    
    // Find the two nearest points for interpolation
    for (size_t i = 0; i < angles.size() - 1; ++i) {
        if (normalized_angle >= angles[i] && normalized_angle <= angles[i + 1]) {
            // Linear interpolation
            double t = (normalized_angle - angles[i]) / (angles[i + 1] - angles[i]);
            return pressures[i] + t * (pressures[i + 1] - pressures[i]);
        }
    }
    
    // If we get here, return the closest value
    if (normalized_angle < angles.front()) {
        return pressures.front();
    } else {
        return pressures.back();
    }
}

void Engine::run_multiple_cycles(int num_cycles) {
    auto overall_start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_cycles; ++i) {
        std::cout << "\n=== Running cycle " << (i + 1) << " of " << num_cycles << " ===" << std::endl;
        run_cycle();
        
        // Show estimated time remaining
        if (i < num_cycles - 1) {
            auto cycle_end = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(cycle_end - overall_start);
            double avg_time_per_cycle = elapsed.count() / static_cast<double>(i + 1);
            double estimated_remaining = avg_time_per_cycle * (num_cycles - i - 1);
            std::cout << "Estimated time remaining: " << static_cast<int>(estimated_remaining) << " s" << std::endl;
        }
    }
    
    auto overall_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(overall_end - overall_start);
    std::cout << "\nAll " << num_cycles << " cycles completed in " << total_duration.count() << " s" << std::endl;
}

void Engine::output_results(const std::string& filename) {
    std::ofstream file(filename);
    
    // Write header
    file << "CrankAngle[deg],Pressure[Pa],Temperature[K],Volume[m3],Mass[kg],HeatReleaseRate[W],HeatTransferRate[W],h_conv[W/(m2K)],SurfaceArea[m2],T_gas[K],FrictionPower[W],Gamma[-],PumpingWork[Pa],ExhaustEnergyLoss[J]";
    
    // Add species mass fraction columns
    for (const auto& species : results.species_names) {
        // Replace spaces and commas with underscores to avoid CSV parsing issues
        std::string clean_species = species;
        std::replace(clean_species.begin(), clean_species.end(), ' ', '_');
        std::replace(clean_species.begin(), clean_species.end(), ',', '_');
        file << ",Y_" << clean_species << "[-]";
    }
    file << "\n";
    
    // Write data
    for (size_t i = 0; i < results.crank_angles.size(); ++i) {
        file << results.crank_angles[i] << ","
             << results.cylinder_pressures[i] << ","
             << results.cylinder_temperatures[i] << ","
             << results.cylinder_volumes[i] << ","
             << results.cylinder_masses[i] << ","
             << results.heat_release_rates[i] << ","
             << results.heat_transfer_rates[i] << ","
             << (i < results.heat_transfer_coefficients.size() ? results.heat_transfer_coefficients[i] : 0.0) << ","
             << (i < results.surface_areas.size() ? results.surface_areas[i] : 0.0) << ","
             << (i < results.gas_temperatures.size() ? results.gas_temperatures[i] : 0.0) << ","
             << results.friction_powers[i] << ","
             << results.gamma_values[i] << ","
             << (i < results.pumping_works.size() ? results.pumping_works[i] : 0.0) << ","
             << (i < results.exhaust_energy_losses.size() ? results.exhaust_energy_losses[i] : 0.0);
        
        // Write species mass fractions
        if (i < results.species_mass_fractions.size()) {
            const auto& mass_fractions = results.species_mass_fractions[i];
            for (size_t j = 0; j < results.species_names.size() && j < mass_fractions.size(); ++j) {
                file << "," << mass_fractions[j];
            }
            // Fill remaining columns with zeros if we have fewer species than expected
            for (size_t j = mass_fractions.size(); j < results.species_names.size(); ++j) {
                file << ",0.0";
            }
        } else {
            // Fill all species columns with zeros if no mass fraction data
            for (size_t j = 0; j < results.species_names.size(); ++j) {
                file << ",0.0";
            }
        }
        file << "\n";
    }
    
    file.close();
    std::cout << "Results written to " << filename << std::endl;
}

double Engine::get_indicated_mean_effective_pressure() const {
    // Calculate IMEP from pressure-volume data (entire cycle)
    double indicated_work = 0.0;
    for (size_t i = 1; i < results.cylinder_pressures.size(); ++i) {
        double dV = results.cylinder_volumes[i] - results.cylinder_volumes[i - 1];
        indicated_work += results.cylinder_pressures[i] * dV;
    }
    
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    return indicated_work / displacement;
}

double Engine::get_brake_mean_effective_pressure() const {
    double imep = get_indicated_mean_effective_pressure();
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    double fmep = friction->get_integrated_fmep(displacement);
    return imep - fmep;
}

double Engine::get_gross_mean_effective_pressure() const {
    // Calculate GMEP from compression and expansion strokes only (180-540°)
    double gross_work = 0.0;
    for (size_t i = 1; i < results.cylinder_pressures.size(); ++i) {
        double angle = results.crank_angles[i];
        if (angle >= 180.0 && angle <= 540.0) {
            double dV = results.cylinder_volumes[i] - results.cylinder_volumes[i - 1];
            gross_work += results.cylinder_pressures[i] * dV;
        }
    }
    
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    return gross_work / displacement;
}

double Engine::get_pumping_mean_effective_pressure() const {
    // Calculate PMEP from intake and exhaust strokes (0-180° and 540-720°)
    double pumping_work_total = 0.0;
    for (size_t i = 1; i < results.cylinder_pressures.size(); ++i) {
        double angle = results.crank_angles[i];
        if ((angle >= 0.0 && angle <= 180.0) || (angle >= 540.0 && angle <= 720.0)) {
            double dV = results.cylinder_volumes[i] - results.cylinder_volumes[i - 1];
            pumping_work_total += results.cylinder_pressures[i] * dV;
        }
    }
    
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    return pumping_work_total / displacement;
}

double Engine::get_pumping_work() const {
    double pmep = get_pumping_mean_effective_pressure();
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    return std::abs(pmep * displacement);  // [J] - absolute value for energy ledger
}

double Engine::get_indicated_power() const {
    double imep = get_indicated_mean_effective_pressure();
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    const double strokes_per_cycle = 2.0;  // 4-stroke engine
    return imep * displacement * engine_speed / (60.0 * strokes_per_cycle);  // [W]
}

double Engine::get_brake_power() const {
    double bmep = get_brake_mean_effective_pressure();
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    const double strokes_per_cycle = 2.0;  // 4-stroke engine
    return bmep * displacement * engine_speed / (60.0 * strokes_per_cycle);  // [W]
}

void Engine::set_spark_timing(double timing_btdc) {
    spark->set_spark_timing(timing_btdc);
}

void Engine::set_spark_energy(double energy) {
    spark->set_spark_energy(energy);
}

void Engine::enable_spark(bool enabled) {
    spark->enable_spark(enabled);
}

double Engine::get_cycle_efficiency() const {
    if (fuel_energy_rate <= 0.0) {
        return 0.0;
    }
    return get_brake_power() / fuel_energy_rate;
}

void Engine::calculate_energy_ledger() {
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    const double strokes_per_cycle = 2.0;  // 4-stroke engine
    double cycles_per_second = engine_speed / (60.0 * strokes_per_cycle);
    double cycle_time = 1.0 / cycles_per_second;
    
    total_lhv_energy = fuel_energy_rate * cycle_time;  // [J] per cycle
    indicated_work = get_indicated_mean_effective_pressure() * displacement;  // [J]
    pumping_work = get_pumping_work();  // [J]
    friction_work = friction->get_integrated_fmep(displacement) * displacement;  // [J]
    
    // exhaust_energy_loss is accumulated during the cycle
}

void Engine::print_energy_ledger() const {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "                    ENERGY LEDGER" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // Convert to kJ for better readability
    double total_kJ = total_lhv_energy / 1000.0;
    double indicated_kJ = indicated_work / 1000.0;
    double pumping_kJ = pumping_work / 1000.0;
    double friction_kJ = friction_work / 1000.0;
    double heat_transfer_kJ = heat_transfer_loss / 1000.0;
    double exhaust_kJ = exhaust_energy_loss / 1000.0;
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::string(60, '-') << std::endl;
    std::cout << "Energy Available (LHV):           " << std::setw(8) << total_kJ 
              << " kJ  (100.0%)" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    if (total_lhv_energy > 0) {
        double indicated_pct = (indicated_work / total_lhv_energy) * 100.0;
        double pumping_pct = (pumping_work / total_lhv_energy) * 100.0;
        double friction_pct = (friction_work / total_lhv_energy) * 100.0;
        double heat_transfer_pct = (heat_transfer_loss / total_lhv_energy) * 100.0;
        double exhaust_pct = (exhaust_energy_loss / total_lhv_energy) * 100.0;
        
        std::cout << "Indicated Work (Useful):          " << std::setw(8) << indicated_kJ 
                  << " kJ  (" << std::setw(5) << indicated_pct << "%)" << std::endl;
        std::cout << "Friction Losses:                  " << std::setw(8) << friction_kJ 
                  << " kJ  (" << std::setw(5) << friction_pct << "%)" << std::endl;
        std::cout << "Pumping Losses:                   " << std::setw(8) << pumping_kJ 
                  << " kJ  (" << std::setw(5) << pumping_pct << "%)" << std::endl;
        std::cout << "Heat Transfer Losses:             " << std::setw(8) << heat_transfer_kJ 
                  << " kJ  (" << std::setw(5) << heat_transfer_pct << "%)" << std::endl;
        std::cout << "Exhaust Energy Losses:            " << std::setw(8) << exhaust_kJ 
                  << " kJ  (" << std::setw(5) << exhaust_pct << "%)" << std::endl;
        
        std::cout << std::string(60, '-') << std::endl;
        
        // Indicated work already includes friction and pumping work
        double total_accounted = indicated_kJ + heat_transfer_kJ + exhaust_kJ;
        double balance_error = total_kJ - total_accounted;
        
        std::cout << "Total Accounted:                  " << std::setw(8) << total_accounted 
                  << " kJ" << std::endl;
        std::cout << "Energy Balance Error:             " << std::setw(8) << balance_error 
                  << " kJ" << std::endl;
        
        std::cout << std::string(60, '=') << std::endl;
        std::cout << "Net Brake Work:                   " << std::setw(8) 
                  << (indicated_kJ - friction_kJ) << " kJ" << std::endl;
        std::cout << "Brake Efficiency:                 " << std::setw(7) 
                  << get_cycle_efficiency() * 100.0 << " %" << std::endl;
        std::cout << "Indicated Efficiency:             " << std::setw(7) 
                  << indicated_pct << " %" << std::endl;
    }
    
    std::cout << std::string(60, '=') << std::endl;
}

void Engine::output_performance_summary(const std::string& filename) {
    std::ofstream file(filename);
    
    // Write header with parameter names and units
    file << "Parameter,Value,Unit\n";
    
    // Engine operating conditions
    file << "Engine Speed," << engine_speed << ",rpm\n";
    file << "Bore," << bore * 1000.0 << ",mm\n";
    file << "Stroke," << stroke * 1000.0 << ",mm\n";
    file << "Compression Ratio," << compression_ratio << ",-\n";
    file << "Displacement," << (M_PI * bore * bore / 4.0 * stroke) * 1000.0 << ",L\n";
    
    // Performance metrics
    file << "IMEP," << get_indicated_mean_effective_pressure() / 1.0e5 << ",bar\n";
    file << "GMEP," << get_gross_mean_effective_pressure() / 1.0e5 << ",bar\n";
    file << "PMEP," << get_pumping_mean_effective_pressure() / 1.0e5 << ",bar\n";
    file << "BMEP," << get_brake_mean_effective_pressure() / 1.0e5 << ",bar\n";
    file << "FMEP," << friction->get_integrated_fmep(M_PI * bore * bore / 4.0 * stroke) / 1.0e5 << ",bar\n";
    
    // Power outputs
    file << "Indicated Power," << get_indicated_power() / 1000.0 << ",kW\n";
    file << "Brake Power," << get_brake_power() / 1000.0 << ",kW\n";
    
    // Efficiencies
    file << "Cycle Efficiency," << get_cycle_efficiency() * 100.0 << ",%\n";
    file << "Indicated Efficiency," << (indicated_work / total_lhv_energy) * 100.0 << ",%\n";
    
    // Energy ledger (convert to kJ for consistency with printed output)
    file << "Total LHV Energy Available," << total_lhv_energy / 1000.0 << ",kJ\n";
    file << "Indicated Work," << indicated_work / 1000.0 << ",kJ\n";
    file << "Friction Work," << friction_work / 1000.0 << ",kJ\n";
    file << "Pumping Work," << pumping_work / 1000.0 << ",kJ\n";
    file << "Heat Transfer Losses," << heat_transfer_loss / 1000.0 << ",kJ\n";
    file << "Exhaust Energy Losses," << exhaust_energy_loss / 1000.0 << ",kJ\n";
    file << "Net Brake Work," << (indicated_work - friction_work) / 1000.0 << ",kJ\n";
    
    // Energy balance
    double total_accounted = indicated_work + heat_transfer_loss + exhaust_energy_loss;
    double balance_error = total_lhv_energy - total_accounted;
    file << "Total Energy Accounted," << total_accounted / 1000.0 << ",kJ\n";
    file << "Energy Balance Error," << balance_error / 1000.0 << ",kJ\n";
    
    // Fuel properties
    file << "Lower Heating Value," << lower_heating_value / 1e6 << ",MJ/kg\n";
    file << "Fuel Energy Rate," << fuel_energy_rate / 1e3 << ",kW\n";
    
    // Boundary conditions
    file << "Intake Pressure," << intake_pressure / 1e3 << ",kPa\n";
    file << "Intake Temperature," << intake_temperature << ",K\n";
    file << "Exhaust Pressure," << exhaust_pressure / 1e3 << ",kPa\n";
    file << "Wall Temperature," << engine_params.operating_conditions.wall_temperature << ",K\n";
    file << "T3 (Post-Combustion Temperature)," << T3 << ",K\n";
    
    file.close();
    std::cout << "Performance summary written to " << filename << std::endl;
}
