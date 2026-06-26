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
      reference_enthalpy(0.0), T3(0.0) {
    
    // Calculate clearance volume from compression ratio
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    clearance_volume = displacement / (compression_ratio - 1.0);
    
    // Use wall temperature from engine parameters
    double T_wall = engine_params.operating_conditions.wall_temperature;
    
    // Initialize components with parameters from YAML files
    cylinder = std::make_unique<Cylinder>(bore, stroke, clearance_volume, mechanism_file, T_wall);
    
    // Update heat transfer model with parameters from YAML file
    cylinder->update_heat_transfer_parameters(sim_params.heat_transfer.woschni_c1,
                                            sim_params.heat_transfer.woschni_c2,
                                            sim_params.heat_transfer.woschni_c3);
    
    intake_valve = std::make_unique<ValveFlow>(engine_params.valves.intake.diameter, 
                                             engine_params.valves.intake.seat_diameter,
                                             engine_params.valves.intake.discharge_coefficient);
    exhaust_valve = std::make_unique<ValveFlow>(engine_params.valves.exhaust.diameter, 
                                              engine_params.valves.exhaust.seat_diameter,
                                              engine_params.valves.exhaust.discharge_coefficient);
    
    // Pass friction parameters from YAML file
    friction = std::make_unique<FrictionModel>(bore, stroke, con_rod_length,
                                             sim_params.friction.cf_c1,
                                             sim_params.friction.cf_c2,
                                             sim_params.friction.cf_c3,
                                             sim_params.friction.cf_c4);
    
    crankshaft = std::make_unique<Crankshaft>(stroke, con_rod_length, bore);
    
    // Use spark parameters from engine parameters
    spark = std::make_unique<SparkModel>(engine_params.spark.timing_btdc, 
                                       engine_params.spark.duration, 
                                       engine_params.spark.energy);
    
    // Initialize persistent intake gas object
    intake_gas = std::make_unique<CanteraInterface>(mechanism_file);
    
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
    // Convert RPM to rad/s and calculate time step
    double omega = engine_speed * 2.0 * M_PI / 60.0;
    time_step = crank_angle_step * M_PI / 180.0 / omega;
}

void Engine::update_cylinder_volume() {
    double volume = crankshaft->calculate_cylinder_volume(crank_angle, clearance_volume);
    cylinder->set_volume(volume);
}

void Engine::calculate_valve_flows() {
    double intake_lift = intake_valve->get_lift(crank_angle);
    double exhaust_lift = exhaust_valve->get_lift(crank_angle);
    
    // Calculate intake mass flow
    if (intake_lift > 0) {
        double mass_flow = intake_valve->calculate_mass_flow(
            intake_lift, intake_pressure, cylinder->get_pressure(),
            intake_temperature, intake_pressure / (287.0 * intake_temperature));
        
        //std::cout << "intake mass flow: " << mass_flow << std::endl;
        if (mass_flow > 0) {
            // Use persistent intake gas object and update its state
            intake_gas->set_state_TPX(intake_temperature, intake_pressure, intake_composition);
            
            // Pass composition and temperature to the mixing function
            cylinder->add_mass_flow(mass_flow, time_step, intake_composition, intake_temperature);
        }else{
            cylinder->remove_mass_flow(-1.0*mass_flow, time_step);
	}
    }
    
    // Calculate exhaust mass flow - FIXED
    if (exhaust_lift > 0 && cylinder->get_pressure() > exhaust_pressure) {
        double mass_flow = exhaust_valve->calculate_mass_flow(
            exhaust_lift, cylinder->get_pressure(), exhaust_pressure,
            cylinder->get_temperature(), cylinder->get_density());
        
        if (mass_flow > 0) {
            // Calculate exhaust energy loss before removing mass
            double exhaust_enthalpy = cylinder->get_enthalpy(); // [J/kg]
            double enthalpy_difference = exhaust_enthalpy - reference_enthalpy; // [J/kg]
            double energy_loss_rate = mass_flow * enthalpy_difference; // [W]
            
            // Accumulate exhaust energy loss (only positive losses)
            if (enthalpy_difference > 0.0) {
                exhaust_energy_loss += energy_loss_rate * time_step; // [J]
            }
            
            cylinder->remove_mass_flow(mass_flow, time_step);
        }
    }
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
    
    // Calculate instantaneous pumping work
    double instantaneous_pumping_work = 0.0;
    double angle = crank_angle;
    // Intake stroke: 0-180°, Exhaust stroke: 540-720°
    if ((angle >= 0 && angle <= 180) || (angle >= 540 && angle <= 720)) {
        double pumping_pressure = cylinder->get_pressure() - intake_pressure;
        instantaneous_pumping_work = pumping_pressure; // Store pressure difference for now
    }
    results.pumping_works.push_back(instantaneous_pumping_work);
    
    // Store current exhaust energy loss rate
    results.exhaust_energy_losses.push_back(exhaust_energy_loss);
    
    // Store valve flows (simplified)
    results.intake_mass_flows.push_back(0.0);
    results.exhaust_mass_flows.push_back(0.0);
}

void Engine::initialize_cycle() {
    crank_angle = sim_params.simulation.start_angle;  // Use parameter instead of hardcoded -360.0
    results = EngineResults(); // Clear previous results
    
    // Reset energy tracking variables
    exhaust_energy_loss = 0.0;
    heat_transfer_loss = 0.0;
    T3 = 0.0;  // Reset T3 for new cycle
    
    // Initialize cylinder state at TDC exhaust (crank_angle = -360°)
    // This corresponds to the end of the previous exhaust stroke
    cylinder->get_gas()->set_state_TPX(intake_temperature, intake_pressure, intake_composition);
    
    // Update cylinder volume for the starting position
    update_cylinder_volume();
    
    // Update cylinder mass based on the initialized state
    double initial_density = cylinder->get_gas()->density();
    double initial_volume = cylinder->get_volume();
    cylinder->set_mass(initial_density * initial_volume);
    
    // Get species names for CSV header (only need to do this once)
    if (results.species_names.empty()) {
        results.species_names = cylinder->get_gas()->species_names();
    }
    
    // Start friction integration
    friction->start_friction_integration(crank_angle);
    
    // Calculate Lower Heating Value (LHV) for efficiency calculation
    try {
        // Create a temporary gas object for LHV calculation
        auto lhv_gas = std::make_unique<CanteraInterface>(mechanism_file);
        
	//To ensure we are getting LHV and not HHV
	//double lhv_temp = std::max(intake_temperature,600.0);
	//double lhv_pressure = std::max(intake_pressure,1.0e5);
	double lhv_temp = intake_temperature;
	double lhv_pressure = intake_pressure;
        // Set to intake conditions
        lhv_gas->set_state_TPX(lhv_temp, lhv_pressure, intake_composition);
        
        // Store initial enthalpy
        double h_initial = lhv_gas->enthalpy();
        
        // Perform equilibrium combustion at constant temperature and pressure
        lhv_gas->get_solution()->thermo()->equilibrate("TP");
        
        // Calculate enthalpy change (LHV per unit mass of mixture)
        double h_final = lhv_gas->enthalpy();
        lower_heating_value = h_initial - h_final; // J/kg of mixture
        
        // Calculate reference enthalpy at intake conditions for exhaust energy loss calculation
        // Set to intake temperature and pressure
        lhv_gas->set_state_TP(intake_temperature, intake_pressure);
        reference_enthalpy = lhv_gas->enthalpy();
        
    } catch (const std::exception& e) {
        std::cerr << "Warning: LHV calculation failed: " << e.what() << std::endl;
        lower_heating_value = 0.0;
        fuel_energy_rate = 0.0;
    }
}

void Engine::run_cycle() {
    PROFILE_FUNCTION();
    
    initialize_cycle();
    
    // Use simulation parameters instead of hardcoded values
    double target_crank_angle = sim_params.simulation.end_angle;
    double adaptive_angle_step = crank_angle_step;
    double min_angle_step = sim_params.time_stepping.min_angle_step;
    double max_angle_step = sim_params.time_stepping.max_angle_step;
    bool adaptive_stepping = sim_params.time_stepping.adaptive_stepping;
    
    // Progress tracking variables
    double start_angle = sim_params.simulation.start_angle;
    double total_angle_range = target_crank_angle - start_angle;
    int last_progress_percent = -1;
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "Starting engine cycle simulation..." << std::endl;
    std::cout << "Progress: [";
    for (int i = 0; i < 50; ++i) std::cout << " ";
    std::cout << "] 0%" << std::flush;
    
    for (crank_angle = start_angle; crank_angle < target_crank_angle; ) {
        PROFILE_SCOPE("main_loop_iteration");
        // Calculate progress percentage
        double progress = (crank_angle - start_angle) / total_angle_range;
        int progress_percent = static_cast<int>(progress * 100);
        
        //// Update progress bar every 2%
        //if (progress_percent > last_progress_percent && progress_percent % 2 == 0) {
        //    last_progress_percent = progress_percent;
        //    
        //    // Clear the line and redraw progress bar
        //    std::cout << "\rProgress: [";
        //    int filled = static_cast<int>(progress * 50);
        //    for (int i = 0; i < 50; ++i) {
        //        if (i < filled) std::cout << "=";
        //        else if (i == filled) std::cout << ">";
        //        else std::cout << " ";
        //    }
        //    std::cout << "] " << progress_percent << "%" << std::flush;
        //    
        //    // Show current crank angle and phase
        //    std::string phase;
        //    if (crank_angle < 0) phase = "Setup";
        //    else if (crank_angle < 180) phase = "Intake";
        //    else if (crank_angle < 360) phase = "Compression";
        //    else if (crank_angle < 540) phase = "Expansion";
        //    else phase = "Exhaust";
        //    
        //    std::cout << " (" << static_cast<int>(crank_angle) << "° - " << phase << ")" << std::flush;
        //}
        
        // Check if we're in a reactive phase (high temperature or spark active)
        double spark_power = spark->get_spark_power(crank_angle);
        bool reactive_phase = (cylinder->get_temperature() > 1000.0) || 
                             (spark_power > 0.0) ||
                             (spark->is_spark_active(crank_angle));
        
        // Adjust angle step based on conditions (only if adaptive stepping is enabled)
        if (adaptive_stepping) {
            if (reactive_phase) {
                adaptive_angle_step = std::max(min_angle_step, adaptive_angle_step * 0.5);
            } else {
                adaptive_angle_step = std::min(max_angle_step, adaptive_angle_step * 1.1);
            }
        }
        
        // Don't overshoot target
        adaptive_angle_step = std::min(adaptive_angle_step, target_crank_angle - crank_angle);
        
        // Calculate corresponding time step
        double omega = engine_speed * 2.0 * M_PI / 60.0;
        double adaptive_time_step = adaptive_angle_step * M_PI / 180.0 / omega;
        
        // Update time_step for valve flow calculations
        time_step = adaptive_time_step;
        
        // 1. Update cylinder volume due to piston motion
        {
            PROFILE_SCOPE("update_cylinder_volume");
            update_cylinder_volume();
        }
        
        // 2. Calculate and apply mass flows BEFORE thermodynamics
        {
            PROFILE_SCOPE("calculate_valve_flows");
            calculate_valve_flows();
        }
        
        // 3. Handle spark ignition
        if (spark_power > 0.0) {
            PROFILE_SCOPE("spark_ignition");
            // Calculate fuel energy input rate
            // This is the rate at which chemical energy enters the cylinder
            double cycles_per_second = engine_speed / (60.0 * 2.0); // 4-stroke engine
            double mass_flow_rate = cylinder->get_mass() * cycles_per_second;

            fuel_energy_rate = mass_flow_rate * lower_heating_value;
            
            std::cout << "LHV calculated: " << lower_heating_value / 1e6 << " MJ/kg" << std::endl;
            std::cout << "Fuel energy rate: " << fuel_energy_rate / 1000.0 << " kW" << std::endl;
        
            cylinder->handle_spark_ignition(spark_power * adaptive_time_step);
            
            // Capture temperature immediately after spark ignition
            T3 = cylinder->get_temperature();
        }
        
        // 4. Update thermodynamics with correct mass and volume
        {
            PROFILE_SCOPE("update_thermodynamics");
            double mean_piston_speed = crankshaft->calculate_mean_piston_speed(engine_speed);
            cylinder->update_thermodynamics(adaptive_time_step, mean_piston_speed, 0.0);
        }
        
        // 5. Calculate and accumulate friction work
        {
            PROFILE_SCOPE("friction_calculation");
            double piston_velocity = crankshaft->calculate_piston_velocity(crank_angle, engine_speed);
            double friction_power = friction->calculate_instantaneous_friction_power(
                engine_speed, cylinder->get_pressure(), crank_angle, piston_velocity);
            friction->accumulate_friction_work(friction_power, adaptive_time_step);
        }
        
        // 6. Track energy flows for ledger
        {
            PROFILE_SCOPE("energy_tracking");
            
            // Track heat transfer losses
            double heat_loss = cylinder->get_heat_transfer_rate() * adaptive_time_step;
            heat_transfer_loss += std::abs(heat_loss); // Always positive for losses
            
            // Track friction work (already accumulated in friction model)
            // friction_work will be calculated at the end from integrated FMEP
        }
        
        crank_angle += adaptive_angle_step;
        
        // Only store results starting from 0° and at specified intervals
        if (crank_angle >= 0.0 && (results.crank_angles.empty() || 
            crank_angle - results.crank_angles.back() >= sim_params.simulation.output_interval)) {
            store_results();
        }
    }
    
    // Complete the progress bar
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "\rProgress: [";
    for (int i = 0; i < 50; ++i) std::cout << "=";
    std::cout << "] 100% - Complete!" << std::endl;
    std::cout << "Simulation completed in " << duration.count() << " ms" << std::endl;
    
    // Calculate final energy balance for ledger
    calculate_energy_ledger();
    
    // Print profiling report at the end
    #ifdef PROFILING_ENABLED
    Profiler::getInstance().printReport();
    #endif
}

void Engine::run_multiple_cycles(int num_cycles) {
    auto overall_start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_cycles; ++i) {
        std::cout << "\n=== Running cycle " << (i + 1) << " of " << num_cycles << " ===" << std::endl;
        run_cycle();
        
        auto cycle_end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(cycle_end - overall_start);
        
        if (i < num_cycles - 1) {
            double avg_time_per_cycle = elapsed.count() / (double)(i + 1);
            double estimated_remaining = avg_time_per_cycle * (num_cycles - i - 1);
            std::cout << "Estimated time remaining: " << static_cast<int>(estimated_remaining) << " seconds" << std::endl;
        }
    }
    
    auto overall_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(overall_end - overall_start);
    std::cout << "\nAll " << num_cycles << " cycles completed in " << total_duration.count() << " seconds" << std::endl;
}

void Engine::output_results(const std::string& filename) {
    std::ofstream file(filename);
    
    // Write header
    file << "CrankAngle[deg],Pressure[Pa],Temperature[K],Volume[m3],Mass[kg],HeatReleaseRate[W],HeatTransferRate[W],h_conv[W/(m2K)],SurfaceArea[m2],T_gas[K],FrictionPower[W],Gamma[-],PumpingWork[Pa],ExhaustEnergyLoss[J]";
    
    // Add species mass fraction columns
    for (const auto& species : results.species_names) {
        file << ",Y_" << species << "[-]";
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
    // Calculate IMEP from pressure-volume data
    double work = 0.0;
    for (size_t i = 1; i < results.cylinder_pressures.size(); ++i) {
        double dV = results.cylinder_volumes[i] - results.cylinder_volumes[i-1];
        work += results.cylinder_pressures[i] * dV;
    }
    
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    return work / displacement;
}

double Engine::get_brake_mean_effective_pressure() const {
    double imep = get_indicated_mean_effective_pressure();
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    double fmep = friction->get_integrated_fmep(displacement);
    return imep - fmep;
}

double Engine::get_gross_mean_effective_pressure() const {
    // Calculate GMEP from pressure-volume data during compression and expansion strokes only
    // GMEP = work during compression (180-360°) + work during expansion (360-540°)
    double work = 0.0;
    for (size_t i = 1; i < results.cylinder_pressures.size(); ++i) {
        double angle = results.crank_angles[i];
        // Only include compression and expansion strokes (180-540°)
        if (angle >= 180.0 && angle <= 540.0) {
            double dV = results.cylinder_volumes[i] - results.cylinder_volumes[i-1];
            work += results.cylinder_pressures[i] * dV;
        }
    }
    
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    return work / displacement;
}

double Engine::get_pumping_mean_effective_pressure() const {
    // Calculate PMEP from pressure-volume data during intake and exhaust strokes
    // PMEP = work during intake (0-180°) + work during exhaust (540-720°)
    double pumping_work_calculated = 0.0;
    for (size_t i = 1; i < results.cylinder_pressures.size(); ++i) {
        double angle = results.crank_angles[i];
        // Include intake and exhaust strokes (0-180° and 540-720°)
        if ((angle >= 0.0 && angle <= 180.0) || (angle >= 540.0 && angle <= 720.0)) {
            double dV = results.cylinder_volumes[i] - results.cylinder_volumes[i-1];
            pumping_work_calculated += results.cylinder_pressures[i] * dV;
        }
    }
    
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    return pumping_work_calculated / displacement;
}

double Engine::get_pumping_work() const {
    // Calculate absolute pumping work in Joules
    double pmep = get_pumping_mean_effective_pressure();
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    return std::abs(pmep * displacement); // Return absolute value for energy ledger
}

double Engine::get_indicated_power() const {
    double imep = get_indicated_mean_effective_pressure();
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    return imep * displacement * engine_speed / (60.0 * 2.0); // 4-stroke engine
}

double Engine::get_brake_power() const {
    double bmep = get_brake_mean_effective_pressure();
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    return bmep * displacement * engine_speed / (60.0 * 2.0); // 4-stroke engine
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
    if (fuel_energy_rate <= 0.0) return 0.0;
    
    double brake_power = get_brake_power();
    return brake_power / fuel_energy_rate;
}

void Engine::calculate_energy_ledger() {
    // Calculate total LHV energy available per cycle
    double displacement = M_PI * bore * bore / 4.0 * stroke;
    double cycles_per_second = engine_speed / (60.0 * 2.0); // 4-stroke engine
    double cycle_time = 1.0 / cycles_per_second;
    total_lhv_energy = fuel_energy_rate * cycle_time; // Energy per cycle
    
    // Use consistent method call for indicated work calculation
    indicated_work = get_indicated_mean_effective_pressure() * displacement;
    
    // Use the new consistent method for pumping work calculation
    pumping_work = get_pumping_work();
    
    // Get friction work from friction model
    friction_work = friction->get_integrated_fmep(displacement) * displacement;
    
    // Exhaust energy loss is now calculated during valve flow calculations
    // as the integral of enthalpy difference above intake conditions
    
    // Note: exhaust_energy_loss is already accumulated during the cycle
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
    std::cout << "Energy Available (LHV):           " << std::setw(8) << total_kJ << " kJ  (100.0%)" << std::endl;
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
        //double total_accounted = indicated_kJ + friction_kJ + pumping_kJ + heat_transfer_kJ + exhaust_kJ;
	//Indicated work includes friction work
        double total_accounted = indicated_kJ + heat_transfer_kJ + exhaust_kJ;
        double balance_error = total_kJ - total_accounted;
        std::cout << "Total Accounted:                  " << std::setw(8) << total_accounted << " kJ" << std::endl;
        std::cout << "Energy Balance Error:             " << std::setw(8) << balance_error << " kJ" << std::endl;
        
        std::cout << std::string(60, '=') << std::endl;
        std::cout << "Net Brake Work:                   " << std::setw(8) << (indicated_kJ - friction_kJ) << " kJ" << std::endl;
        std::cout << "Brake Efficiency:                 " << std::setw(7) << get_cycle_efficiency() * 100.0 << " %" << std::endl;
        std::cout << "Indicated Efficiency:             " << std::setw(7) << indicated_pct << " %" << std::endl;
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
    file << "Lower Heating Value," << lower_heating_value / 1.0e6 << ",MJ/kg\n";
    file << "Fuel Energy Rate," << fuel_energy_rate / 1000.0 << ",kW\n";
    
    // Boundary conditions
    file << "Intake Pressure," << intake_pressure / 1000.0 << ",kPa\n";
    file << "Intake Temperature," << intake_temperature << ",K\n";
    file << "Exhaust Pressure," << exhaust_pressure / 1000.0 << ",kPa\n";
    file << "Wall Temperature," << engine_params.operating_conditions.wall_temperature << ",K\n";
    file << "T3 (Post-Spark Temperature)," << T3 << ",K\n";
    
    file.close();
    std::cout << "Performance summary written to " << filename << std::endl;
}
