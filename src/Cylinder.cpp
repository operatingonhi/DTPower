#include "Cylinder.h"
#include "Profiler.h"
#include "DebugLogger.h"
#include <cmath>
#include <iostream>
#include <iomanip>

Cylinder::Cylinder(double cylinder_bore, double cylinder_stroke, double clearance_vol,
                   const std::string& mechanism_file, double wall_temp)
    : bore(cylinder_bore), stroke(cylinder_stroke), clearance_volume(clearance_vol),
      wall_temperature(wall_temp), mass(0.0), volume(clearance_vol),
      previous_volume(clearance_vol), combustion_active(false), combustion_initiated(false),
      current_heat_transfer_rate(0.0), instantaneous_heat_release_rate(0.0),
      adiabatic_flame_temp(0.0), target_flame_temp(0.0), flame_radius(0.001),
      mass_unburned(0.0), mass_flame(0.0), mass_burned(0.0),
      mechanism_file(mechanism_file), logging_enabled(false) {
    
    gas = std::make_unique<CanteraInterface>(mechanism_file);
    heat_transfer = std::make_unique<WoschniHeatTransfer>(bore, stroke);
    
    // Initialize persistent mixing gas object
    mixing_gas = std::make_unique<CanteraInterface>(mechanism_file);
    
    // Initialize combustion zone reactors
    unburned_zone = std::make_unique<CanteraInterface>(mechanism_file);
    flame_zone = std::make_unique<CanteraInterface>(mechanism_file);
    burned_zone = std::make_unique<CanteraInterface>(mechanism_file);
    
    // Mass will be set when initialize_cycle() is called
    mass = 0.0;
}

void Cylinder::set_volume(double new_volume) {
    volume = new_volume;
}

void Cylinder::set_wall_temperature(double T_wall) {
    wall_temperature = T_wall;
}

void Cylinder::set_mass(double new_mass) {
    mass = new_mass;
}

void Cylinder::update_heat_transfer_parameters(double c1, double c2, double c3) {
    heat_transfer = std::make_unique<WoschniHeatTransfer>(bore, stroke, c1, c2, c3);
}

void Cylinder::enable_combustion_logging(const std::string& filename) {
    if (combustion_log.is_open()) {
        combustion_log.close();
    }
    
    combustion_log.open(filename);
    if (combustion_log.is_open()) {
        logging_enabled = true;
        // Write CSV header
        combustion_log << "time_step,mass_flow_in,mass_flow_out,rho_flame,flame_radius,"
                       << "flame_shell_thickness,T_flame,target_flame_temp,temp_error\n";
        combustion_log << std::scientific << std::setprecision(6);
    } else {
        std::cerr << "Warning: Could not open combustion log file: " << filename << std::endl;
        logging_enabled = false;
    }
}

void Cylinder::disable_combustion_logging() {
    if (combustion_log.is_open()) {
        combustion_log.close();
    }
    logging_enabled = false;
}

void Cylinder::initiate_spark_ignition(double spark_energy) {
    if (!combustion_initiated && spark_energy > 0.0) {
        // Initialize the three-zone combustion model
        initialize_combustion_zones();
        
        // Calculate adiabatic flame temperature for current conditions
        calculate_adiabatic_flame_temperature();
        
        // Set target flame temperature
        target_flame_temp = flame_temp_fraction * adiabatic_flame_temp;
        
        // Calculate initial flame radius from initial flame kernel mass
        // Assume flame kernel is a sphere with density of flame zone
        double rho_flame = flame_zone->density();
        double flame_volume = mass_flame / rho_flame;
        flame_radius = std::pow(3.0 * flame_volume / (4.0 * M_PI), 1.0/3.0);
        
        // Mark combustion as initiated and active
        combustion_initiated = true;
        combustion_active = true;
        
        std::cout << "Combustion initiated: Tad = " << adiabatic_flame_temp 
                  << " K, Target flame temp = " << target_flame_temp << " K"
                  << ", Initial flame radius = " << flame_radius * 1000.0 << " mm" << std::endl;
    }
}

void Cylinder::initialize_combustion_zones() {
    // Initialize all zones with current cylinder state
    double T = gas->temperature();
    double T_spark = 2000; 
    double T_ign = T_spark;
    double P = gas->pressure();
    std::vector<double> Y = gas->mass_fractions();
    
    // Unburned zone starts with all the mass
    unburned_zone->set_state_TPY(T, P, Y);
    mass_unburned = mass;
    
    // Flame zone starts with a small initial kernel
    flame_zone->set_state_TPY(T_ign, P, Y);
    mass_flame = 1e-9;  // Very small initial mass
    
    // Burned zone starts empty
    burned_zone->set_state_TPY(T, P, Y);
    mass_burned = 0.0;
}

void Cylinder::calculate_adiabatic_flame_temperature() {
    // Create a temporary gas object for adiabatic flame calculation
    auto temp_gas = std::make_unique<CanteraInterface>(mechanism_file);
    
    // Set to current unburned conditions at current cylinder pressure
    double T_u = unburned_zone->temperature();
    double P = gas->pressure();  // Use current cylinder pressure
    std::vector<double> Y_u = unburned_zone->mass_fractions();
    temp_gas->set_state_TPY(T_u, P, Y_u);
    
    // Equilibrate at constant enthalpy and pressure to get adiabatic flame temperature
    try {
        temp_gas->get_solution()->thermo()->equilibrate("HP");
        adiabatic_flame_temp = temp_gas->temperature();
    } catch (const std::exception& e) {
        std::cerr << "Warning: Adiabatic flame temperature calculation failed: " 
                  << e.what() << std::endl;
        adiabatic_flame_temp = T_u * 3.0;  // Rough estimate
    }
}

double Cylinder::calculate_mass_flow_into_flame() {
    // Calculate mass flow rate into flame zone based on temperature control
    // Target: maintain flame zone at target_flame_temp (0.75 * Tad)
    
    if (mass_unburned <= 1e-12) {
        return 0.0;
    }
    
    // Recalculate adiabatic flame temperature at current pressure
    calculate_adiabatic_flame_temperature();
    target_flame_temp = flame_temp_fraction * adiabatic_flame_temp;
    
    double T_flame = flame_zone->temperature();
    double rho_flame = flame_zone->density();
    
    // Control mass flow based on temperature deviation
    // If flame is too cold, reduce mass flow; if too hot, increase mass flow
    double temp_error = T_flame - target_flame_temp;
    
    // Base mass flow rate from flame shell geometry
    // Inner radius of flame shell
    double r_inner = flame_radius;
    // Outer radius of flame shell
    double r_outer = flame_radius + flame_shell_thickness;
    
    // Volume of flame shell
    double V_shell = (4.0/3.0) * M_PI * (r_outer*r_outer*r_outer - r_inner*r_inner*r_inner);
    
    // Estimate characteristic time for flame to traverse shell thickness
    // Use a characteristic velocity based on thermal diffusivity
    double alpha = 1e-4;  // Thermal diffusivity [m²/s], typical for gases
    double characteristic_time = flame_shell_thickness * flame_shell_thickness / alpha;
    
    // Base mass flow rate to maintain shell
    double rho_u = unburned_zone->density();
    double base_mass_flow = rho_u * V_shell / characteristic_time;  // [kg/s]
    
    // Modulate mass flow based on temperature error
    // Use proportional control with gain
    double control_gain = 0.1;  // Tuning parameter
    double temp_factor = 1.0 - control_gain * (temp_error / target_flame_temp);
    
    // Clamp temperature factor to reasonable range
    temp_factor = std::max(0.1, std::min(temp_factor, 2.0));
    
    double mass_flow_in = base_mass_flow * temp_factor;
    
    // Limit mass flow to available unburned mass
    double max_flow = mass_unburned / 1e-4;  // Don't consume all mass in 0.1 ms
    mass_flow_in = std::min(mass_flow_in, max_flow);
    
    // Log combustion variables if logging is enabled
    static int time_step = 0;
    if (logging_enabled && combustion_log.is_open()) {
        // Calculate mass_flow_out here for logging (will be recalculated in update_combustion_zones)
        double mass_flow_out_estimate = calculate_mass_flow_out_of_flame(mass_flow_in);
        
        combustion_log << time_step << ","
                       << mass_flow_in << ","
                       << mass_flow_out_estimate << ","
                       << rho_flame << ","
                       << flame_radius << ","
                       << flame_shell_thickness << ","
                       << T_flame << ","
                       << target_flame_temp << ","
                       << temp_error << "\n";
        time_step++;
    }
    
    return mass_flow_in;
}

double Cylinder::calculate_mass_flow_out_of_flame(double mass_flow_in) {
    // For a spherical shell with uniform thickness, the mass flow out
    // is determined by geometry and the requirement to maintain shell thickness
    
    if (mass_flame <= 1e-12) {
        return mass_flow_in;  // All incoming mass flows through
    }
    
    // Inner radius
    double r_inner = flame_radius;
    // Outer radius
    double r_outer = flame_radius + flame_shell_thickness;
    
    // Volume of shell
    double V_shell = (4.0/3.0) * M_PI * (r_outer*r_outer*r_outer - r_inner*r_inner*r_inner);
    
    // Density in flame zone
    double rho_flame = flame_zone->density();
    
    // Target mass in shell to maintain uniform thickness
    double target_mass_shell = rho_flame * V_shell;
    
    // Mass flow out needed to maintain target mass
    // If current mass > target, increase outflow; if < target, decrease outflow
    double mass_error = mass_flame - target_mass_shell;
    
    // Base outflow equals inflow (steady state)
    double mass_flow_out = mass_flow_in;
    
    // Add correction term to drive toward target mass
    double correction_gain = 10.0;  // Tuning parameter
    mass_flow_out += correction_gain * mass_error;
    
    // Ensure non-negative and don't remove more mass than available
    mass_flow_out = std::max(0.0, std::min(mass_flow_out, mass_flame / 1e-4));
    
    return mass_flow_out;
}

void Cylinder::update_combustion_zones(double dt) {
    std::cout << "[DEBUG] update_combustion_zones: mass_unburned=" << mass_unburned 
              << ", mass_flame=" << mass_flame 
              << ", mass_burned=" << mass_burned << std::endl;
    
    if (!combustion_active || mass_unburned <= 1e-12) {
        std::cout << "[DEBUG] Combustion not active or no unburned mass remaining" << std::endl;
        combustion_active = false;
        instantaneous_heat_release_rate = 0.0;
        return;
    }
    
    // Calculate change in specific internal energy due to compression/expansion
    // This represents the work done on/by the gas during volume change
    double u_gas_before = gas->internal_energy();  // [J/kg]
    double P_gas_before = gas->pressure();
    
    // Get current pressure after volume change (from main gas object)
    double P = gas->pressure();
    
    // Calculate change in specific internal energy
    // For an isentropic process: du = -P*dv where v is specific volume
    double specific_volume_before = previous_volume / mass;  // [m³/kg]
    double specific_volume_now = volume / mass;  // [m³/kg]
    double delta_u = u_gas_before - gas->internal_energy();  // [J/kg]
    
    // Update all zones with new pressure AND internal energy change
    double T_u = unburned_zone->temperature();
    std::vector<double> Y_u = unburned_zone->mass_fractions();
    double u_u = unburned_zone->internal_energy();
    double v_u = 1.0 / unburned_zone->density();
    unburned_zone->set_state_UV(u_u + delta_u, specific_volume_now);
    
    double T_f = flame_zone->temperature();
    std::vector<double> Y_f = flame_zone->mass_fractions();
    double u_f = flame_zone->internal_energy();
    double v_f = 1.0 / flame_zone->density();
    flame_zone->set_state_UV(u_f + delta_u, specific_volume_now);
    
    double T_b = burned_zone->temperature();
    std::vector<double> Y_b = burned_zone->mass_fractions();
    double u_b = burned_zone->internal_energy();
    double v_b = 1.0 / burned_zone->density();
    burned_zone->set_state_UV(u_b + delta_u, specific_volume_now);
    
    // Calculate mass flow rates
    std::cout << "[DEBUG] Calculating mass flow rates..." << std::endl;
    double mass_flow_in = calculate_mass_flow_into_flame();  // [kg/s]
    std::cout << "[DEBUG] mass_flow_in=" << mass_flow_in << " kg/s" << std::endl;
    
    double mass_flow_out = calculate_mass_flow_out_of_flame(mass_flow_in);  // [kg/s]
    std::cout << "[DEBUG] mass_flow_out=" << mass_flow_out << " kg/s" << std::endl;
    
    // Calculate mass transfers
    double dm_in = mass_flow_in * dt;  // [kg]
    double dm_out = mass_flow_out * dt;  // [kg]
    std::cout << "[DEBUG] dm_in=" << dm_in << " kg, dm_out=" << dm_out << " kg" << std::endl;
    
    // Limit transfers to available mass
    dm_in = std::min(dm_in, mass_unburned);
    dm_out = std::min(dm_out, mass_flame);
    
    if (dm_in <= 0.0 && dm_out <= 0.0) {
    	std::cout << "[DEBUG] combustion stoped by negatvie flow rate" << std::endl;
        combustion_active = false;
        instantaneous_heat_release_rate = 0.0;
        return;
    }
    
    // Store initial enthalpies for heat release calculation
    double h_u_initial = unburned_zone->enthalpy();
    double h_f_initial = flame_zone->enthalpy();
    double h_b_initial = burned_zone->enthalpy();
    double total_enthalpy_initial = h_u_initial * mass_unburned + 
                                    h_f_initial * mass_flame + 
                                    h_b_initial * mass_burned;
    
    // Transfer mass from unburned to flame zone
    if (dm_in > 1e-15) {
        mass_unburned -= dm_in;
        
        // Mix incoming unburned gas into flame zone
        double mass_flame_old = mass_flame;
        mass_flame += dm_in;
        
        if (mass_flame > 1e-12) {
            // Mass-weighted mixing of compositions
            std::vector<double> Y_flame_new(Y_f.size());
            for (size_t i = 0; i < Y_f.size(); ++i) {
                Y_flame_new[i] = (Y_f[i] * mass_flame_old + Y_u[i] * dm_in) / mass_flame;
            }
            
            // Mass-weighted mixing of enthalpies
            double h_flame_mixed = (h_f_initial * mass_flame_old + h_u_initial * dm_in) / mass_flame;
            
            // Set flame zone state with mixed composition at current pressure
            flame_zone->set_state_TPY(T_f, P, Y_flame_new);
            
            // Iterate to match mixed enthalpy
            for (int iter = 0; iter < 5; ++iter) {
                double h_current = flame_zone->enthalpy();
                double T_current = flame_zone->temperature();
                double cp = flame_zone->cp();
                double T_new = T_current + (h_flame_mixed - h_current) / cp;
                T_new = std::max(300.0, std::min(T_new, 5000.0));
                flame_zone->set_state_TPY(T_new, P, Y_flame_new);
            }
        }
    }
    
    // Advance chemistry in flame zone (reactions progress here)
    if (mass_flame > 1e-12) {
        std::cout << "[DEBUG] Advancing flame zone chemistry: T=" << flame_zone->temperature() 
                  << " K, P=" << flame_zone->pressure() / 1e5 << " bar" << std::endl;
        flame_zone->advance_chemistry(dt);
        std::cout << "[DEBUG] Flame zone chemistry complete: T=" << flame_zone->temperature() 
                  << " K" << std::endl;
    }
    
    // Transfer mass from flame to burned zone
    if (dm_out > 1e-15) {
        // Get current flame zone state after chemistry
        double h_f_after_chem = flame_zone->enthalpy();
        std::vector<double> Y_f_after_chem = flame_zone->mass_fractions();
        
        mass_flame -= dm_out;
        
        // Mix outgoing flame gas into burned zone
        double mass_burned_old = mass_burned;
        mass_burned += dm_out;
        
        if (mass_burned > 1e-12) {
            // Mass-weighted mixing of compositions
            std::vector<double> Y_burned_new(Y_b.size());
            for (size_t i = 0; i < Y_b.size(); ++i) {
                Y_burned_new[i] = (Y_b[i] * mass_burned_old + Y_f_after_chem[i] * dm_out) / mass_burned;
            }
            
            // Mass-weighted mixing of enthalpies
            double h_burned_mixed = (h_b_initial * mass_burned_old + h_f_after_chem * dm_out) / mass_burned;
            
            // Set burned zone state with mixed composition at current pressure
            burned_zone->set_state_TPY(T_b, P, Y_burned_new);
            
            // Iterate to match mixed enthalpy
            for (int iter = 0; iter < 5; ++iter) {
                double h_current = burned_zone->enthalpy();
                double T_current = burned_zone->temperature();
                double cp = burned_zone->cp();
                double T_new = T_current + (h_burned_mixed - h_current) / cp;
                T_new = std::max(300.0, std::min(T_new, 5000.0));
                burned_zone->set_state_TPY(T_new, P, Y_burned_new);
            }
        }
    }
    
    // Calculate heat release rate from enthalpy change
    double h_u_final = unburned_zone->enthalpy();
    double h_f_final = flame_zone->enthalpy();
    double h_b_final = burned_zone->enthalpy();
    double total_enthalpy_final = h_u_final * mass_unburned + 
                                  h_f_final * mass_flame + 
                                  h_b_final * mass_burned;
    
    double heat_release = total_enthalpy_final - total_enthalpy_initial;  // [J]
    instantaneous_heat_release_rate = heat_release / dt;  // [W]
    
    // Update flame radius based on burned volume
    double total_volume = volume;
    double burned_volume_fraction = mass_burned / mass;
    double burned_volume = total_volume * burned_volume_fraction;
    flame_radius = std::pow(3.0 * burned_volume / (4.0 * M_PI), 1.0/3.0);
    
    // Check if combustion is complete
    if (mass_unburned < 1e-12 && mass_flame < 1e-12) {
        combustion_active = false;
        std::cout << "[DEBUG] Combustion complete: all mass burned" << std::endl;
    }
    
    std::cout << "[DEBUG] update_combustion_zones complete" << std::endl;
}

void Cylinder::mix_zones_to_cylinder() {
    // Combine all three zones back into the main cylinder gas object
    // This maintains the overall cylinder state
    
    if (!combustion_initiated) {
        return;  // No combustion zones to mix
    }
    
    double total_mass = mass_unburned + mass_flame + mass_burned;
    
    if (total_mass < 1e-12) {
        return;
    }
    
    // Calculate mass-weighted average composition
    std::vector<double> Y_avg(gas->n_species(), 0.0);
    std::vector<double> Y_u = unburned_zone->mass_fractions();
    std::vector<double> Y_f = flame_zone->mass_fractions();
    std::vector<double> Y_b = burned_zone->mass_fractions();
    
    for (size_t i = 0; i < Y_avg.size(); ++i) {
        Y_avg[i] = (Y_u[i] * mass_unburned + Y_f[i] * mass_flame + Y_b[i] * mass_burned) / total_mass;
    }
    
    // Calculate mass-weighted average internal energy
    double u_avg = (unburned_zone->internal_energy() * mass_unburned + 
                    flame_zone->internal_energy() * mass_flame + 
                    burned_zone->internal_energy() * mass_burned) / total_mass;
    
    // Calculate specific volume
    double v_specific = volume / total_mass;
    
    // Set cylinder state using internal energy and specific volume
    // This allows pressure to naturally emerge from the thermodynamic state
    gas->set_state_UV(u_avg, v_specific);
    
    // Update composition to match mass-weighted average
    // We need to do this after setting UV to ensure correct species distribution
    double T_current = gas->temperature();
    double P_current = gas->pressure();
    gas->set_state_TPY(T_current, P_current, Y_avg);
    
    // Re-apply UV constraint to maintain energy and volume with new composition
    gas->set_state_UV(u_avg, v_specific);
}

void Cylinder::update_thermodynamics(double dt, double mean_piston_speed, double motored_pressure, bool valves_flowing) {
    PROFILE_FUNCTION();
    
    // Store initial state (from PREVIOUS time step)
    double initial_volume = previous_volume;  // Use previous volume, not current
    double initial_T = gas->temperature();
    double initial_P = gas->pressure();
    double initial_S = gas->entropy();  // Store entropy for isentropic process
    std::vector<double> initial_Y = gas->mass_fractions();
    
    DEBUG_CHECK_TEMPERATURE(initial_T, "Initial state");
    DEBUG_CHECK_PRESSURE(initial_P, "Initial state");
    DEBUG_LOG("Initial: T=" + std::to_string(initial_T) + " K, P=" + std::to_string(initial_P / 1e5) + " bar, V_ratio=" + std::to_string(volume / initial_volume));
    
    // Advance chemistry first
    advance_chemistry(dt);
    
    // Get state after chemistry
    double T_after_chemistry = gas->temperature();
    double P_after_chemistry = gas->pressure();
    double S_after_chemistry = gas->entropy();
    
    DEBUG_CHECK_TEMPERATURE(T_after_chemistry, "After chemistry");
    DEBUG_CHECK_PRESSURE(P_after_chemistry, "After chemistry");
    
    // Apply isentropic compression/expansion due to volume change
    double volume_ratio = volume / initial_volume;
    
    if (std::abs(volume_ratio - 1.0) > 1e-6) {
        PROFILE_SCOPE("isentropic_process");
        
        // Apply isentropic compression/expansion from the post-chemistry state
        // The entropy after chemistry already includes chemical reaction effects
        double new_density = mass / volume;
        double specific_volume = 1.0 / new_density;
        
        // Use Cantera's setState_SV (entropy, specific volume) to maintain isentropic process
        gas->set_state_SV(S_after_chemistry, specific_volume);
    }
    
    // Get state after isentropic process
    double T_after_isentropic = gas->temperature();
    double P_after_isentropic = gas->pressure();
    
    DEBUG_CHECK_TEMPERATURE(T_after_isentropic, "After isentropic");
    DEBUG_CHECK_PRESSURE(P_after_isentropic, "After isentropic");
    
    // Calculate heat transfer losses
    {
        PROFILE_SCOPE("heat_transfer");
        
        double h_conv = heat_transfer->calculate_heat_transfer_coefficient(
            P_after_isentropic, T_after_isentropic, mean_piston_speed, motored_pressure);
        double surface_area = heat_transfer->calculate_surface_area(volume, bore, stroke, clearance_volume);
        double heat_flux = heat_transfer->calculate_heat_flux(
            T_after_isentropic, wall_temperature, surface_area, h_conv);
        
        // Apply heat transfer effects at constant volume
        double net_heat_rate = -heat_flux; // J/s (negative = heat loss)
        
        // Store the heat transfer rate for output
        current_heat_transfer_rate = net_heat_rate;
        
        // Apply heat transfer at constant volume by updating internal energy
        double current_internal_energy = gas->internal_energy() * mass;  // Total internal energy [J]
        double new_internal_energy = current_internal_energy + net_heat_rate * dt;  // [J]
        double specific_internal_energy = new_internal_energy / mass;  // [J/kg]
        double specific_volume = volume / mass;  // [m³/kg]
        
        // Use setState_UV (internal energy, specific volume) for constant volume heat transfer
        gas->set_state_UV(specific_internal_energy, specific_volume);
            
    }
    
    // Final checks
    double final_T = gas->temperature();
    double final_P = gas->pressure();
    
    DEBUG_CHECK_TEMPERATURE(final_T, "Final state");
    DEBUG_CHECK_PRESSURE(final_P, "Final state");
    
    if (final_T > 5000.0) {
        DEBUG_LOG("EXTREME TEMPERATURE DETECTED!");
        DEBUG_LOG("Volume ratio: " + std::to_string(volume_ratio));
    }
    
    // Update previous_volume for next time step
    previous_volume = volume;
}

void Cylinder::add_mass_flow(double mass_flow_rate, double dt, const std::string& incoming_composition, double incoming_temperature) {
    double dm = mass_flow_rate * dt;
    if (dm <= 0.0) return;
    
    try {
        // Use persistent mixing gas object instead of creating new one
        double current_pressure = gas->pressure();
        mixing_gas->set_state_TPX(incoming_temperature, current_pressure, incoming_composition);
        
        // Now mix the two gas streams at constant pressure
        double current_mass = mass;
        double new_total_mass = current_mass + dm;
        
        // Get current cylinder gas properties
        double current_enthalpy = gas->enthalpy();
        auto Y_current = gas->mass_fractions();
        auto Y_incoming = mixing_gas->mass_fractions();
        
        // Mass-weighted mixing for energy conservation
        double incoming_enthalpy = mixing_gas->enthalpy();
        double mixed_enthalpy = (current_mass * current_enthalpy + dm * incoming_enthalpy) / new_total_mass;
        
        // Mass-weighted mixing of species compositions
        std::vector<double> Y_mixed(Y_current.size());
        for (size_t i = 0; i < Y_current.size(); ++i) {
            Y_mixed[i] = (current_mass * Y_current[i] + dm * Y_incoming[i]) / new_total_mass;
        }
        
        // Normalize mass fractions to ensure they sum to 1.0
        double sum_Y = 0.0;
        for (double& yi : Y_mixed) sum_Y += yi;
        if (sum_Y > 0.0) {
            for (double& yi : Y_mixed) yi /= sum_Y;
        }
        
        // Update mass
        mass = new_total_mass;
        
        // Find temperature that gives mixed enthalpy with mixed composition at constant pressure
        double T_guess = gas->temperature();
        double P_mixing = current_pressure;
        
        // Iterate to find temperature that gives the target mixed enthalpy
        const int max_iterations = 20;
        for (int iter = 0; iter < max_iterations; iter++) {
            gas->set_state_TPY(T_guess, P_mixing, Y_mixed);
            double h_current = gas->enthalpy();
            double error = h_current - mixed_enthalpy;
            
            const double convergence_tolerance = 1e-6;
            if (std::abs(error) < convergence_tolerance * std::abs(mixed_enthalpy)) {
                break;  // Converged
            }
            
            // Calculate derivative dh/dT at constant pressure
            double dT_test = 1.0;
            gas->set_state_TPY(T_guess + dT_test, P_mixing, Y_mixed);
            double h_plus = gas->enthalpy();
            double dh_dT = (h_plus - h_current) / dT_test;
            
            // Newton-Raphson update
            const double min_derivative = 1e-10;
            const double min_temperature = 200.0;   // [K]
            const double max_temperature = 4000.0;  // [K]
            
            if (std::abs(dh_dT) > min_derivative) {
                T_guess -= error / dh_dT;
                T_guess = std::max(min_temperature, std::min(T_guess, max_temperature));
            } else {
                // Fallback: proportional adjustment
                T_guess *= mixed_enthalpy / h_current;
                T_guess = std::max(min_temperature, std::min(T_guess, max_temperature));
                break;
            }
        }
        
        // Final pressure adjustment for new density in the cylinder volume
        double new_density = mass / volume;
        double final_temperature = gas->temperature();
        double R_specific = 8314.0 / gas->mean_molecular_weight();
        double final_pressure = new_density * R_specific * final_temperature;
        
        // Set final mixed state
        gas->set_state_TPY(final_temperature, final_pressure, Y_mixed);
        
        // Reset combustion state when new mass is added
        combustion_initiated = false;
        combustion_active = false;
        
    } catch (const std::exception& e) {
        std::cerr << "Warning: Advanced mixing failed: " << e.what() << std::endl;
        
        // Fallback: assume incoming gas has same composition as current gas
        double current_mass = mass;
        double new_total_mass = current_mass + dm;
        double current_enthalpy = gas->enthalpy();
        double mixed_enthalpy = (current_mass * current_enthalpy + dm * current_enthalpy) / new_total_mass;
        
        mass = new_total_mass;
        
        // Estimate temperature from mixed enthalpy
        const double min_temperature = 200.0;   // [K]
        const double max_temperature = 4000.0;  // [K]
        double mixed_T = mixed_enthalpy / gas->cp();
        mixed_T = std::max(min_temperature, std::min(mixed_T, max_temperature));
        
        // Calculate new pressure from ideal gas law
        double new_density = mass / volume;
        double R_specific = 8314.0 / gas->mean_molecular_weight();
        double new_P = new_density * R_specific * mixed_T;
        
        gas->set_state_TPY(mixed_T, new_P, gas->mass_fractions());
    }
}

void Cylinder::remove_mass_flow(double mass_flow_rate, double dt) {
    double dm = mass_flow_rate * dt;
    if (dm <= 0.0) return;
    
    // Store current state
    double current_mass = mass;
    double current_enthalpy = gas->enthalpy();  // Specific enthalpy [J/kg]
    auto current_Y = gas->mass_fractions();
    
    // Calculate new mass (prevent excessive mass removal)
    const double min_mass_fraction = 0.1;
    double new_mass = std::max(current_mass * min_mass_fraction, current_mass - dm);
    
    // Calculate enthalpy balance
    double total_enthalpy_before = current_mass * current_enthalpy;  // [J]
    double enthalpy_removed = dm * current_enthalpy;  // [J]
    double total_enthalpy_after = total_enthalpy_before - enthalpy_removed;  // [J]
    
    // Calculate new specific enthalpy of remaining gas
    double new_specific_enthalpy = total_enthalpy_after / new_mass;  // [J/kg]
    
    // Update mass
    mass = new_mass;
    
    try {
        // Set state using enthalpy and pressure
        double current_pressure = gas->pressure();
        gas->set_state_HP(new_specific_enthalpy, current_pressure);
        
    } catch (const std::exception& e) {
        std::cerr << "Warning: Mass removal failed: " << e.what() << std::endl;
        
        // Fallback: maintain temperature, recalculate pressure for new density
        double current_T = gas->temperature();
        double new_density = mass / volume;
        double R_specific = 8314.0 / gas->mean_molecular_weight();
        double new_P = new_density * R_specific * current_T;
        
        gas->set_state_TPY(current_T, new_P, current_Y);
    }
}

void Cylinder::advance_chemistry(double dt) {
    std::cout << "[DEBUG] Cylinder::advance_chemistry called: dt=" << dt 
              << ", combustion_active=" << combustion_active 
              << ", combustion_initiated=" << combustion_initiated << std::endl;
    
    // Update combustion zones if combustion is active
    if (combustion_active) {
        std::cout << "[DEBUG] Updating combustion zones..." << std::endl;
        update_combustion_zones(dt);
        std::cout << "[DEBUG] Mixing zones to cylinder..." << std::endl;
        mix_zones_to_cylinder();
        std::cout << "[DEBUG] Combustion zones updated successfully" << std::endl;
    } else if (!combustion_initiated) {
        // Only advance chemistry if combustion hasn't started yet
        std::cout << "[DEBUG] Advancing chemistry (no combustion)..." << std::endl;
        gas->advance_chemistry(dt);
        std::cout << "[DEBUG] Chemistry advanced successfully" << std::endl;
    } else {
        std::cout << "[DEBUG] Skipping chemistry (combustion initiated but not active)" << std::endl;
    }
}


double Cylinder::get_temperature() const {
    return gas->temperature();
}

double Cylinder::get_pressure() const {
    return gas->pressure();
}

double Cylinder::get_density() const {
    return gas->density();
}

double Cylinder::get_mass() const {
    return mass;
}

double Cylinder::get_volume() const {
    return volume;
}

double Cylinder::get_internal_energy() const {
    return gas->internal_energy();
}

double Cylinder::get_enthalpy() const {
    return gas->enthalpy();
}

double Cylinder::get_heat_transfer_rate() const {
    return current_heat_transfer_rate;
}

double Cylinder::get_heat_release_rate() {
    // Return and reset instantaneous heat release rate
    double heat_release = instantaneous_heat_release_rate;
    instantaneous_heat_release_rate = 0.0;
    return heat_release;
}

double Cylinder::get_gamma() const {
    return gas->gamma();
}
