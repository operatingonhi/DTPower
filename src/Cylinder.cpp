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
    heat_transfer = std::make_unique<WoschniHeatTransfer>(bore, stroke, 3.26, 2.28, 3.24e-3);
    
    // Initialize persistent mixing gas object
    mixing_gas = std::make_unique<CanteraInterface>(mechanism_file);
    
    // Initialize persistent adiabatic gas object for flame temperature calculations
    adiabatic_gas = std::make_unique<CanteraInterface>(mechanism_file);
    
    // Initialize combustion zone reactors
    unburned_zone = std::make_unique<CanteraInterface>(mechanism_file);
    flame_zone = std::make_unique<CanteraInterface>(mechanism_file);
    burned_zone = std::make_unique<CanteraInterface>(mechanism_file);
    
    // Mass will be set when initialize_cycle() is called
    mass = 0.0;
}

void Cylinder::set_volume(double new_volume) {
    
    // Apply isentropic compression/expansion due to volume change
    double volume_ratio = new_volume / volume;
    
    if (std::abs(volume_ratio - 1.0) > 0.0) {
        PROFILE_SCOPE("isentropic_process");
        
        // Apply isentropic compression/expansion from the post-chemistry state
        // The entropy after chemistry already includes chemical reaction effects
        double new_density = mass / new_volume;
        double specific_volume = 1.0 / new_density;
        
        // Use Cantera's setState_SV (entropy, specific volume) to maintain isentropic process
        gas->set_state_SV(gas->entropy(), specific_volume);
    }
    previous_volume = volume;
    volume = new_volume;
    
}

void Cylinder::set_wall_temperature(double T_wall) {
    wall_temperature = T_wall;
}

void Cylinder::set_mass(double new_mass) {
    mass = new_mass;
}

void Cylinder::update_heat_transfer_parameters(double c1, double c2, double c3) {
    // Store current reference values if they were set
    bool had_reference_values = false;
    double old_Tr = 300.0, old_Pr = 100.0, old_Vr = 0.0;
    
    if (heat_transfer) {
        // Note: We can't directly access private members, so we'll lose reference values
        // This is acceptable since update_heat_transfer_parameters is typically called
        // before the simulation starts
    }
    
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
                       << "flame_shell_thickness,T_flame,target_flame_temp,temp_error,"
                       << "mass_unburned,mass_burned\n";
        combustion_log << std::scientific << std::setprecision(6);
        combustion_log.flush();  // Ensure header is written immediately
        std::cout << "Combustion logging enabled: " << filename << std::endl;
    } else {
        std::cerr << "Warning: Could not open combustion log file: " << filename << std::endl;
        logging_enabled = false;
    }
}

void Cylinder::disable_combustion_logging() {
    if (combustion_log.is_open()) {
        combustion_log.flush();  // Ensure all data is written before closing
        combustion_log.close();
        std::cout << "Combustion logging disabled and file closed" << std::endl;
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
    }
}

void Cylinder::initialize_combustion_zones() {
    // Initialize all zones with current cylinder state
    double T = gas->temperature();
    double T_spark = 2000; 
    double T_ign = T_spark;
    double sparkm_ratio = 0.001;
    double P = gas->pressure();
    std::vector<double> Y = gas->mass_fractions();
    
    // Unburned zone starts with all the mass
    unburned_zone->set_state_TPY(T, P, Y);
    double rho_u = unburned_zone->density();
    mass_unburned = mass*(1-sparkm_ratio);
    
    // Flame zone starts with a small initial kernel
    flame_zone->set_state_TPY(T_ign, P, Y);
    mass_flame = mass*sparkm_ratio;  // Very small initial mass
    
    // Burned zone starts empty
    burned_zone->set_state_TPY(T, P, Y);
    mass_burned = 0.0;
}

void Cylinder::calculate_adiabatic_flame_temperature() {
    // Reuse persistent gas object for adiabatic flame calculation
    // Set to current unburned conditions at current cylinder pressure
    double T_u = unburned_zone->temperature();
    double P = gas->pressure();  // Use current cylinder pressure
    std::vector<double> Y_u = unburned_zone->mass_fractions();
    adiabatic_gas->set_state_TPY(T_u, P, Y_u);
    
    // Equilibrate at constant enthalpy and pressure to get adiabatic flame temperature
    adiabatic_gas->get_solution()->thermo()->equilibrate("HP");
    adiabatic_flame_temp = adiabatic_gas->temperature();
}

double Cylinder::calculate_mass_flow_into_flame(double dt) {
    // Calculate mass flow rate into flame zone based on thermal mass balance
    // Target: maintain flame zone at target_flame_temp (0.75 * Tad)
    
    if (mass_unburned <= 1e-12) {
        return 0.0;
    }
    
    // Recalculate adiabatic flame temperature at current pressure
    calculate_adiabatic_flame_temperature();
    target_flame_temp = flame_temp_fraction * adiabatic_flame_temp;
    
    double T_flame = flame_zone->temperature();
    double T_unburned = unburned_zone->temperature();
    
    // Temperature error: positive means flame is too hot, negative means too cold
    double temp_error = T_flame - target_flame_temp;
    
    // Get thermal properties
    double cp_flame = flame_zone->cp();  // Specific heat of flame zone [J/(kg·K)]
    double cp_unburned = unburned_zone->cp();  // Specific heat of unburned zone [J/(kg·K)]
    
    // Thermal mass of flame zone [J/K]
    double thermal_mass_flame = mass_flame * cp_flame;
    
    // Energy needed to correct temperature error [J]
    // Negative energy_needed means we need to cool the flame (add cold gas)
    // Positive energy_needed means we need to heat the flame (add less gas or hotter gas)
    double energy_needed = -thermal_mass_flame * temp_error;
    
    // Energy per unit mass that incoming unburned gas can provide [J/kg]
    // This is the enthalpy difference between unburned and flame temperatures
    double energy_per_kg = cp_unburned * (T_unburned - T_flame);
    
    // Characteristic time scale for mass flow [s]
    // Based on thermal diffusivity and flame shell thickness
    double alpha = 1e-4;  // Thermal diffusivity [m²/s]
    double characteristic_time = flame_shell_thickness * flame_shell_thickness / alpha;
    
    // Calculate mass flow rate needed to provide the required energy over characteristic time
    double mass_flow_in = 0.0;
    
    if (std::abs(energy_per_kg) > 1e-6) {
        // Mass flow rate [kg/s] = Energy rate [J/s] / Energy per kg [J/kg]
        double energy_rate = energy_needed / dt;  // [J/s]
        mass_flow_in = energy_rate / energy_per_kg;  // [kg/s]
    } else {
        // Fallback: if temperature difference is negligible, use geometric estimate
        double r_inner = flame_radius;
        double r_outer = flame_radius + flame_shell_thickness;
        double V_shell = (4.0/3.0) * M_PI * (r_outer*r_outer*r_outer - r_inner*r_inner*r_inner);
        double rho_u = unburned_zone->density();
        mass_flow_in = rho_u * V_shell / characteristic_time;
    }
    
    // Ensure non-negative mass flow
    mass_flow_in = std::max(0.0, mass_flow_in);
    
    // Limit mass flow to available unburned mass or half of flame zone
    double max_flow = mass_unburned / 1e-4;  // Don't consume all mass in 0.1 ms
    mass_flow_in = std::min(mass_flow_in, max_flow);
    mass_flow_in = std::min(mass_flow_in, mass_flame/(2.0*dt));
    
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
    if (!combustion_active || mass_unburned <= 1e-12) {
        combustion_active = false;
        instantaneous_heat_release_rate = 0.0;
        return;
    }
    
    // Get current pressure after volume change (from main gas object)
    double P = gas->pressure();
    
    // Cache densities before state changes
    double rho_u = unburned_zone->density();
    double rho_f = flame_zone->density();
    double rho_b = burned_zone->density();
    double old_volume = mass_unburned/rho_u + mass_burned/rho_b + mass_flame/rho_f;
    if ((old_volume-previous_volume)/previous_volume>0.01){
        std::cout << "Warning: unburned_density " << rho_u << " density " << gas->density() << " calculation " << mass/volume << std::endl;
    } else{
	old_volume=previous_volume;
    }
    
    // Get mass fractions once before state updates
    std::vector<double> Y_u = unburned_zone->mass_fractions();
    std::vector<double> Y_f = flame_zone->mass_fractions();
    std::vector<double> Y_b = burned_zone->mass_fractions();
    
    // Update all zones with new pressure AND specific volume change
    unburned_zone->set_state_DP((rho_u/volume*old_volume),P);
    flame_zone->set_state_DP((rho_f/volume*old_volume),P);
    burned_zone->set_state_DP((rho_b/volume*old_volume),P);
    
    // Update densities after DP state change
    rho_u = unburned_zone->density();
    rho_f = flame_zone->density();
    rho_b = burned_zone->density();

    // Cache internal energies
    double u_u = unburned_zone->internal_energy();
    double u_f = flame_zone->internal_energy();
    double u_b = burned_zone->internal_energy();
    
    double total_mass = mass_unburned + mass_burned + mass_flame;
    double cumulative_internal_energy = (mass_unburned*u_u + mass_burned*u_b + mass_flame*u_f)/total_mass;
    double u_diff = gas->internal_energy() - cumulative_internal_energy;
    
    // Apply corrections in one pass
    unburned_zone->set_state_UV((u_u + u_diff), 1.0/rho_u);
    flame_zone->set_state_UV((u_f + u_diff), 1.0/rho_f);
    burned_zone->set_state_UV((u_b + u_diff), 1.0/rho_b);


    // Calculate mass flow rates
    double mass_flow_in = calculate_mass_flow_into_flame(dt);  // [kg/s]
    
    double mass_flow_out = calculate_mass_flow_out_of_flame(mass_flow_in);  // [kg/s]

    // Log combustion variables if logging is enabled
    static int time_step = 0;
    if (logging_enabled && combustion_log.is_open()) {
        double T_flame = flame_zone->temperature();
        double temp_error = T_flame - target_flame_temp;
        
        combustion_log << time_step << ","
                       << mass_flow_in << ","
                       << mass_flow_out << ","
                       << rho_f << ","
                       << flame_radius << ","
                       << flame_shell_thickness << ","
                       << T_flame << ","
                       << target_flame_temp << ","
                       << temp_error << ","
                       << mass_unburned << ","
                       << mass_burned << "\n";
        combustion_log.flush();  // Flush after each write to ensure data is saved
        time_step++;
        
    }
    
    
    // Calculate mass transfers
    double dm_in = mass_flow_in * dt;  // [kg]
    double dm_out = mass_flow_out * dt;  // [kg]
    
    // Limit transfers to available mass
    dm_in = std::min(dm_in, mass_unburned);
    dm_out = std::min(dm_out, mass_flame);
    
    if (dm_in <= 0.0 && dm_out <= 0.0) {
        combustion_active = false;
        instantaneous_heat_release_rate = 0.0;
        return;
    }
    
    // Store initial enthalpies for heat release calculation (cache these values)
    double h_u_initial = unburned_zone->enthalpy();
    double h_f_initial = flame_zone->enthalpy();
    double h_b_initial = burned_zone->enthalpy();
    double total_enthalpy_initial = h_u_initial * mass_unburned + 
                                    h_f_initial * mass_flame + 
                                    h_b_initial * mass_burned;
    double u_u_initial = unburned_zone->internal_energy();
    double u_f_initial = flame_zone->internal_energy();
    double u_b_initial = burned_zone->internal_energy();
    
    // Transfer mass from flame to burned zone
    if (dm_out > 1e-15) {
        
        mass_flame -= dm_out;
        
        // Mix outgoing flame gas into burned zone
        double mass_burned_old = mass_burned;
        mass_burned += dm_out;
        
        if (mass_burned > 1e-12) {
            // Mass-weighted mixing of compositions (reuse Y_f and Y_b from earlier)
            std::vector<double> Y_burned_new(Y_b.size());
            for (size_t i = 0; i < Y_b.size(); ++i) {
                Y_burned_new[i] = (Y_b[i] * mass_burned_old + Y_f[i] * dm_out) / mass_burned;
            }
            
            // Mass-weighted mixing of enthalpies
            double u_burned_mixed = (u_b_initial * mass_burned_old + u_f * dm_out) / mass_burned;
            
            // Set composition first, then state with enthalpy and pressure
            burned_zone->set_mass_fractions(Y_burned_new);
            burned_zone->set_state_UV(u_burned_mixed, (1/rho_b*mass_burned_old + 1/rho_f*dm_out)/mass_burned);
        }
    }

    // Transfer mass from unburned to flame zone
    if (dm_in > 1e-15) {
        mass_unburned -= dm_in;
        
        // Mix incoming unburned gas into flame zone
        double mass_flame_old = mass_flame;
        mass_flame += dm_in;
        
        if (mass_flame > 1e-12) {
            // Mass-weighted mixing of compositions (reuse Y_f and Y_u from earlier)
            std::vector<double> Y_flame_new(Y_f.size());
            for (size_t i = 0; i < Y_f.size(); ++i) {
                Y_flame_new[i] = (Y_f[i] * mass_flame_old + Y_u[i] * dm_in) / mass_flame;
            }
            
            // Mass-weighted mixing of internal energies
            double u_flame_mixed = (u_f_initial * mass_flame_old + u_u_initial * dm_in) / mass_flame;
            
            // Set composition first, then state with enthalpy and pressure
            flame_zone->set_mass_fractions(Y_flame_new);
            flame_zone->set_state_UV(u_flame_mixed, (1/rho_u*dm_in + 1/rho_f*mass_flame_old)/mass_flame);
        }
    }
    
    // Advance chemistry in all zones (reactions progress here)
    if (mass_flame > 1e-12) {
        flame_zone->advance_chemistry(dt);
    }
    if (mass_burned > 1e-12) {
        burned_zone->advance_chemistry(dt);
    }
    if (mass_unburned > 1e-12) {
        unburned_zone->advance_chemistry(dt);
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
    
    // Update flame radius based on burned mass
    double burned_volume = mass_burned / burned_zone->density();
    flame_radius = std::pow(3.0 * burned_volume / (4.0 * M_PI), 1.0/3.0);
    
    // Check if combustion is complete
    if (mass_unburned < 1e-12 && mass_flame < 1e-12) {
        combustion_active = false;
    }
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
    
    // Get mass fractions once for each zone
    std::vector<double> Y_u = unburned_zone->mass_fractions();
    std::vector<double> Y_f = flame_zone->mass_fractions();
    std::vector<double> Y_b = burned_zone->mass_fractions();
    
    // Calculate mass-weighted average composition
    size_t n_species = gas->n_species();
    std::vector<double> Y_avg(n_species, 0.0);
    
    for (size_t i = 0; i < n_species; ++i) {
        Y_avg[i] = (Y_u[i] * mass_unburned + Y_f[i] * mass_flame + Y_b[i] * mass_burned) / total_mass;
    }
    
    // Calculate mass-weighted average internal energy (cache individual values)
    double u_u = unburned_zone->internal_energy();
    double u_f = flame_zone->internal_energy();
    double u_b = burned_zone->internal_energy();
    double u_avg = (u_u * mass_unburned + u_f * mass_flame + u_b * mass_burned) / total_mass;
    
    // Calculate specific volume
    double v_specific = volume / total_mass;
    
    // Update composition to match mass-weighted average, then set state
    gas->set_mass_fractions(Y_avg);
    gas->set_state_UV(u_avg, v_specific);
}

void Cylinder::heat_loss(double dt, double mean_piston_speed, double motored_pressure, bool valves_flowing) {
    PROFILE_FUNCTION();
    
    // Get state after isentropic process
    double T_after_isentropic = gas->temperature();
    double P_after_isentropic = gas->pressure();
    
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
}

void Cylinder::add_mass_flow(double mass_flow_rate, double dt, const std::string& incoming_composition, double incoming_temperature) {
    double dm = mass_flow_rate * dt;
    if (dm <= 0.0) return;
    
    // Use persistent mixing gas object instead of creating new one
    double current_pressure = gas->pressure();
    mixing_gas->set_state_TPX(incoming_temperature, current_pressure, incoming_composition);
    
    // Cache values before mixing
    double current_mass = mass;
    double new_total_mass = current_mass + dm;
    double current_internal_energy = gas->internal_energy();
    double incoming_enthalpy = mixing_gas->enthalpy();
    
    // Get mass fractions once
    std::vector<double> Y_current = gas->mass_fractions();
    std::vector<double> Y_incoming = mixing_gas->mass_fractions();

    // Mass-weighted mixing of species compositions
    size_t n_species = Y_current.size();
    std::vector<double> Y_mixed(n_species);
    double sum_Y = 0.0;
    
    for (size_t i = 0; i < n_species; ++i) {
        Y_mixed[i] = (current_mass * Y_current[i] + dm * Y_incoming[i]) / new_total_mass;
        sum_Y += Y_mixed[i];
    }
    
    // Normalize mass fractions to ensure they sum to 1.0
    if (sum_Y > 0.0) {
        double inv_sum = 1.0 / sum_Y;
        for (double& yi : Y_mixed) yi *= inv_sum;
    }
    
    // Update composition first, then state
    gas->set_mass_fractions(Y_mixed);
    gas->set_state_UV((current_internal_energy*current_mass + incoming_enthalpy*dm)/new_total_mass, volume/new_total_mass);
    mass = new_total_mass;
    
    // Reset combustion state when new mass is added
    combustion_initiated = false;
    combustion_active = false;
}

void Cylinder::remove_mass_flow(double mass_flow_rate, double dt) {
    double dm = mass_flow_rate * dt;
    if (dm <= 0.0) return;
    
    // Store current state
    double current_mass = mass;
    double current_internal_energy = gas->internal_energy();  // Specific enthalpy [J/kg]
    double current_enthalpy = gas->enthalpy();  // Specific enthalpy [J/kg]
    auto current_Y = gas->mass_fractions();
    
    // Calculate new mass (prevent excessive mass removal)
    const double min_mass_fraction = 0.1;
    double new_mass = std::max(current_mass * min_mass_fraction, current_mass - dm);
    
    // Calculate enthalpy balance
    double total_internal_energy_before = current_mass * current_internal_energy;  // [J]
    double energy_removed = dm * current_enthalpy;  // [J]
    double total_energy_after = total_internal_energy_before - energy_removed;  // [J]
    
    // Calculate new specific energy of remaining gas
    double new_specific_energy = total_energy_after / new_mass;  // [J/kg]
    
    // Update mass
    mass = new_mass;
    
    // Reset combustion state when mass is removed (exhaust valve open)
    combustion_initiated = false;
    combustion_active = false;
    
    // Set state using enthalpy and pressure
    double current_pressure = gas->pressure();
    gas->set_state_UV(new_specific_energy, volume/mass);
        
}

void Cylinder::advance_chemistry(double dt) {
    // Update combustion zones if combustion is active
    if (combustion_active) {
        
        update_combustion_zones(dt);
        
        mix_zones_to_cylinder();
    } else if (!combustion_initiated) {
        // Only advance chemistry if combustion hasn't started yet
        gas->advance_chemistry(dt);
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
