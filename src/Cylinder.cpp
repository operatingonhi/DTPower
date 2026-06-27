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
      mechanism_file(mechanism_file) {
    
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

void Cylinder::initiate_spark_ignition(double spark_energy) {
    if (!combustion_initiated && spark_energy > 0.0) {
        // Initialize the three-zone combustion model
        initialize_combustion_zones();
        
        // Calculate adiabatic flame temperature for current conditions
        calculate_adiabatic_flame_temperature();
        
        // Set target flame temperature
        target_flame_temp = flame_temp_fraction * adiabatic_flame_temp;
        
        // Reset flame radius to initial kernel size
        flame_radius = initial_flame_radius;
        
        // Mark combustion as initiated and active
        combustion_initiated = true;
        combustion_active = true;
        
        std::cout << "Combustion initiated: Tad = " << adiabatic_flame_temp 
                  << " K, Target flame temp = " << target_flame_temp << " K" << std::endl;
    }
}

void Cylinder::initialize_combustion_zones() {
    // Initialize all zones with current cylinder state
    double T = gas->temperature();
    double P = gas->pressure();
    std::vector<double> Y = gas->mass_fractions();
    
    // Unburned zone starts with all the mass
    unburned_zone->set_state_TPY(T, P, Y);
    mass_unburned = mass;
    
    // Flame zone starts with a small initial kernel
    flame_zone->set_state_TPY(T, P, Y);
    mass_flame = 1e-9;  // Very small initial mass
    
    // Burned zone starts empty
    burned_zone->set_state_TPY(T, P, Y);
    mass_burned = 0.0;
}

void Cylinder::calculate_adiabatic_flame_temperature() {
    // Create a temporary gas object for adiabatic flame calculation
    auto temp_gas = std::make_unique<CanteraInterface>(mechanism_file);
    
    // Set to current unburned conditions
    double T = gas->temperature();
    double P = gas->pressure();
    std::vector<double> Y = gas->mass_fractions();
    temp_gas->set_state_TPY(T, P, Y);
    
    // Equilibrate at constant enthalpy and pressure to get adiabatic flame temperature
    try {
        temp_gas->get_solution()->thermo()->equilibrate("HP");
        adiabatic_flame_temp = temp_gas->temperature();
    } catch (const std::exception& e) {
        std::cerr << "Warning: Adiabatic flame temperature calculation failed: " 
                  << e.what() << std::endl;
        adiabatic_flame_temp = T * 3.0;  // Rough estimate
    }
}

double Cylinder::calculate_burning_rate() {
    // Calculate mass burning rate based on flame surface area and burning velocity
    // Surface area of spherical flame front
    double flame_area = 4.0 * M_PI * flame_radius * flame_radius;  // [m²]
    
    // Unburned gas density
    double rho_u = unburned_zone->density();  // [kg/m³]
    
    // Mass burning rate = density * velocity * area
    double burning_rate = rho_u * burning_velocity * flame_area;  // [kg/s]
    
    return burning_rate;
}

void Cylinder::update_combustion_zones(double dt) {
    if (!combustion_active || mass_unburned <= 1e-12) {
        combustion_active = false;
        instantaneous_heat_release_rate = 0.0;
        return;
    }
    
    // Get current pressure (same for all zones)
    double P = gas->pressure();
    
    // Update unburned zone to current pressure
    double T_u = unburned_zone->temperature();
    std::vector<double> Y_u = unburned_zone->mass_fractions();
    unburned_zone->set_state_TPY(T_u, P, Y_u);
    
    // Calculate burning rate
    double dm_burn = calculate_burning_rate() * dt;  // [kg]
    dm_burn = std::min(dm_burn, mass_unburned);
    
    if (dm_burn <= 0.0) {
        combustion_active = false;
        instantaneous_heat_release_rate = 0.0;
        return;
    }
    
    // Transfer mass from unburned to flame zone
    mass_unburned -= dm_burn;
    mass_flame += dm_burn;
    
    // Update flame zone composition (mix incoming unburned gas)
    if (mass_flame > 1e-12) {
        // Simple mixing: weighted average of compositions
        std::vector<double> Y_flame = flame_zone->mass_fractions();
        for (size_t i = 0; i < Y_flame.size(); ++i) {
            Y_flame[i] = (Y_flame[i] * (mass_flame - dm_burn) + Y_u[i] * dm_burn) / mass_flame;
        }
        flame_zone->set_state_TPY(flame_zone->temperature(), P, Y_flame);
    }
    
    // Equilibrate flame zone at constant HP to get products
    double h_before_flame = flame_zone->enthalpy() * mass_flame;
    try {
        flame_zone->get_solution()->thermo()->equilibrate("HP");
    } catch (const std::exception& e) {
        std::cerr << "Warning: Flame zone equilibration failed" << std::endl;
    }
    double h_after_flame = flame_zone->enthalpy() * mass_flame;
    
    // Calculate heat release from flame zone
    double heat_release = h_after_flame - h_before_flame;  // [J]
    instantaneous_heat_release_rate = heat_release / dt;   // [W]
    
    // Control flame zone temperature by transferring mass to burned zone
    double T_flame = flame_zone->temperature();
    
    // Calculate mass to transfer to maintain target temperature
    double dm_to_burned = 0.0;
    if (T_flame > target_flame_temp && mass_flame > 1e-12) {
        // Transfer enough mass to cool flame zone to target temperature
        // This is a simplified approach - transfer a fraction of flame mass
        double temp_ratio = (T_flame - target_flame_temp) / T_flame;
        dm_to_burned = mass_flame * temp_ratio * 0.5;  // Transfer 50% of excess
        dm_to_burned = std::min(dm_to_burned, mass_flame * 0.9);  // Keep some mass in flame
    }
    
    if (dm_to_burned > 1e-12) {
        // Transfer mass from flame to burned zone
        mass_flame -= dm_to_burned;
        mass_burned += dm_to_burned;
        
        // Update burned zone (mix with incoming products)
        if (mass_burned > 1e-12) {
            std::vector<double> Y_flame_products = flame_zone->mass_fractions();
            std::vector<double> Y_burned = burned_zone->mass_fractions();
            
            for (size_t i = 0; i < Y_burned.size(); ++i) {
                Y_burned[i] = (Y_burned[i] * (mass_burned - dm_to_burned) + 
                              Y_flame_products[i] * dm_to_burned) / mass_burned;
            }
            
            // Mix enthalpies to get burned zone temperature
            double h_burned_old = burned_zone->enthalpy() * (mass_burned - dm_to_burned);
            double h_flame_transfer = flame_zone->enthalpy() * dm_to_burned;
            double h_burned_new = (h_burned_old + h_flame_transfer) / mass_burned;
            
            burned_zone->set_state_TPY(burned_zone->temperature(), P, Y_burned);
            // Iterate to match enthalpy
            for (int iter = 0; iter < 5; ++iter) {
                double h_current = burned_zone->enthalpy();
                double T_current = burned_zone->temperature();
                double cp = burned_zone->cp();
                double T_new = T_current + (h_burned_new - h_current) / cp;
                burned_zone->set_state_TPY(T_new, P, Y_burned);
            }
        }
    }
    
    // Update flame radius based on burned volume
    double total_volume = volume;
    double burned_volume_fraction = mass_burned / mass;
    double burned_volume = total_volume * burned_volume_fraction;
    flame_radius = std::pow(3.0 * burned_volume / (4.0 * M_PI), 1.0/3.0);
    
    // Check if combustion is complete
    if (mass_unburned < 1e-12 && mass_flame < 1e-12) {
        combustion_active = false;
        std::cout << "Combustion complete: all mass burned" << std::endl;
    }
}

void Cylinder::mix_zones_to_cylinder() {
    // Combine all three zones back into the main cylinder gas object
    // This maintains the overall cylinder state
    
    if (!combustion_initiated) {
        return;  // No combustion zones to mix
    }
    
    double P = gas->pressure();  // Maintain current pressure
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
    
    // Calculate mass-weighted average enthalpy
    double h_avg = (unburned_zone->enthalpy() * mass_unburned + 
                    flame_zone->enthalpy() * mass_flame + 
                    burned_zone->enthalpy() * mass_burned) / total_mass;
    
    // Set cylinder state with averaged composition and enthalpy
    gas->set_state_TPY(gas->temperature(), P, Y_avg);
    
    // Iterate to match enthalpy
    for (int iter = 0; iter < 5; ++iter) {
        double h_current = gas->enthalpy();
        double T_current = gas->temperature();
        double cp = gas->cp();
        double T_new = T_current + (h_avg - h_current) / cp;
        gas->set_state_TPY(T_new, P, Y_avg);
    }
}

void Cylinder::update_thermodynamics(double dt, double mean_piston_speed, double motored_pressure) {
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
    
    // Get state after chemistry (chemistry is now handled by the reactor)
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
