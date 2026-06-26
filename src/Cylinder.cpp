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
      previous_volume(clearance_vol), combustion_occurred(false), spark_conditions_met(false),
      current_heat_transfer_rate(0.0), instantaneous_heat_release_rate(0.0) {
    
    gas = std::make_unique<CanteraInterface>(mechanism_file);
    heat_transfer = std::make_unique<WoschniHeatTransfer>(bore, stroke);
    
    // Initialize persistent mixing gas object
    mixing_gas = std::make_unique<CanteraInterface>(mechanism_file);
    
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

void Cylinder::handle_spark_ignition(double spark_energy) {
    if (!combustion_occurred && spark_energy > 0.0) {
        // Store state before combustion for heat release calculation
        double h_before = gas->enthalpy() * mass;  // Total enthalpy before [J]
        
        // Perform equilibrium combustion at constant internal energy and volume
        gas->get_solution()->thermo()->equilibrate("UV");
        
        // Calculate heat release from enthalpy change
        double h_after = gas->enthalpy() * mass;   // Total enthalpy after [J]
        double heat_released = h_before - h_after; // Heat released (positive value) [J]
        
        // Store the instantaneous heat release rate for this time step
        // Assume combustion occurs over 1 microsecond
        const double combustion_duration = 1e-6;  // [s]
        instantaneous_heat_release_rate = heat_released / combustion_duration;  // [J/s]
        
        std::cout << "Instantaneous combustion completed!" << std::endl;
        std::cout << "  Temperature after: " << gas->temperature() << " K" << std::endl;
        std::cout << "  Pressure after: " << gas->pressure() / 1e5 << " bar" << std::endl;
        std::cout << "  Heat released: " << heat_released / 1e6 << " MJ" << std::endl;
        std::cout << "  Heat release rate: " << instantaneous_heat_release_rate / 1e6 << " MW" << std::endl;
        
        combustion_occurred = true;
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
        DEBUG_LOG("Combustion occurred: " + std::to_string(combustion_occurred));
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
        combustion_occurred = false;
        
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
        
        // Reset combustion state
        combustion_occurred = false;
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
    if (!combustion_occurred) {
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
