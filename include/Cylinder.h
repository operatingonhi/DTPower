#ifndef CYLINDER_H                                                                
#define CYLINDER_H                                                                
                                                                                  
#include "CanteraInterface.h"                                                     
#include "WoschniHeatTransfer.h"                                                  
#include <memory>                                                                 
                                                                                  
class Cylinder {                                                                  
private:                                                                          
    std::unique_ptr<CanteraInterface> gas;                                        
    std::unique_ptr<WoschniHeatTransfer> heat_transfer;                           

    // Add persistent incoming gas object for mixing calculations
    std::unique_ptr<CanteraInterface> mixing_gas;
    
    // Persistent gas object for adiabatic flame temperature calculations
    std::unique_ptr<CanteraInterface> adiabatic_gas;
    
    // Three-zone combustion model
    std::unique_ptr<CanteraInterface> unburned_zone;    // Reactants
    std::unique_ptr<CanteraInterface> flame_zone;       // Flame front
    std::unique_ptr<CanteraInterface> burned_zone;      // Products
    
    double mass_unburned;      // Mass in unburned zone [kg]
    double mass_flame;         // Mass in flame zone [kg]
    double mass_burned;        // Mass in burned zone [kg]
    
    bool combustion_active;    // True when combustion is ongoing
    bool combustion_initiated; // True after spark has occurred
    double adiabatic_flame_temp; // Calculated adiabatic flame temperature [K]
    double target_flame_temp;    // Target flame zone temperature (0.75 * Tad) [K]
    
    // Combustion model parameters
    const double flame_temp_fraction = 0.1;  // Fraction of adiabatic flame temp
    const double flame_shell_thickness = 1.0e-3; // Flame shell thickness [m]
    double flame_radius;                      // Current flame radius [m]
    
    std::string mechanism_file;  // Store mechanism file path for temp gas objects
    
    // Helper methods for combustion model
    void initialize_combustion_zones();
    void calculate_adiabatic_flame_temperature(CanteraInterface* zone);
    void update_combustion_zones(double dt);
    double calculate_mass_flow_into_flame(double dt);
    double calculate_mass_flow_out_of_flame(double mass_flow_in);
    void mix_zones_to_cylinder();
                                                                                  
    double bore;              // Cylinder bore [m]                                
    double stroke;            // Cylinder stroke [m]                              
    double clearance_volume;  // Clearance volume [m³]                            
    double wall_temperature;  // Cylinder wall temperature [K]                    
    double mass;              // Gas mass in cylinder [kg]                        
    double volume;            // Current cylinder volume [m³]                     
    double previous_volume;   // Volume from previous time step [m³]
    double current_heat_transfer_rate;  // Current heat transfer rate [W]
    double instantaneous_heat_release_rate;  // Heat release rate from instantaneous combustion [W]
    
    // Store previous state for combustion zone updates
    double previous_u;  // Previous specific internal energy [J/kg]
    double previous_P;  // Previous pressure [Pa]
    
    // CSV logging for combustion diagnostics
    std::ofstream combustion_log;
    bool logging_enabled;
                                                                                  
public:                                                                           
    Cylinder(double cylinder_bore, double cylinder_stroke, double clearance_vol,  
             const std::string& mechanism_file, double wall_temp = 450.0);        
                                                                                  
    ~Cylinder() = default;
    
    void enable_combustion_logging(const std::string& filename);
    void disable_combustion_logging();
                                                                                  
    void set_volume(double new_volume);                                           
    void set_wall_temperature(double T_wall);
    void set_mass(double new_mass);
    void update_heat_transfer_parameters(double c1, double c2, double c3);
    
    // Set reference values for Woschni heat transfer (call at IVC during motored run)
    void set_heat_transfer_reference_values(double T_ivc, double P_ivc, double V_ivc) {
        heat_transfer->set_reference_values(T_ivc, P_ivc, V_ivc);
    }
    
    void initiate_spark_ignition(double spark_energy);
                                                                                     
    void heat_loss(double dt, double mean_piston_speed, double motored_pressure = 0.0, bool valves_flowing = false);
                                                                                  
    void add_mass_flow(double mass_flow_rate, double dt, const std::string& incoming_composition, double incoming_temperature);     
    void remove_mass_flow(double mass_flow_rate, double dt);                      
                                                                                  
    void advance_chemistry(double dt);                                            
                                                                                  
    // Getters                                                                    
    double get_temperature() const;                                               
    double get_pressure() const;                                                  
    double get_density() const;                                                   
    double get_mass() const;                                                      
    double get_volume() const;                                                    
    double get_internal_energy() const;                                           
    double get_enthalpy() const;                                                  
                                                                                  
    // Heat transfer                                                              
    double get_heat_transfer_rate() const;                                        
    double get_heat_release_rate();
    double get_gamma() const;

    WoschniHeatTransfer* get_heat_transfer() const { return heat_transfer.get(); }
                                                                                  
    CanteraInterface* get_gas() { return gas.get(); }                             
};                                                                                
                                                                                  
#endif
