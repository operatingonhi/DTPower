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
    
    // Simple combustion state
    bool combustion_occurred;
    bool spark_conditions_met;
                                                                                  
    double bore;              // Cylinder bore [m]                                
    double stroke;            // Cylinder stroke [m]                              
    double clearance_volume;  // Clearance volume [m³]                            
    double wall_temperature;  // Cylinder wall temperature [K]                    
    double mass;              // Gas mass in cylinder [kg]                        
    double volume;            // Current cylinder volume [m³]                     
    double previous_volume;   // Volume from previous time step [m³]
    double current_heat_transfer_rate;  // Current heat transfer rate [W]
    double instantaneous_heat_release_rate;  // Heat release rate from instantaneous combustion [W]
                                                                                  
public:                                                                           
    Cylinder(double cylinder_bore, double cylinder_stroke, double clearance_vol,  
             const std::string& mechanism_file, double wall_temp = 450.0);        
                                                                                  
    ~Cylinder() = default;                                                        
                                                                                  
    void set_volume(double new_volume);                                           
    void set_wall_temperature(double T_wall);
    void set_mass(double new_mass);
    void update_heat_transfer_parameters(double c1, double c2, double c3);
    
    void handle_spark_ignition(double spark_energy);
                                                                                     
    void update_thermodynamics(double dt, double mean_piston_speed, double motored_pressure = 0.0);
                                                                                  
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
