#ifndef ENGINE_H                                                                  
#define ENGINE_H                                                                  
                                                                                  
#include "Cylinder.h"                                                             
#include "ValveFlow.h"                                                            
#include "FrictionModel.h"                                                        
#include "Crankshaft.h"
#include "SparkModel.h"                                                           
#include "CanteraInterface.h"
#include "ConfigLoader.h"
#include <memory>                                                                 
#include <vector>                                                                 
#include <fstream>                                                                
                                                                                  
struct EngineResults {                                                            
    std::vector<double> crank_angles;                                             
    std::vector<double> cylinder_pressures;                                       
    std::vector<double> cylinder_temperatures;                                    
    std::vector<double> cylinder_volumes;                                         
    std::vector<double> cylinder_masses;                                          
    std::vector<double> heat_release_rates;                                       
    std::vector<double> heat_transfer_rates;                                      
    std::vector<double> intake_mass_flows;                                        
    std::vector<double> exhaust_mass_flows;                                       
    std::vector<double> friction_powers;
    std::vector<double> gamma_values;
    std::vector<double> pumping_works;                    // Add this line
    std::vector<double> exhaust_energy_losses;           // Add this line
    std::vector<double> heat_transfer_coefficients;      // Woschni h_conv
    std::vector<double> surface_areas;                   // Cylinder surface area      
    std::vector<double> gas_temperatures;                // Gas temperature for Woschni
    std::vector<std::vector<double>> species_mass_fractions;
    std::vector<std::string> species_names;                                       
};                                               
                                                                                  
class Engine {                                                                    
private:                                                                          
    std::unique_ptr<Cylinder> cylinder;                                           
    std::unique_ptr<ValveFlow> intake_valve;                                      
    std::unique_ptr<ValveFlow> exhaust_valve;                                     
    std::unique_ptr<FrictionModel> friction;                                      
    std::unique_ptr<Crankshaft> crankshaft;                                       
    std::unique_ptr<SparkModel> spark;
    
    // Add persistent intake gas object
    std::unique_ptr<CanteraInterface> intake_gas;
                                                                                  
    // Engine parameters                                                          
    double bore;              // [m]                                              
    double stroke;            // [m]                                              
    double con_rod_length;    // [m]                                              
    double clearance_volume;  // [m³]                                             
    double compression_ratio;                                                     
    double engine_speed;      // [rpm]                                            
                                                                                  
    // Simulation parameters                                                      
    double crank_angle;       // [degrees]                                        
    double crank_angle_step;  // [degrees]                                        
    double time_step;         // [seconds]                                        
                                                                                  
    // Boundary conditions                                                        
    double intake_pressure;   // [Pa]                                             
    double intake_temperature; // [K]                                             
    double exhaust_pressure;  // [Pa]
    std::string intake_composition;  // Intake gas composition string
    std::string mechanism_file;  // Path to the mechanism file
    
    // For efficiency calculation
    double fuel_energy_rate;  // [W] - fuel energy input rate
    double lower_heating_value; // [J/kg] - LHV of fuel
    
    // Energy ledger tracking
    double total_lhv_energy;           // [J] - Total LHV energy available
    double indicated_work;             // [J] - Useful work from combustion
    double pumping_work;               // [J] - Work lost to pumping
    double friction_work;              // [J] - Work lost to friction
    double heat_transfer_loss;         // [J] - Energy lost to heat transfer
    double exhaust_energy_loss;        // [J] - Energy lost in exhaust
    double reference_enthalpy;         // [J/kg] - Reference enthalpy at intake conditions
    double T3;         // [] - Temperature just after combustion
                                                                                  
    // Results storage                                                            
    EngineResults results;
    
    // Store parameter structures
    SimulationParameters sim_params;
    EngineParameters engine_params;
                                                                                      
    void calculate_time_step();                                                   
    void update_cylinder_volume();                                                
    void calculate_valve_flows();                                                 
    void initialize_flame_model();                                           
    void store_results();                                             
    void calculate_energy_ledger();
                                                                                  
public:                                                                           
    Engine(double cylinder_bore, double cylinder_stroke, double                   
connecting_rod_length, double compression_ratio, const std::string& mechanism_file, const SimulationParameters& sim_params, const EngineParameters& engine_params);          
    ~Engine() = default;                                                          
                                                                                  
    void set_engine_speed(double rpm);                                            
    void set_boundary_conditions(double P_intake, double T_intake, double         
P_exhaust);                                                                       
    void set_intake_composition(const std::string& composition);
    void set_valve_profiles(const std::vector<double>& intake_angles,             
                           const std::vector<double>& intake_lifts,               
                           const std::vector<double>& exhaust_angles,             
                           const std::vector<double>& exhaust_lifts);             
                                                                                  
    void initialize_cycle();                                                      
    void run_cycle();                                                             
    void run_multiple_cycles(int num_cycles);                                     
                                                                                  
    void output_results(const std::string& filename = "engine_results.csv");      
    void output_performance_summary(const std::string& filename = "engine_performance.csv");
                                                                                  
    // Getters                                                                    
    double get_indicated_mean_effective_pressure() const;                         
    double get_brake_mean_effective_pressure() const;
    double get_gross_mean_effective_pressure() const;
    double get_pumping_mean_effective_pressure() const;
    double get_pumping_work() const;
    double get_indicated_power() const;                                           
    double get_brake_power() const;
    double get_cycle_efficiency() const;
    
    // Spark control methods
    void set_spark_timing(double timing_btdc);
    void set_spark_energy(double energy);
    void enable_spark(bool enabled);
    void print_energy_ledger() const;
};                                                                                
                                                                                  
#endif
