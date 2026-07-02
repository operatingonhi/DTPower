#ifndef CANTERA_INTERFACE_H                                                       
#define CANTERA_INTERFACE_H                                                       
                                                                                  
#include "cantera/core.h"                                                       
#include "cantera/zerodim.h"
#include <memory>                                                                 
#include <vector>                                                                 
#include <string>                                                                 
                                                                                  
class CanteraInterface {                                                          
private:                                                                          
    std::shared_ptr<Cantera::Solution> gas;                                       
    Cantera::ThermoPhase* thermo;                                                 
    Cantera::Kinetics* kinetics;                                                  
    std::unique_ptr<Cantera::IdealGasReactor> reactor;
    std::unique_ptr<Cantera::ReactorNet> reactor_net;
                                                                                  
public:                                                                           
    CanteraInterface(const std::string& mechanism_file);                          
    ~CanteraInterface() = default;                                                
                                                                                  
    void set_state_TPY(double T, double P, const std::vector<double>& Y);         
    void set_state_TPX(double T, double P, const std::string& X);                 
    void set_state_TP(double T, double P);                 
    void set_state_DP(double D, double P);                 
    void advance_chemistry(double dt);  

    // Set only composition (mass fractions or mole fractions)                                                              
    void set_mass_fractions(const std::vector<double>& Y);                                                                  
    void set_mole_fractions(const std::string& X);
                                                                                      
    // Property getters                                                           
    double temperature() const;                                                   
    double pressure() const;                                                      
    double density() const;                                                       
    double cp() const;                                                            
    double cv() const;                                                            
    double enthalpy() const;                                                      
    double internal_energy() const;                                               
    std::vector<double> mass_fractions() const;                                   
    double mean_molecular_weight() const;                                         
    double entropy() const;
    double gamma() const;  // Specific heat ratio (cp/cv)
                                                                                  
    // State setters using different property pairs
    void set_state_SV(double s, double v);
    void set_state_SD(double s, double rho);
    void set_state_UV(double u, double v);
    void set_state_HP(double h, double p);
                                                                                  
    // Kinetics                                                                   
    std::vector<double> net_production_rates() const;                             
    double heat_release_rate() const;                                             
    
    // Species information
    std::vector<std::string> species_names() const;
    size_t n_species() const;
    
    // Access to underlying Solution object
    std::shared_ptr<Cantera::Solution> get_solution() { return gas; }
};                                                                                
                                                                                  
#endif
