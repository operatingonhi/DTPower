#ifndef VALVE_FLOW_H                                                              
#define VALVE_FLOW_H                                                              
                                                                                  
#include <vector>                                                                 
#include <map>                                                                    
                                                                                  
class ValveFlow {                                                                 
private:                                                                          
    double valve_diameter;     // Valve head diameter [m]                         
    double seat_diameter;      // Valve seat diameter [m]                         
    double discharge_coeff;    // Discharge coefficient                           
                                                                                  
    std::vector<double> lift_angles;  // Crank angles for lift profile            
    std::vector<double> lift_values;  // Corresponding lift values [m]            
                                                                                  
    double interpolate_lift(double crank_angle) const;                            
    double calculate_effective_area(double lift) const;                           
    bool is_choked(double pressure_ratio) const;                                  
                                                                                  
public:                                                                           
    ValveFlow(double valve_dia, double seat_dia, double Cd = 0.6);                
                                                                                  
    void set_lift_profile(const std::vector<double>& crank_angles,                
                         const std::vector<double>& lifts);                       
                                                                                  
    double get_lift(double crank_angle) const;                                    
                                                                                  
    double calculate_mass_flow(double lift, double P_upstream, double             
P_downstream,                                                                     
                              double T_upstream, double rho_upstream,             
                              double gamma = 1.4) const;                          
}; 
                                                                                  
#endif
