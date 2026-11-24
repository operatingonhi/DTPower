#ifndef WOSCHNI_HEAT_TRANSFER_H                                                   
#define WOSCHNI_HEAT_TRANSFER_H                                                   
#include <iostream>
                                                                                  
class WoschniHeatTransfer {                                                       
private:                                                                          
    double C1;         // Woschni constant (typically 3.26)                       
    double C2;         // Woschni constant (typically 2.28)                    
    double C3;         // Woschni constant (typically 3.24e-3)                    
    double bore;       // Cylinder bore [m]                                       
    double stroke;     // Cylinder stroke [m]
    
    // Reference values at IVC for combustion term
    double Tr;         // Reference temperature at IVC [K]
    double Pr;         // Reference pressure at IVC [kPa]
    double Vr;         // Reference volume at IVC [m³]
    bool reference_values_set;  // Flag to track if reference values have been set
                                                                                  
public:                                                                           
    WoschniHeatTransfer(double cylinder_bore, double cylinder_stroke, double c1 = 3.26, double c2 = 2.28, double c3 = 3.24e-3);
    
    // Set reference values at IVC (should be called during motored run)
    void set_reference_values(double T_ivc, double P_ivc, double V_ivc);
    
    double calculate_heat_transfer_coefficient(double pressure, double            
temperature,                                                                      
    double mean_piston_speed,            
    double motored_pressure = 0.0);      
                                                                                  
    double calculate_heat_flux(double T_gas, double T_wall, double surface_area, double h_conv);                                     
                                                                                  
    double calculate_surface_area(double volume, double bore, double stroke, double clearance_volume);     
};                                                                                
                                                                                  
#endif    
