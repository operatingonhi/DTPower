#include "WoschniHeatTransfer.h"
#include <cmath>

WoschniHeatTransfer::WoschniHeatTransfer(double cylinder_bore, double cylinder_stroke,
                                       double c1, double c2, double c3)
    : C1(c1), C2(c2), C3(c3), bore(cylinder_bore), stroke(cylinder_stroke) {
}

double WoschniHeatTransfer::calculate_heat_transfer_coefficient(double pressure, 
                                                              double temperature,
                                                              double mean_piston_speed,
                                                              double motored_pressure) {
    // Woschni correlation: h = C1 * B^(-0.2) * P^0.8 * T^(-0.53) * w^0.8
    // where w = C2 * Sp + C3 * (Vd*Tr)/(Pr*Vr) * (P - Pm)
    
    double w = C2 * mean_piston_speed;
    
    // Add combustion term if motored pressure is provided
    // Simplified form: w += C3 * Vd/Vr * Tr/Pr * (P - Pm)
    // Using Vd/Vr ≈ 1 and Tr/Pr ≈ 300K/1e5Pa = 0.003 K/Pa
    if (motored_pressure > 0.0 && pressure > motored_pressure) {
        double Vd = M_PI * bore * bore * stroke / 4.0;  // Displacement volume
        double Tr = 300.0;  // Reference temperature [K] (typical intake)
        double Pr = 1.0e5;  // Reference pressure [Pa] (typical intake)
        double Vr = Vd;     // Reference volume (displacement volume)
        
        w += C3 * (Vd * Tr) / (Pr * Vr) * (pressure - motored_pressure);
    }
    
    // Woschni correlation
    // h in W/(m²K) when: B in meters, P in bar, T in K, w in m/s
    double h = C1 * std::pow(bore, -0.2) * std::pow(pressure / 1e5, 0.8) * 
               std::pow(temperature, -0.53) * std::pow(w, 0.8);
    
    return h;
}

double WoschniHeatTransfer::calculate_heat_flux(double T_gas, double T_wall, 
                                               double surface_area, double h_conv) {
    return h_conv * surface_area * (T_gas - T_wall);
}

double WoschniHeatTransfer::calculate_surface_area(double volume, double bore, double stroke, double clearance_volume) {
    // Surface area = cylinder head + exposed cylinder walls + piston top
    double piston_area = M_PI * bore * bore / 4.0;
    double clearance_height = clearance_volume / piston_area;
    
    // Calculate actual piston position from volume
    double piston_height = (volume - clearance_volume) / piston_area;
    piston_height = std::max(0.0, std::min(piston_height, stroke));
    
    // Exposed cylinder wall area = circumference × exposed height
    double exposed_height = clearance_height + piston_height;
    double wall_area = M_PI * bore * exposed_height;
    
    // Total surface area: head + walls + piston top
    return 2.0 * piston_area + wall_area;
}
