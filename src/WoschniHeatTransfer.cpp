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
    // where w = C2 * mean_piston_speed + additional terms for combustion
    
    double w = C2 * mean_piston_speed;
    
    // Add combustion term if motored pressure is provided
    if (motored_pressure > 0.0 && pressure > motored_pressure) {
	double Tr = 300; //K
	double Pr = 1.0e5; //K
        w += C3 * (pressure - motored_pressure) * Tr / Pr; // Convert to m/s
    }
    
    // Woschni correlation with CORRECT units
    // Original Woschni: h = 3.26 * B^(-0.2) * P^0.8 * T^(-0.53) * w^0.8
    // But this gives h in W/(m²K) when:
    // B is in meters, P is in bar, T is in K, w is in m/s
    double h = C1 * std::pow(bore, -0.2) * std::pow(pressure / 1e5, 0.8) * 
               std::pow(temperature, -0.53) * std::pow(w, 0.8);
    //std::cout << " h=" << h << " W/(m2*s)" << std::endl;
    
    return h;
}

double WoschniHeatTransfer::calculate_heat_flux(double T_gas, double T_wall, 
                                               double surface_area, double h_conv) {
    double heat_flux = h_conv * surface_area * (T_gas - T_wall);
    
    // Debug output
    //std::cout << "Heat transfer debug: h=" << h_conv << " W/(m²K), A=" << surface_area << " m², ΔT=" << (T_gas - T_wall) << " K, Q=" << heat_flux/1000.0 << " kW" << std::endl;
    
    return heat_flux;
}

double WoschniHeatTransfer::calculate_surface_area(double volume, double bore, double stroke, double clearance_volume) {
    // Surface area = cylinder head + exposed cylinder walls + piston top
    double head_area = M_PI * bore * bore / 4.0;    // Cylinder head (constant)
    double piston_area = M_PI * bore * bore / 4.0;  // Piston top (constant)
    double clearance_height = clearance_volume / piston_area;
    
    // Calculate actual piston position from volume
    double piston_height = (volume - clearance_volume) / piston_area;
    piston_height = std::max(0.0, std::min(piston_height, stroke)); // Clamp to stroke limits
    
    // Exposed cylinder wall area = circumference × exposed height
    // At TDC: piston_height ≈ 0, exposed_height ≈ stroke (maximum wall area)
    // At BDC: piston_height ≈ stroke, exposed_height ≈ 0 (minimum wall area)
    double exposed_height = clearance_height + piston_height;
    double wall_area = M_PI * bore * exposed_height;
    
    return head_area + wall_area + piston_area;
}
