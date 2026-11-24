#include "WoschniHeatTransfer.h"
#include <cmath>

WoschniHeatTransfer::WoschniHeatTransfer(double cylinder_bore, double cylinder_stroke,
                                       double c1, double c2, double c3)
    : C1(c1), C2(c2), C3(c3), bore(cylinder_bore), stroke(cylinder_stroke),
      Tr(300.0), Pr(100.0), Vr(0.0), reference_values_set(false) {
    // Initialize with default values; should be updated with actual IVC values
    double Vd = M_PI * bore * bore * stroke / 4.0;  // Displacement volume
    Vr = Vd;  // Default reference volume
}

void WoschniHeatTransfer::set_reference_values(double T_ivc, double P_ivc, double V_ivc) {
    Tr = T_ivc;
    Pr = P_ivc / 1.0e3;  // Convert Pa to kPa
    Vr = V_ivc;
    reference_values_set = true;
}

double WoschniHeatTransfer::calculate_heat_transfer_coefficient(double pressure, 
                                                              double temperature,
                                                              double mean_piston_speed,
                                                              double motored_pressure) {
    // Woschni correlation: h = C1 * B^(-0.2) * P^0.8 * T^(-0.53) * w^0.8
    // where w = C2 * Sp + C3 * (Vd*Tr)/(Pr*Vr) * (P - Pm)
    
    double w = C2 * mean_piston_speed;

    double P_kPa = pressure / 1.0e3;   // Pa -> kPa
    double Pm_kPa = motored_pressure / 1.0e3;  // Pa -> kPa
    
    // Add combustion term if motored pressure is provided and reference values are set
    if (motored_pressure > 0.0 && pressure > motored_pressure && reference_values_set) {
        double Vd = M_PI * bore * bore * stroke / 4.0;  // Displacement volume [m³]
        
        // Use reference values from IVC
        w += C3 * (Vd * Tr) / (Pr * Vr) * (P_kPa - Pm_kPa);
    }
    
    // Woschni correlation
    // h in W/(m²K) when: B in meters, P in kPa, T in K, w in m/s
    double h = C1 * std::pow(bore, -0.2) * std::pow(P_kPa, 0.8) * 
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
