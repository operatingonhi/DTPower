#include "ValveFlow.h"
#include <cmath>
#include <algorithm>

ValveFlow::ValveFlow(double valve_dia, double seat_dia, double Cd)
    : valve_diameter(valve_dia), seat_diameter(seat_dia), discharge_coeff(Cd) {
}

void ValveFlow::set_lift_profile(const std::vector<double>& crank_angles,
                                const std::vector<double>& lifts) {
    lift_angles = crank_angles;
    lift_values = lifts;
}

double ValveFlow::interpolate_lift(double crank_angle) const {
    if (lift_angles.empty() || lift_values.empty()) {
        return 0.0;
    }
    
    // Handle negative angles by finding equivalent positive angle in the 0-720 range
    double normalized_angle = crank_angle;
    while (normalized_angle < 0) normalized_angle += 720.0;
    while (normalized_angle >= 720.0) normalized_angle -= 720.0;
    
    // Find interpolation points
    for (size_t i = 0; i < lift_angles.size() - 1; i++) {
        if (normalized_angle >= lift_angles[i] && normalized_angle <= lift_angles[i + 1]) {
            double t = (normalized_angle - lift_angles[i]) / (lift_angles[i + 1] - lift_angles[i]);
            return lift_values[i] + t * (lift_values[i + 1] - lift_values[i]);
        }
    }
    
    return 0.0; // Default if not found
}

double ValveFlow::get_lift(double crank_angle) const {
    return interpolate_lift(crank_angle);
}

double ValveFlow::calculate_effective_area(double lift) const {
    if (lift <= 0.0) return 0.0;
    
    // Effective area calculation
    double curtain_area = M_PI * seat_diameter * lift;
    double valve_area = M_PI * valve_diameter * valve_diameter / 4.0;
    
    // Use minimum of curtain area and valve area
    return std::min(curtain_area, valve_area);
}

bool ValveFlow::is_choked(double pressure_ratio) const {
    // Critical pressure ratio for air (gamma = 1.4)
    double critical_ratio = std::pow(2.0 / (1.4 + 1.0), 1.4 / (1.4 - 1.0));
    return pressure_ratio <= critical_ratio;
}

double ValveFlow::calculate_mass_flow(double lift, double P_upstream, double P_downstream,
                                     double T_upstream, double rho_upstream, double gamma) const {
    if (lift <= 0.0) {
        return 0.0;
    }
    double forward = 1.0;
    if (P_upstream <= P_downstream) {
        double holder = P_downstream;
	P_downstream = P_upstream;
	P_upstream = holder;
	forward = -1.0;
    }
    
    double effective_area = calculate_effective_area(lift);
    if (effective_area <= 0.0) return 0.0;
    
    double pressure_ratio = P_downstream / P_upstream;
    double mass_flow = 0.0;
    
    if (is_choked(pressure_ratio)) {
        // Choked flow
        double critical_ratio = std::pow(2.0 / (gamma + 1.0), gamma / (gamma - 1.0));
        mass_flow = discharge_coeff * effective_area * P_upstream * 
                   std::sqrt(gamma / (287.0 * T_upstream)) * 
                   std::pow(2.0 / (gamma + 1.0), (gamma + 1.0) / (2.0 * (gamma - 1.0)));
    } else {
        // Subsonic flow
        //mass_flow = discharge_coeff * effective_area * P_upstream * 
        //           std::sqrt(2.0 * gamma / ((gamma - 1.0) * 287.0 * T_upstream)) *
        //           std::sqrt(std::pow(pressure_ratio, 2.0/gamma) - std::pow(pressure_ratio, (gamma + 1.0)/gamma));
	mass_flow = discharge_coeff * effective_area * P_upstream *
            std::sqrt(gamma / (287.0 * T_upstream)) *
            std::sqrt( (2.0 / (gamma - 1.0)) *
                       (std::pow(pressure_ratio, 2.0/gamma) -
                        std::pow(pressure_ratio, (gamma + 1.0)/gamma)) );
    }
    
    return mass_flow*forward;
}
