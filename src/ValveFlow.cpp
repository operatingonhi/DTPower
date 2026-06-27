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
    
    // Normalize angle to 0-720 degree range (4-stroke cycle)
    double normalized_angle = crank_angle;
    while (normalized_angle < 0.0) {
        normalized_angle += 720.0;
    }
    while (normalized_angle >= 720.0) {
        normalized_angle -= 720.0;
    }
    
    // Linear interpolation between lift profile points
    for (size_t i = 0; i < lift_angles.size() - 1; i++) {
        if (normalized_angle >= lift_angles[i] && normalized_angle <= lift_angles[i + 1]) {
            double angle_fraction = (normalized_angle - lift_angles[i]) / 
                                   (lift_angles[i + 1] - lift_angles[i]);
            return lift_values[i] + angle_fraction * (lift_values[i + 1] - lift_values[i]);
        }
    }
    
    return 0.0;  // Outside defined lift profile range
}

double ValveFlow::get_lift(double crank_angle) const {
    return interpolate_lift(crank_angle);
}

double ValveFlow::calculate_effective_area(double lift) const {
    if (lift <= 0.0) {
        return 0.0;
    }
    
    // Curtain area: flow area around the valve seat
    double curtain_area = M_PI * seat_diameter * lift;
    
    // Valve head area: maximum possible flow area
    double valve_head_area = M_PI * valve_diameter * valve_diameter / 4.0;
    
    // Effective area is limited by the smaller of the two
    return std::min(curtain_area, valve_head_area);
}

bool ValveFlow::is_choked(double pressure_ratio) const {
    // Critical pressure ratio for air (gamma = 1.4): approximately 0.528
    const double gamma_air = 1.4;
    double critical_ratio = std::pow(2.0 / (gamma_air + 1.0), gamma_air / (gamma_air - 1.0));
    return pressure_ratio <= critical_ratio;
}

double ValveFlow::calculate_mass_flow(double lift, double P_upstream, double P_downstream,
                                     double T_upstream, double rho_upstream, double gamma) const {
    if (lift <= 0.0) {
        return 0.0;
    }
    
    // Determine flow direction and actual upstream/downstream pressures
    double flow_direction = 1.0;
    double P_high = P_upstream;
    double P_low = P_downstream;
    double T_high = T_upstream;
    
    if (P_upstream < P_downstream) {
        // Flow is reversed (from downstream to upstream)
        P_high = P_downstream;
        P_low = P_upstream;
        flow_direction = -1.0;
    }
    
    double effective_area = calculate_effective_area(lift);
    if (effective_area <= 0.0) return 0.0;
    
    double pressure_ratio = P_low / P_high;
    
    // Specific gas constant for air [J/(kg·K)]
    const double R_air = 287.0;
    
    double mass_flow = 0.0;
    
    if (is_choked(pressure_ratio)) {
        // Choked flow (sonic conditions at throat)
        mass_flow = discharge_coeff * effective_area * P_high * 
                   std::sqrt(gamma / (R_air * T_high)) * 
                   std::pow(2.0 / (gamma + 1.0), (gamma + 1.0) / (2.0 * (gamma - 1.0)));
    } else {
        // Subsonic flow
        mass_flow = discharge_coeff * effective_area * P_high *
                   std::sqrt(gamma / (R_air * T_high)) *
                   std::sqrt((2.0 / (gamma - 1.0)) *
                            (std::pow(pressure_ratio, 2.0 / gamma) -
                             std::pow(pressure_ratio, (gamma + 1.0) / gamma)));
    }
    
    return mass_flow * flow_direction;
}
