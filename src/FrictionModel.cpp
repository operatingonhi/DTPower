#include "FrictionModel.h"
#include <cmath>

FrictionModel::FrictionModel(double cylinder_bore, double cylinder_stroke, double con_rod_len,
                           double c1, double c2, double c3, double c4)
    : bore(cylinder_bore), stroke(cylinder_stroke), con_rod_length(con_rod_len),
      cf_c1(c1), cf_c2(c2), cf_c3(c3), cf_c4(c4),
      accumulated_friction_work(0.0), cycle_start_angle(0.0), integration_active(false) {
    // Parameters are now passed in from YAML file
}

double FrictionModel::calculate_piston_friction(double engine_speed, double cylinder_pressure,
                                               double crank_angle) {
    // Convert RPM to rad/s
    double omega = engine_speed * M_PI / 30.0;
    
    // Mean piston speed [m/s]
    double mean_piston_speed = 2.0 * stroke * engine_speed / 60.0;
    
    // Instantaneous piston velocity
    double crank_angle_rad = crank_angle * M_PI / 180.0;
    double lambda = stroke / (2.0 * con_rod_length); // Connecting rod ratio
    
    double piston_velocity = omega * stroke / 2.0 * 
                           (std::sin(crank_angle_rad) + lambda * std::sin(2.0 * crank_angle_rad));
    
    // Friction force based on Chen-Flynn model
    double friction_pressure = cylinder_pressure*(cf_c1 + cf_c2 * mean_piston_speed + cf_c3 * mean_piston_speed * mean_piston_speed);
    
    // Convert to friction force
    double piston_area = M_PI * bore * bore / 4.0;
    double friction_force = friction_pressure * piston_area; // Convert bar to Pa
    
    // Apply sign based on piston velocity direction
    if (piston_velocity != 0.0) {
        friction_force *= (piston_velocity > 0.0) ? 1.0 : -1.0;
    }
    
    return friction_force;
}

double FrictionModel::calculate_bearing_friction(double engine_speed, double load) {
    // Simple bearing friction model
    double bearing_friction_coeff = 0.01; // Typical value
    return bearing_friction_coeff * load * engine_speed / 1000.0; // Simplified
}

double FrictionModel::calculate_fmep(double engine_speed, double mean_cylinder_pressure) {
    // Friction Mean Effective Pressure using Chen-Flynn correlation
    double mean_piston_speed = 2.0 * stroke * engine_speed / 60.0;
    
    double fmep = cf_c1 + cf_c2 * mean_piston_speed + cf_c3 * mean_piston_speed * mean_piston_speed;
    
    // Add pressure-dependent term
    fmep += cf_c4 * mean_cylinder_pressure / 1e5; // Convert Pa to bar
    
    return fmep * 1e5; // Convert bar to Pa
}

double FrictionModel::calculate_friction_force(double engine_speed, double cylinder_pressure,
                                              double crank_angle) {
    return calculate_piston_friction(engine_speed, cylinder_pressure, crank_angle);
}

double FrictionModel::calculate_instantaneous_friction_power(double engine_speed, double cylinder_pressure,
                                                           double crank_angle, double piston_velocity) {
    // Calculate instantaneous friction force
    double friction_force = calculate_friction_force(engine_speed, cylinder_pressure, crank_angle);
    
    // Friction power = friction force × piston velocity
    double friction_power = std::abs(friction_force * piston_velocity);
    
    return friction_power;
}

void FrictionModel::start_friction_integration(double start_angle) {
    cycle_start_angle = start_angle;
    accumulated_friction_work = 0.0;
    integration_active = true;
}

void FrictionModel::accumulate_friction_work(double friction_power, double dt) {
    if (integration_active) {
        accumulated_friction_work += friction_power * dt;
    }
}

double FrictionModel::get_integrated_fmep(double displacement) {
    if (!integration_active || displacement <= 0.0) {
        return 0.0;
    }
    
    // FMEP = Total friction work per cycle / displacement volume
    return accumulated_friction_work / displacement;
}

void FrictionModel::reset_friction_integration() {
    accumulated_friction_work = 0.0;
    integration_active = false;
}
