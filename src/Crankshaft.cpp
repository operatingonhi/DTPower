#include "Crankshaft.h"
#include <cmath>

Crankshaft::Crankshaft(double stroke_length, double connecting_rod_length, double cylinder_bore)
    : stroke(stroke_length), con_rod_length(connecting_rod_length), bore(cylinder_bore) {
    crank_radius = stroke / 2.0;
}

double Crankshaft::calculate_piston_position(double crank_angle) const {
    double theta = crank_angle * M_PI / 180.0; // Convert to radians
    double lambda = crank_radius / con_rod_length;
    
    // Piston position from TDC
    double position = crank_radius * (1.0 - std::cos(theta)) + 
                     con_rod_length * (1.0 - std::sqrt(1.0 - lambda * lambda * std::sin(theta) * std::sin(theta)));
    
    return position;
}

double Crankshaft::calculate_cylinder_volume(double crank_angle, double clearance_volume) const {
    double piston_position = calculate_piston_position(crank_angle);
    double bore_area = M_PI * bore * bore / 4.0; // Correct bore area calculation
    
    return clearance_volume + bore_area * piston_position;
}

double Crankshaft::calculate_piston_velocity(double crank_angle, double engine_speed) const {
    double theta = crank_angle * M_PI / 180.0;
    double omega = engine_speed * M_PI / 30.0; // Convert RPM to rad/s
    double lambda = crank_radius / con_rod_length;
    
    double velocity = omega * crank_radius * (std::sin(theta) + 
                     lambda * std::sin(2.0 * theta) / (2.0 * std::sqrt(1.0 - lambda * lambda * std::sin(theta) * std::sin(theta))));
    
    return velocity;
}

double Crankshaft::calculate_piston_acceleration(double crank_angle, double engine_speed) const {
    double theta = crank_angle * M_PI / 180.0;
    double omega = engine_speed * M_PI / 30.0;
    double lambda = crank_radius / con_rod_length;
    
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);
    double sqrt_term = std::sqrt(1.0 - lambda * lambda * sin_theta * sin_theta);
    
    double acceleration = omega * omega * crank_radius * 
                         (cos_theta + lambda * (std::cos(2.0 * theta) / sqrt_term + 
                          lambda * lambda * sin_theta * sin_theta * cos_theta / (sqrt_term * sqrt_term * sqrt_term)));
    
    return acceleration;
}

double Crankshaft::calculate_mean_piston_speed(double engine_speed) const {
    return 2.0 * stroke * engine_speed / 60.0; // m/s
}

double Crankshaft::calculate_volume_derivative(double crank_angle, double engine_speed) const {
    double piston_velocity = calculate_piston_velocity(crank_angle, engine_speed);
    double bore_area = M_PI * bore * bore / 4.0; // Correct bore area calculation
    
    return bore_area * piston_velocity;
}
