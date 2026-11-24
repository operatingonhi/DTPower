#ifndef CRANKSHAFT_H                                                              
#define CRANKSHAFT_H                                                              
                                                                                  
class Crankshaft {                                                                
private:                                                                          
    double stroke;         // Stroke length [m]                                   
    double con_rod_length; // Connecting rod length [m]                           
    double crank_radius;   // Crank radius [m]                                    
    double bore;   // Crank radius [m]                                    
                                                                                  
public:                                                                           
    Crankshaft(double stroke_length, double connecting_rod_length, double cylinder_bore);               
                                                                                  
    // Kinematic calculations                                                     
    double calculate_piston_position(double crank_angle) const;                   
    double calculate_cylinder_volume(double crank_angle, double clearance_volume) const;                                                                            
    double calculate_piston_velocity(double crank_angle, double engine_speed) const;                                                                            
    double calculate_piston_acceleration(double crank_angle, double engine_speed) const;                                                                            
                                                                                  
    // Mean piston speed                                                          
    double calculate_mean_piston_speed(double engine_speed) const;                
                                                                                  
    // Volume derivative                                                          
    double calculate_volume_derivative(double crank_angle, double engine_speed) const;                                                                            
};                                                                                
                                                                                  
#endif     
