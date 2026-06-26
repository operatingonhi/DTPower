#ifndef FRICTION_MODEL_H                                                          
#define FRICTION_MODEL_H                                                          
                                                                                  
class FrictionModel {                                                             
private:                                                                          
    double bore;           // Cylinder bore [m]                                   
    double stroke;         // Cylinder stroke [m]                                 
    double con_rod_length; // Connecting rod length [m]                           
                                                                                  
    // Chen-Flynn model parameters                                                
    double cf_c1, cf_c2, cf_c3, cf_c4;                                            

        // For FMEP integration
    double accumulated_friction_work;  // [J] - Total friction work over cycle
    double cycle_start_angle;          // [deg] - Starting angle for integration
    bool integration_active;           // Flag to track if we're integrating
                                                                                  
public:                                                                           
    FrictionModel(double cylinder_bore, double cylinder_stroke, double con_rod_len, double c1 = 0.97, double c2 = 0.15, double c3 = 0.05, double c4 = 0.0);                                                                     
                                                                                  
    // Piston assembly friction (Chen-Flynn model)                                
    double calculate_piston_friction(double engine_speed, double                  
cylinder_pressure,                                                                
                                   double crank_angle);                           
                                                                                  
    // Journal bearing friction                                                   
    double calculate_bearing_friction(double engine_speed, double load);          
                                                                                  
    // Total friction mean effective pressure                                     
    double calculate_fmep(double engine_speed, double mean_cylinder_pressure);    
                                                                  
        // Instantaneous friction force                                               
    double calculate_friction_force(double engine_speed, double cylinder_pressure,
                                  double crank_angle);                            
         
    // New methods for time-varying friction
    double calculate_instantaneous_friction_power(double engine_speed, double 
cylinder_pressure,
                                                 double crank_angle, double 
piston_velocity);
     
    void start_friction_integration(double start_angle);
    void accumulate_friction_work(double friction_power, double dt);
    double get_integrated_fmep(double displacement);
    void reset_friction_integration();

};                                                                                
                                                                                  
#endif
