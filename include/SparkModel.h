#ifndef SPARK_MODEL_H
#define SPARK_MODEL_H

class SparkModel {
private:
    double spark_timing;        // Spark timing [degrees BTDC]
    double spark_duration;      // Duration of spark energy release [degrees]
    double spark_energy;        // Total spark energy [J]
    bool spark_enabled;         // Enable/disable spark
    
public:
    SparkModel(double timing_btdc = 30.0, double duration = 2.0, double energy = 0.03);
    
    void set_spark_timing(double timing_btdc);
    void set_spark_energy(double energy);
    void set_spark_duration(double duration);
    void enable_spark(bool enabled);
    
    double get_spark_power(double crank_angle) const;
    bool is_spark_active(double crank_angle) const;
};

#endif
