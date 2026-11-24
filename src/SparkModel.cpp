#include "SparkModel.h"
#include <cmath>

SparkModel::SparkModel(double timing_btdc, double duration, double energy)
    : spark_timing(timing_btdc), spark_duration(duration), 
      spark_energy(energy), spark_enabled(true) {
}

void SparkModel::set_spark_timing(double timing_btdc) {
    spark_timing = timing_btdc;
}

void SparkModel::set_spark_energy(double energy) {
    spark_energy = energy;
}

void SparkModel::set_spark_duration(double duration) {
    spark_duration = duration;
}

void SparkModel::enable_spark(bool enabled) {
    spark_enabled = enabled;
}

double SparkModel::get_spark_power(double crank_angle) const {
    if (!spark_enabled || !is_spark_active(crank_angle)) {
        return 0.0;
    }
    
    // Convert timing to crank angle (spark occurs before TDC of compression)
    double spark_start = get_start_angle();
    double spark_end = spark_start + spark_duration;  // End angle
    
    // Gaussian energy distribution over spark duration
    double angle_center = (spark_start + spark_end) / 2.0;
    double sigma = spark_duration / 4.0;  // Standard deviation
    
    double power = spark_energy * std::exp(-0.5 * std::pow((crank_angle - angle_center) / sigma, 2)) 
                   / (sigma * std::sqrt(2.0 * M_PI));
    
    return power;
}

bool SparkModel::is_spark_active(double crank_angle) const {
    if (!spark_enabled || spark_energy <= 0.0 || spark_duration <= 0.0) return false;
    
    double spark_start = get_start_angle();
    double spark_end = spark_start + spark_duration;
    
    return (crank_angle >= spark_start && crank_angle <= spark_end);
}
