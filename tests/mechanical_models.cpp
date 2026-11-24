#include "Crankshaft.h"
#include "FrictionModel.h"
#include "SparkModel.h"
#include "ValveFlow.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const char* message) {
    require(std::abs(actual - expected) <= tolerance, message);
}

void check_crankshaft() {
    constexpr double stroke = 0.086, bore = 0.086, clearance = 5e-5, rpm = 2000.0;
    Crankshaft crankshaft(stroke, 0.143, bore);
    const double piston_area = std::acos(-1.0) * bore * bore / 4.0;
    near(crankshaft.calculate_cylinder_volume(0.0, clearance), clearance, 1e-14,
         "Wrong cylinder volume at top dead center");
    near(crankshaft.calculate_cylinder_volume(180.0, clearance), clearance + piston_area * stroke,
         1e-14, "Wrong cylinder volume at bottom dead center");
    near(crankshaft.calculate_piston_position(67.0), crankshaft.calculate_piston_position(427.0),
         1e-14, "Piston position is not periodic");
    const double step = 1e-3;
    const double numerical_velocity = (crankshaft.calculate_piston_position(67.0 + step) -
                                       crankshaft.calculate_piston_position(67.0 - step)) *
                                      (6.0 * rpm) / (2.0 * step);
    near(crankshaft.calculate_piston_velocity(67.0, rpm), numerical_velocity, 1e-7,
         "Piston velocity disagrees with the derivative of position");
}

void check_valves() {
    ValveFlow valve(0.034, 0.030, 0.6);
    valve.set_lift_profile({0.0, 90.0, 180.0, 720.0}, {0.0, 0.008, 0.0, 0.0});
    near(valve.get_lift(45.0), 0.004, 1e-14, "Wrong interpolated valve lift");
    near(valve.get_lift(-675.0), 0.004, 1e-14, "Negative angle did not wrap correctly");
    near(valve.get_lift(765.0), 0.004, 1e-14, "Valve timing did not repeat");
    require(valve.get_lift(540.0) == 0.0, "Intake valve opened in the exhaust stroke");

    const auto flow = [&valve](double lift, double inlet, double outlet) {
        return valve.calculate_mass_flow(lift, inlet, outlet, 300.0, 1.0);
    };
    require(flow(0.0, 200000.0, 100000.0) == 0.0, "Closed valve passed mass");
    require(flow(0.008, 100000.0, 100000.0) == 0.0, "Equal pressures produced mass flow");
    const double choked = flow(0.008, 200000.0, 100000.0);
    require(choked > 0.0 && std::isfinite(choked), "Choked flow is invalid");
    near(flow(0.008, 200000.0, 50000.0), choked, 1e-12, "Choked flow depends on outlet pressure");
    near(flow(0.008, 100000.0, 200000.0), -choked, 1e-12, "Backflow has the wrong sign");
    const double subsonic = flow(0.008, 200000.0, 160000.0);
    require(subsonic > 0.0 && subsonic < choked, "Subsonic flow outside physical range");

    for (const auto& angles : {std::vector<double>{0.0}, std::vector<double>{1.0, 1.0}}) {
        bool rejected = false;
        try {
            valve.set_lift_profile(angles, {0.0, 0.008});
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "Invalid valve lift profile was accepted");
    }
}

void check_spark_and_friction() {
    SparkModel spark(10.0, 2.0, 1.0);
    require(spark.get_spark_power(349.9) == 0.0 && spark.get_spark_power(352.1) == 0.0,
            "Spark fired outside the selected window");
    require(spark.get_spark_power(351.0) > spark.get_spark_power(350.0),
            "Spark power did not peak at the center of the window");
    spark.enable_spark(false);
    require(spark.get_spark_power(351.0) == 0.0, "Disabled spark delivered power");

    FrictionModel friction(0.086, 0.086, 0.143);
    require(friction.calculate_instantaneous_friction_power(2000.0, 1e6, 0.0, 0.0) == 0.0,
            "Stationary piston dissipated friction power");
    require(friction.calculate_instantaneous_friction_power(2000.0, 1e6, 270.0, -5.0) > 0.0,
            "Friction failed to dissipate energy on the return stroke");
    friction.start_friction_integration(0.0);
    friction.accumulate_friction_work(20.0, 0.5);
    near(friction.get_integrated_fmep(5e-4), 20000.0, 1e-8, "Incorrect friction integration");
    friction.reset_friction_integration();
    require(friction.get_integrated_fmep(5e-4) == 0.0, "Friction integration did not reset");
}
}

int main() {
    try {
        check_crankshaft();
        check_valves();
        check_spark_and_friction();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
