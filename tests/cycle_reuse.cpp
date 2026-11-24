#include "Engine.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::invalid_argument("Expected the example configuration directory");
        const std::string example_dir = argv[1];
        const auto parameters = ConfigLoader::loadEngineParameters(example_dir + "/engine_parameters.yaml");
        const auto simulation =
            ConfigLoader::loadSimulationParameters(example_dir + "/simulation_parameters.yaml");
        const auto& geometry = parameters.geometry;
        const auto& operating = parameters.operating_conditions;
        Engine engine(geometry.bore, geometry.stroke, geometry.connecting_rod_length,
                      geometry.compression_ratio, parameters.mechanism_file, simulation, parameters);
        engine.set_engine_speed(operating.engine_speed);
        engine.set_boundary_conditions(operating.intake_pressure, operating.intake_temperature,
                                       operating.exhaust_pressure);
        engine.set_intake_composition(operating.intake_composition);
        engine.set_valve_profiles(
            {0.0, 90.0, 180.0, 720.0}, {0.0, parameters.valves.intake.max_lift, 0.0, 0.0},
            {0.0, 540.0, 630.0, 720.0}, {0.0, 0.0, parameters.valves.exhaust.max_lift, 0.0});

        engine.run_cycle();
        const double fired_power = engine.get_brake_power();
        const double fired_efficiency = engine.get_cycle_efficiency();
        if (!(fired_power > 0.0 && fired_efficiency > 0.0)) {
            throw std::runtime_error("Initial fired cycle produced invalid performance");
        }

        engine.enable_spark(false);
        engine.run_cycle();
        if (engine.get_cycle_efficiency() != 0.0 || !(engine.get_brake_power() < 0.0)) {
            throw std::runtime_error("Spark-off cycle carried over fired-cycle state or fuel energy");
        }

        engine.enable_spark(true);
        engine.run_cycle();
        if (!(std::abs(engine.get_brake_power() - fired_power) < 1e-6 &&
              std::abs(engine.get_cycle_efficiency() - fired_efficiency) < 1e-9)) {
            throw std::runtime_error("Repeated fired cycle did not reproduce the original result");
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
