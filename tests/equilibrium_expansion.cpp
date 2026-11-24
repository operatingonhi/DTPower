#include "Cylinder.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, const char* message) {
    require(std::abs(actual - expected) < 2e-7 * std::max(1.0, std::abs(expected)), message);
}

void initialize(Cylinder& cylinder) {
    auto* gas = cylinder.get_gas();
    gas->set_state_TPX(650.0, 4e6, "CH4:0.3,O2:0.6,AR:0.1");
    cylinder.set_mass(gas->density() * cylinder.get_volume());
}

double composition_error(const std::vector<double>& a, const std::vector<double>& b) {
    double error = 0.0;
    for (size_t k = 0; k < a.size(); ++k) error = std::max(error, std::abs(a[k] - b[k]));
    return error;
}

void check_equilibrium(Cylinder& cylinder) {
    auto* gas = cylinder.get_gas();
    auto y = gas->mass_fractions();
    const double u = gas->internal_energy(), rho = gas->density();
    gas->get_solution()->thermo()->equilibrate("UV");
    require(composition_error(y, gas->mass_fractions()) < 2e-7,
            "Burned gas drifted away from UV equilibrium");
    near(gas->internal_energy(), u, "Equilibration changed internal energy");
    near(gas->density(), rho, "Equilibration changed density");
}
}

int main() {
    try {
        Cylinder eq(0.086, 0.086, 5e-5, "gri30.yaml", 470.0);
        initialize(eq);
        require(eq.get_combustion_model() == CombustionModel::EQUILIBRIUM,
                "Default model must remain equilibrium");
        require(eq.get_equilibrium_expansion(), "Continuing equilibrium must default on");
        const double ignition_u = eq.get_gas()->internal_energy();
        eq.initiate_spark_ignition(1.0);
        near(eq.get_gas()->internal_energy(), ignition_u, "Ignition did not conserve U");
        const auto ignition_y = eq.get_gas()->mass_fractions();
        const auto thermo = eq.get_gas()->get_solution()->thermo();
        const auto c_index = thermo->elementIndex("C"), h_index = thermo->elementIndex("H");
        const double carbon = thermo->elementalMassFraction(c_index);
        const double hydrogen = thermo->elementalMassFraction(h_index);
        const double entropy = eq.get_gas()->entropy();
        eq.set_volume(4e-4);
        near(eq.get_gas()->entropy(), entropy, "Adiabatic expansion did not conserve entropy");
        near(thermo->elementalMassFraction(c_index), carbon, "Expansion changed carbon inventory");
        near(thermo->elementalMassFraction(h_index), hydrogen, "Expansion changed hydrogen inventory");
        require(composition_error(ignition_y, eq.get_gas()->mass_fractions()) > 0.01,
                "Expansion did not allow appreciable recombination");
        check_equilibrium(eq);

        const double m = eq.get_mass(), before_heat = m * eq.get_gas()->internal_energy();
        const double dt = 1e-5;
        eq.heat_loss(dt, 5.0);
        near(m * eq.get_gas()->internal_energy(), before_heat + dt * eq.get_heat_transfer_rate(),
             "Wall heat plus recombination did not conserve energy");
        check_equilibrium(eq);

        const double before_outflow = m * eq.get_gas()->internal_energy();
        const double h = eq.get_gas()->enthalpy(), dm = 0.01 * m;
        eq.remove_mass_flow(dm / dt, dt);
        near(eq.get_mass(), m - dm, "Incorrect discharged mass");
        near(eq.get_mass() * eq.get_gas()->internal_energy(), before_outflow - dm * h,
             "Equilibrium discharge did not conserve the enthalpy-flow energy balance");
        check_equilibrium(eq);
        eq.set_volume(4.1e-4);
        check_equilibrium(eq);  // Outflow must not disable subsequent equilibration.

        Cylinder frozen(0.086, 0.086, 5e-5, "gri30.yaml", 470.0);
        initialize(frozen);
        frozen.set_equilibrium_expansion(false);
        frozen.initiate_spark_ignition(1.0);
        const auto frozen_y = frozen.get_gas()->mass_fractions();
        frozen.set_volume(4e-4);
        frozen.heat_loss(dt, 5.0);
        frozen.remove_mass_flow(0.01 * frozen.get_mass() / dt, dt);
        require(composition_error(frozen_y, frozen.get_gas()->mass_fractions()) < 1e-12,
                "Frozen comparison mode changed composition");

        // A cycle reset must prevent equilibrium chemistry in unignited charge.
        eq.reset_combustion_state();
        initialize(eq);
        const auto fresh_y = eq.get_gas()->mass_fractions();
        eq.set_volume(4.5e-4);
        eq.heat_loss(dt, 5.0);
        require(composition_error(fresh_y, eq.get_gas()->mass_fractions()) < 1e-12,
                "Fresh charge reacted before ignition after cycle reset");
        std::cout << "Equilibrium expansion, heat, discharge, reset, and frozen checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
