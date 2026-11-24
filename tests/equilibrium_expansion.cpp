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
    auto thermo = gas->get_solution()->thermo();
    std::vector<double> state;
    thermo->saveState(state);
    thermo->equilibrate("UV");
    require(composition_error(y, gas->mass_fractions()) < 2e-7,
            "Burned gas drifted away from UV equilibrium");
    near(gas->internal_energy(), u, "Equilibration changed internal energy");
    near(gas->density(), rho, "Equilibration changed density");
    thermo->restoreState(state);
}

void check_cold_charge_balances() {
    Cylinder cylinder(0.086, 0.086, 5e-5, "gri30.yaml", 550.0);
    cylinder.initialize_charge(3e-4, 300.0, 101325.0, "N2:1");
    near(cylinder.get_temperature(), 300.0, "Initial temperature changed with piston position");
    near(cylinder.get_pressure(), 101325.0, "Initial pressure changed with piston position");

    const double m0 = cylinder.get_mass();
    const double u0 = m0 * cylinder.get_internal_energy();
    const double s0 = cylinder.get_gas()->entropy();
    cylinder.set_volume(1e-4);
    near(cylinder.get_gas()->entropy(), s0, "Fresh-charge compression changed entropy");
    require(cylinder.get_internal_energy() * m0 > u0, "Compression must do work on the gas");
    cylinder.set_volume(3e-4);
    near(m0 * cylinder.get_internal_energy(), u0, "Reversible compression and expansion changed U");

    const double dt = 1e-4;
    cylinder.heat_loss(dt, 5.0);
    const double heat_rate = cylinder.get_heat_transfer_rate();
    require(heat_rate > 0.0, "Hotter wall must add heat to the gas");
    near(m0 * cylinder.get_internal_energy(), u0 + heat_rate * dt, "Wall heating changed U incorrectly");

    const double before_inflow = m0 * cylinder.get_internal_energy();
    const double dm = 0.15 * m0;
    CanteraInterface reservoir("gri30.yaml");
    reservoir.set_state_TPX(430.0, 180000.0, "AR:1");
    const double expected_enthalpy = dm * reservoir.enthalpy();
    const auto inflow = cylinder.add_mass_flow(dm / dt, dt, "AR:1", 430.0, 180000.0);
    near(inflow.mass / dm, 1.0, "Inflow returned the wrong mass");
    near(inflow.enthalpy, expected_enthalpy, "Inflow returned the wrong enthalpy");
    const double mixed_mass = m0 + dm;
    near(cylinder.get_mass() / mixed_mass, 1.0, "Inflow changed mass incorrectly");
    near(mixed_mass * cylinder.get_internal_energy(), before_inflow + expected_enthalpy,
         "Inflow changed U incorrectly");
    const auto thermo = cylinder.get_gas()->get_solution()->thermo();
    near(thermo->elementalMassFraction(thermo->elementIndex("Ar")), dm / mixed_mass,
         "Inflow failed to conserve species mass");

    const auto no_flow = cylinder.remove_mass_flow(123.0, 0.0);
    require(no_flow.mass == 0.0 && no_flow.enthalpy == 0.0, "Zero-time flow removed mass or energy");

    const double before_discharge = cylinder.get_mass() * cylinder.get_internal_energy();
    const auto outflow = cylinder.remove_mass_flow(2.0 * mixed_mass / dt, dt);
    near(cylinder.get_mass() / mixed_mass, 0.1, "Large discharge violated the mass floor");
    near(outflow.mass / mixed_mass, 0.9, "Discharge reported more mass than it removed");
    near(cylinder.get_mass() * cylinder.get_internal_energy(), before_discharge - outflow.enthalpy,
         "Capped discharge failed the enthalpy balance");
    require(std::isfinite(cylinder.get_temperature()) && cylinder.get_temperature() > 0.0,
            "Capped discharge produced an invalid state");
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
        const double ignition_density = eq.get_density();
        eq.initiate_spark_ignition(1.0);
        near(eq.get_gas()->internal_energy(), ignition_u, "Ignition did not conserve U");
        near(eq.get_density(), ignition_density, "Ignition did not conserve density");
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
        require(eq.get_heat_transfer_rate() < 0.0, "Cooler wall must remove heat from the gas");
        near(m * eq.get_gas()->internal_energy(), before_heat + dt * eq.get_heat_transfer_rate(),
             "Wall heat plus recombination did not conserve energy");
        check_equilibrium(eq);

        const double before_outflow = m * eq.get_gas()->internal_energy();
        const double h = eq.get_gas()->enthalpy(), dm = 0.01 * m;
        const auto exhaust = eq.remove_mass_flow(dm / dt, dt);
        near(eq.get_mass() / m, 0.99, "Incorrect discharged mass");
        near(exhaust.mass / dm, 1.0, "Incorrect reported discharge mass");
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
        check_cold_charge_balances();
        std::cout << "Ignition, expansion, wall heat, mass flow, and cycle reset checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
