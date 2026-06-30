#include "CanteraInterface.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/base/Solution.h"
#include <stdexcept>
#include <iostream>
#include <cmath>

CanteraInterface::CanteraInterface(const std::string& mechanism_file) {
    try {
        gas = Cantera::newSolution(mechanism_file);
        thermo = gas->thermo().get();
        kinetics = gas->kinetics().get();
        
        // Initialize reactor for chemistry integration with clone=true to avoid deprecation warning
        reactor = std::make_unique<Cantera::IdealGasReactor>(gas, true);
        
        // Initialize reactor net with shared_ptr to reactor to avoid deprecation warning
        // Create a shared_ptr that doesn't own the reactor (non-owning shared_ptr)
        std::shared_ptr<Cantera::ReactorBase> reactor_shared(reactor.get(), [](Cantera::ReactorBase*){});
        reactor_net = std::make_unique<Cantera::ReactorNet>(reactor_shared);
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to initialize Cantera: " + std::string(e.what()));
    }
}

void CanteraInterface::set_state_TPY(double T, double P, const std::vector<double>& Y) {
    thermo->setState_TPY(T, P, Y.data());
}

void CanteraInterface::set_state_TPX(double T, double P, const std::string& X) {
    thermo->setState_TPX(T, P, X);
}

void CanteraInterface::set_state_TP(double T, double P) {
    thermo->setState_TP(T, P);
}

void CanteraInterface::advance_chemistry(double dt) {
    if (!kinetics || dt <= 0.0) {
        return;
    }
    
    try {
        // Set reactor state to current thermo state
        reactor->syncState();
        
        // Advance the reactor network
        double t0 = reactor_net->time();
        //thermo->equilibrate("HP");
	auto& precontents = reactor->contents();
	std::cout << "\n[DEBUG] reactor contents before T=" << precontents.temperature() 
	          << " Y[0]=" << precontents.massFraction(0) 
	          << " Y[2]=" << precontents.massFraction(2) 
	          << " Y[3]=" << precontents.massFraction(3) 
	          << " Y[4]=" << precontents.massFraction(4) 
	          << " Y[5]=" << precontents.massFraction(5) 
		  << std::endl;
	std::cout << "[DEBUG] external thermo before T=" << thermo->temperature() 
	          << " Y[0]=" << thermo->massFraction(0) 
	          << " Y[2]=" << thermo->massFraction(2) 
	          << " Y[3]=" << thermo->massFraction(3) 
	          << " Y[4]=" << thermo->massFraction(4) 
	          << " Y[5]=" << thermo->massFraction(5) 
		  << std::endl;
        reactor_net->advance(t0 + dt);
	auto& contents = reactor->contents();
	std::cout << "[DEBUG] reactor contents T=" << contents.temperature() 
	          << " Y[0]=" << contents.massFraction(0) 
	          << " Y[2]=" << contents.massFraction(2) 
	          << " Y[3]=" << contents.massFraction(3) 
	          << " Y[4]=" << contents.massFraction(4) 
	          << " Y[5]=" << contents.massFraction(5) 
		  << std::endl;
	std::cout << "[DEBUG] external thermo T=" << thermo->temperature() 
	          << " Y[0]=" << thermo->massFraction(0) 
	          << " Y[2]=" << thermo->massFraction(2) 
	          << " Y[3]=" << thermo->massFraction(3) 
	          << " Y[4]=" << thermo->massFraction(4) 
	          << " Y[5]=" << thermo->massFraction(5) 
		  << std::endl;
        
        // The reactor automatically updates the gas state
        
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Cantera reactor integration failed: " << e.what() << std::endl;
        std::cerr << "[ERROR] State at failure: T=" << thermo->temperature() 
                  << " K, P=" << thermo->pressure() / 1e5 << " bar" << std::endl;
        // Fallback: no chemistry advancement
    }
}

double CanteraInterface::temperature() const {
    return thermo->temperature();
}

double CanteraInterface::pressure() const {
    return thermo->pressure();
}

double CanteraInterface::density() const {
    return thermo->density();
}

double CanteraInterface::cp() const {
    return thermo->cp_mass();
}

double CanteraInterface::cv() const {
    return thermo->cv_mass();
}

double CanteraInterface::enthalpy() const {
    return thermo->enthalpy_mass();
}

double CanteraInterface::internal_energy() const {
    return thermo->intEnergy_mass();
}

std::vector<double> CanteraInterface::mass_fractions() const {
    std::vector<double> Y(thermo->nSpecies());
    thermo->getMassFractions(Y.data());
    return Y;
}

double CanteraInterface::mean_molecular_weight() const {
    return thermo->meanMolecularWeight();
}

double CanteraInterface::entropy() const {
    return thermo->entropy_mass();
}

void CanteraInterface::set_state_SV(double s, double v) {
    thermo->setState_SV(s, v);
}

void CanteraInterface::set_state_SD(double s, double rho) {
    // Cantera uses setState_SP (entropy-pressure), not setState_SD
    // We need to calculate pressure from density first
    double T_current = thermo->temperature();
    double R_specific = 8314.0 / thermo->meanMolecularWeight();
    double P_estimated = rho * R_specific * T_current;
    
    // Set state with entropy and estimated pressure, then iterate to correct density
    thermo->setState_SP(s, P_estimated);
    
    // Iterate to get the correct density
    for (int iter = 0; iter < 10; iter++) {
        double current_rho = thermo->density();
        if (std::abs(current_rho - rho) / rho < 1e-6) break;
        
        // Adjust pressure
        P_estimated *= rho / current_rho;
        thermo->setState_SP(s, P_estimated);
    }
}

void CanteraInterface::set_state_HP(double u, double v) {
    thermo->setState_HP(u, v);
}
void CanteraInterface::set_state_UV(double u, double v) {
    thermo->setState_UV(u, v);
}

double CanteraInterface::gamma() const {
    double cp_val = thermo->cp_mass();
    double cv_val = thermo->cv_mass();
    
    if (cv_val > 0.0) {
        return cp_val / cv_val;
    } else {
        return 1.4; // Default value for air if cv is invalid
    }
}

std::vector<double> CanteraInterface::net_production_rates() const {
    std::vector<double> wdot(thermo->nSpecies());
    if (kinetics) {
        kinetics->getNetProductionRates(wdot.data());
    }
    return wdot;
}

double CanteraInterface::heat_release_rate() const {
    if (!kinetics) return 0.0;
    
    std::vector<double> wdot(thermo->nSpecies());
    std::vector<double> h_species(thermo->nSpecies());
    
    kinetics->getNetProductionRates(wdot.data());
    thermo->getPartialMolarEnthalpies(h_species.data());
    
    double heat_release = 0.0;
    for (size_t k = 0; k < thermo->nSpecies(); k++) {
        heat_release -= wdot[k] * h_species[k];
    }
    
    return heat_release;
}

std::vector<std::string> CanteraInterface::species_names() const {
    std::vector<std::string> names;
    for (size_t k = 0; k < thermo->nSpecies(); k++) {
        names.push_back(thermo->speciesName(k));
    }
    return names;
}

size_t CanteraInterface::n_species() const {
    return thermo->nSpecies();
}
