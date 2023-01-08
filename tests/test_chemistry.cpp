
#include <iostream>

// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>
// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>

#include <configure.hpp>

using namespace Cantera;

int main() {
  auto vapor = newSolution("H2O_kessler.yaml", "vapor");
  auto cloud = newSolution("H2O_kessler.yaml", "cloud");
  auto precip = newSolution("H2O_kessler.yaml", "precipitation");

  Real ymass[2] = {0.9, 0.1};
  Real xmol[2];
  Real temperature = 300.;
  Real density = 1.;

  // check out number of vapor species
  vapor->thermo()->setMassFractions_NoNorm(ymass);
  vapor->thermo()->setTemperature(temperature);
  vapor->thermo()->setDensity(density);
  vapor->thermo()->getMoleFractions(xmol);
  std::cout << "Number of vapor species = " 
            << vapor->thermo()->nSpecies() << std::endl;
  std::cout << "Number of vapor reactions = "
            << vapor->kinetics()->nReactions() << std::endl;
  std::cout << "Temperature = " << vapor->thermo()->temperature() << std::endl;
  std::cout << "Pressure = " << vapor->thermo()->pressure() << std::endl;
  std::cout << "density = " << vapor->thermo()->density() << std::endl;
  vapor->thermo()->getMoleFractions(xmol);
  std::cout << "Mole fractions = ";
  for (int i = 0; i < vapor->thermo()->nSpecies(); ++i)
    std::cout << xmol[i] << " ";
  std::cout << std::endl;
  std::cout << vapor->thermo()->phaseOfMatter() << std::endl;

  // check out number of cloud species
  std::cout << "Number of cloud species = " 
            << cloud->thermo()->nSpecies() << std::endl;
  std::cout << "Number of cloud reactions = " 
            << cloud->kinetics()->nReactions() << std::endl;
  std::cout << vapor->thermo()->phaseOfMatter() << std::endl;

  // check out number of precip species
  std::cout << "Number of precipitation species = "
            << precip->thermo()->nSpecies() << std::endl;
  std::cout << "Number of precipitation species = "
            << precip->kinetics()->nReactions() << std::endl;
  std::cout << vapor->thermo()->phaseOfMatter() << std::endl;
}
