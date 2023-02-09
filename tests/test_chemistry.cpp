
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
  auto vapor_cloud = newSolution("H2O_kessler.yaml", "vapor:cloud");

  Real ymass[2] = {0.9, 0.1};
  Real xmol[4];
  Real temperature = 300.;
  Real density = 1.;

  // check out vapor species
  vapor->thermo()->setMassFractions_NoNorm(ymass);
  vapor->thermo()->setTemperature(temperature);
  vapor->thermo()->setDensity(density);
  vapor->thermo()->getMoleFractions(xmol);
  std::cout << "Name of the phase = "
            << vapor->thermo()->name() << std::endl;
  std::cout << "Number of elements = "
            << vapor->thermo()->nElements() << std::endl;
  std::cout << "Element names = ";
  for (int i = 0; i < vapor->thermo()->nElements(); ++i)
    std::cout << vapor->thermo()->elementName(i) << " ";
  std::cout << std::endl;
  std::cout << "Number of vapor species = " 
            << vapor->thermo()->nSpecies() << std::endl;
  std::cout << "Species names = ";
  for (int i = 0; i < vapor->thermo()->nSpecies(); ++i)
    std::cout << vapor->thermo()->speciesName(i) << " ";
  std::cout << std::endl;

  std::cout << "Temperature = " << vapor->thermo()->temperature() << std::endl;
  std::cout << "Pressure = " << vapor->thermo()->pressure() << std::endl;
  std::cout << "density = " << vapor->thermo()->density() << std::endl;

  vapor->thermo()->getMoleFractions(xmol);
  std::cout << "Mole fractions = ";
  for (int i = 0; i < vapor->thermo()->nSpecies(); ++i)
    std::cout << xmol[i] << " ";
  std::cout << std::endl;

  std::cout << "Number of vapor reactions = "
            << vapor->kinetics()->nReactions() << std::endl;
  std::cout << std::endl;

  // check out cloud species
  cloud->thermo()->setMassFractions_NoNorm(ymass);
  cloud->thermo()->setTemperature(temperature*2);
  cloud->thermo()->getMoleFractions(xmol);
  std::cout << "Name of the phase = "
            << cloud->thermo()->name() << std::endl;
  std::cout << "Number of elements = "
            << cloud->thermo()->nElements() << std::endl;
  std::cout << "Element names = ";
  for (int i = 0; i < cloud->thermo()->nElements(); ++i)
    std::cout << cloud->thermo()->elementName(i) << " ";
  std::cout << std::endl;
  std::cout << "Number of cloud species = " 
            << cloud->thermo()->nSpecies() << std::endl;
  std::cout << "Species names = ";
  for (int i = 0; i < cloud->thermo()->nSpecies(); ++i)
    std::cout << cloud->thermo()->speciesName(i) << " ";
  std::cout << std::endl;

  std::cout << "Temperature = " << cloud->thermo()->temperature() << std::endl;
  std::cout << "Pressure = " << cloud->thermo()->pressure() << std::endl;
  std::cout << "density = " << cloud->thermo()->density() << std::endl;

  cloud->thermo()->getMoleFractions(xmol);
  std::cout << "Mole fractions = ";
  for (int i = 0; i < cloud->thermo()->nSpecies(); ++i)
    std::cout << xmol[i] << " ";
  std::cout << std::endl;
  
  std::cout << "Number of cloud reactions = " 
            << cloud->kinetics()->nReactions() << std::endl;
  std::cout << std::endl;

  // check out vapor:cloud phase
  std::cout << "Number of vapor-cloud reactions = "
            << vapor_cloud->kinetics()->nReactions() << std::endl;
  std::cout << "Total species = "
            << vapor_cloud->kinetics()->nTotalSpecies() << std::endl;
  std::cout << "Number of Phases = " << vapor_cloud->kinetics()->nPhases() << std::endl;
  std::cout << vapor_cloud->kinetics()->phaseIndex("vapor:cloud") << std::endl;
  std::cout << vapor_cloud->kinetics()->phaseIndex("vapor") << std::endl;
  std::cout << vapor_cloud->kinetics()->phaseIndex("cloud") << std::endl;
  std::cout << "Species names = ";
  for (int i = 0; i < vapor_cloud->kinetics()->nTotalSpecies(); ++i)
    std::cout << vapor_cloud->kinetics()->kineticsSpeciesName(i) << " ";
  std::cout << std::endl;
}
