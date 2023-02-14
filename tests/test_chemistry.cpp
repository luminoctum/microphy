
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
  //auto vapor = newSolution("H2O_kessler.yaml", "vapor");
  //auto cloud = newSolution("H2O_kessler.yaml", "cloud");
  auto sol = newSolution("H2O_kessler.yaml", "vapor:cloud");
  auto vapor = sol->adjacent("vapor");
  auto cloud = sol->adjacent("cloud");

  Real ymass[2] = {0.5, 0.5};
  Real xmol[20];
  Real temperature = 300.;
  Real pressure = 1.E5;
  Real density = 1.;

  // check out vapor species
  vapor->thermo()->setMassFractions_NoNorm(ymass);
  vapor->thermo()->setTemperature(temperature);
  vapor->thermo()->setPressure(pressure);
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

  std::cout << "Mole fractions = ";
  for (int i = 0; i < cloud->thermo()->nSpecies(); ++i)
    std::cout << xmol[i] << " ";
  std::cout << std::endl;
  
  std::cout << "Number of cloud reactions = " 
            << cloud->kinetics()->nReactions() << std::endl;
  std::cout << std::endl;

  // check out vapor:cloud phase
  //sol->thermo()->setMassFractions_NoNorm(ymass);
  //sol->thermo()->setTemperature(temperature);
  //sol->thermo()->getMoleFractions(xmol);

  std::cout << "Name of the solution= "
            << sol->name() << std::endl;

  std::cout << "Number of adjacent phases = "
            << sol->nAdjacent() << std::endl;
  std::cout << "Adjacent solution names = ";
  for (int i = 0; i < sol->nAdjacent(); ++i)
    std::cout << sol->adjacent(i)->name() << " ";
  std::cout << std::endl;


  std::cout << "Number of kinetics reactions = "
            << sol->kinetics()->nReactions() << std::endl;
  for (int i = 0; i < sol->kinetics()->nReactions(); ++i) {
    std::cout << sol->kinetics()->reactionString(i) << std::endl;
  }

  std::vector<double> kfwd(sol->kinetics()->nReactions());

  sol->kinetics()->getFwdRateConstants(kfwd.data());
  std::cout << "Forward rate constants = " << std::endl;
  for (int i = 0; i < sol->kinetics()->nReactions(); ++i) {
    std::cout << kfwd[i] << std::endl;
  }

  sol->kinetics()->getFwdRatesOfProgress(kfwd.data());
  std::cout << "Forward rate of progress = " << std::endl;
  for (int i = 0; i < sol->kinetics()->nReactions(); ++i) {
    std::cout << kfwd[i] << std::endl;
  }

  sol->kinetics()->getRevRatesOfProgress(kfwd.data());
  std::cout << "Reverse rate of progress = " << std::endl;
  for (int i = 0; i < sol->kinetics()->nReactions(); ++i) {
    std::cout << kfwd[i] << std::endl;
  }

  sol->kinetics()->getNetRatesOfProgress(kfwd.data());
  std::cout << "Net rate of progress = " << std::endl;
  for (int i = 0; i < sol->kinetics()->nReactions(); ++i) {
    std::cout << kfwd[i] << std::endl;
  }

  std::cout << "Total species = "
            << sol->kinetics()->nTotalSpecies() << std::endl;
  std::cout << "Number of Phases = " << sol->kinetics()->nPhases() << std::endl;
  std::cout << sol->kinetics()->phaseIndex("vapor:cloud") << std::endl;
  std::cout << sol->kinetics()->phaseIndex("vapor") << std::endl;
  std::cout << sol->kinetics()->phaseIndex("cloud") << std::endl;
  std::cout << "Species names = ";
  for (int i = 0; i < sol->kinetics()->nTotalSpecies(); ++i)
    std::cout << sol->kinetics()->kineticsSpeciesName(i) << " ";
  std::cout << std::endl;
}
