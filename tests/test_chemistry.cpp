
// Solution class describes a phase consists of a mixture of chemical species
#include <cantera/base/Solution.h>
// ThermoPhase object stores the thermodynamic state
#include <cantera/thermo.h>

// Kinetics object stores the chemical kinetics information
#include <cantera/kinetics.h>

#include <iostream>

using namespace Cantera;

int main() {
  auto sol = newSolution("H2O_kessler.yaml", "Kessler");
  auto gas = sol->thermo();

  // check out number of species
  std::cout << "Number of species = "
            << gas->nSpecies() << std::endl;
}
