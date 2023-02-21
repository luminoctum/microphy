#ifndef CONSTRUCT_ATMOSPHERE_HPP
#define CONSTRUCT_ATMOSPHERE_HPP

// C/C++ header
#include <vector>
#include <cmath>
#include <memory>

#include <debugger.hpp>
#include <configure.hpp>
#include "suspension.hpp"
#include "rk4_integrate_z.hpp"
#include "rk4_integrate_lnp.hpp"

template<typename VaporPhase, typename CloudPhase>
void ConstructAtmosphere(Real **w, size_t len,
  AtmConstructionMethod method,
  Vapor<VaporPhase> *vapor,
  Suspension<VaporPhase, CloudPhase> *cloud,
  Cantera::Kinetics *kin)
{
  pdebug->Call("Constructing atmosphere");

  // set small amount of trial condensates
  Real ym[2] = {1., 1.};
  cloud->setRelativeMole(1.E-6);
  cloud->setMassFractions(ym);

  std::vector<Real> conc(vapor->nSpecies());
  std::vector<Real> act_conc(kin->nTotalSpecies());
  vapor->getConcentrations(conc.data());
  kin->getActivityConcentrations(act_conc.data());

  // remove supersaturated vapors
  for (int i = 1; i <= NumVapors; ++i) {
    std::string name = vapor->vaporName(i);
    //std::cout << conc[vapor->speciesIndex(name)] << " "
    //          << act_conc[kin->kineticsSpeciesIndex(name)] << std::endl;
    conc[vapor->speciesIndex(name)] -= act_conc[kin->kineticsSpeciesIndex(name)];
  }
  vapor->setConcentrations(conc.data());

  std::cout << vapor->temperature() << std::endl;
  std::cout << vapor->pressure() << std::endl;

  // set density, vapor mass fractions and pressure
  vapor->getMassFractions(conc.data());
  w[0][idn] = vapor->density();
  for (int i = 1; i <= NumVapors; ++i)
    w[0][i] = conc[vapor->nDryMixtures() + i - 1];
  w[0][ipr] = vapor->pressure();

  for (int n = 0; n < NumHydros; ++n) {
    std::cout << w[0][n] << " ";
  }
  std::cout << std::endl;

  for (int i = 1; i < len; ++i) {
    // RK4 integration 
    int niter = 0;

    if (method.dz != 0) {  // height coordinate
      rk4_integrate_z(method, vapor, cloud, kin);
    } else if (method.dlnp != 0) { // pressure coordinate
      rk4_integrate_lnp(method, vapor, cloud, kin);
    } else {
      Debugger::Fatal("construct_atmosphere", "dz and dlnp are both zero");
    }

    // set vapor mass fractions, density and pressure
    vapor->getMassFractions(conc.data());
    w[i][idn] = vapor->density();
    for (int j = 1; j <= NumVapors; ++j)
      w[i][j] = conc[vapor->nDryMixtures() + j - 1];
    w[i][ipr] = vapor->pressure();
  }
}

#endif
