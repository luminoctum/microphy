#include <iostream>
#include <cantera/kinetics/NucleationKinetics.h>
#include <cantera/thermo/ThermoPhase.h>

namespace Cantera
{

void NucleationKinetics::updateROP()
{
  std::cout << "######### I'm updating ROP ##########" << std::endl;

  // set temperature of the surface surface to be the vapor phase temperature
  m_thermo[surfacePhaseIndex()]->setTemperature(m_vapor_phase->temperature());

  // regular concentration and activity concentration
  for (size_t n = 0; n < nPhases(); ++n) {
    ThermoPhase const* tp = m_thermo[n];
    tp->getConcentrations(m_conc.data() + m_start[n]);
    tp->getActivityConcentrations(m_actConc.data() + m_start[n]);
  }

  std::cout << "concentrations = " << std::endl;
  for (int i = 0; i < m_kk; ++i)
    std::cout << m_conc[i] << " ";
  std::cout << std::endl;


  // loop over interface rates evaluator
  for (auto& rates : m_interfaceRates) {
    bool changed = rates->update(thermo(surfacePhaseIndex()), *this);
    if (changed) {
      rates->getRateConstants(m_rfn.data());
    }
  }

  // revise activity concentration
  for (auto& rates : m_interfaceRates) {
    rates->reviseActivityConcentration(m_actConc.data(), m_conc.data());
  }

  std::cout << "activity concentrations = " << std::endl;
  for (int i = 0; i < m_kk; ++i)
    std::cout << m_actConc[i] << " ";
  std::cout << std::endl;

  for (size_t i = 0; i < nReactions(); i++) {
    // Scale the base forward rate coefficient by the perturbation factor
    m_ropf[i] = m_rfn[i] * m_perturb[i];
    // Multiply the scaled forward rate coefficient by the reciprocal of the
    // equilibrium constant
    m_ropr[i] = m_ropf[i] * m_rkcn[i];
  }

  // multiply ropf by the activity concentration reaction orders to obtain
  // the forward rates of progress.
  m_reactantStoich.multiply(m_actConc.data(), m_ropf.data());

  // For reversible reactions, multiply ropr by the activity concentration
  // products
  m_revProductStoich.multiply(m_actConc.data(), m_ropr.data());

  for (size_t j = 0; j != nReactions(); ++j) {
    m_ropnet[j] = m_ropf[j] - m_ropr[j];
  }

  m_ROP_ok = true;
}

}
