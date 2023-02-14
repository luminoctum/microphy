#include <iostream>
#include <cantera/kinetics/CondensationRate.h>
#include <cantera/kinetics/Reaction.h>
#include <cantera/kinetics/Kinetics.h>
#include <cantera/thermo/ThermoPhase.h>

namespace Cantera
{

// different phase may have different temperature
// phase 0 is the current interface phase (dummy)
// phase 1 is the vapor phase
// phase 2,3,4 ... are the condensed phases
bool CondensationData::update(const ThermoPhase& phase, const Kinetics& kin)
{
  double T = phase.temperature();
  std::cout << "==== CondensationData::update ====" << std::endl;

  if (T != temperature) {
    update(T);
    std::cout << "temperature = " << T << std::endl;
    std::cout << "==================================" << std::endl;
    return true;
  }


  return false;
}

// [conc] = kmol/m^3
void CondensationRate::getActivityConcentration(double *actConc,
  double const *conc, CondensationData const& shared_data)
{
  double svp = saturation_vapor_pressure(shared_data.temperature);
  if ((shared_data.temperature < m_T3 && m_ice) ||
      (shared_data.temperature > m_T3 && !m_ice))
  {
    actConc[m_vapor_index] = conc[m_vapor_index] - svp/(GasConstant*shared_data.temperature);
    if (actConc[m_vapor_index] > 0.) {
      actConc[m_solid_index] = 0;
    } else {
      actConc[m_vapor_index] = 0;
    }
  }
}

void CondensationRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
  ReactionRate::setParameters(node, rate_units);
  if (!node.hasKey("rate-constant")) {
    return;
  }

  auto& rateNode = node["rate-constant"];
  m_A = rateNode["A"].asDouble();
  m_P3 = rateNode["P3"].asDouble();
  m_T3 = rateNode["T3"].asDouble();
  m_beta = rateNode["beta"].asDouble();
  m_gamma = rateNode["gamma"].asDouble();
  m_ice = rateNode["ice"].asBool();

  std::cout << "==== CondensationRate::setParameter ====" << std::endl;
  std::cout << "m_A = " << m_A << std::endl;
  std::cout << "m_P3 = " << m_P3 << std::endl;
  std::cout << "m_T3 = " << m_T3 << std::endl;
  std::cout << "m_beta = " << m_beta << std::endl;
  std::cout << "m_gamma = " << m_gamma << std::endl;
  std::cout << "m_ice = " << m_ice << std::endl;
  std::cout << "========================================" << std::endl;
}

void CondensationRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
  std::string vapor_name = rxn.reactants.begin()->first;
  std::string solid_name = rxn.products.begin()->first;
  m_vapor_index = kin.kineticsSpeciesIndex(vapor_name);
  m_solid_index = kin.kineticsSpeciesIndex(solid_name);
  std::cout << "==== CondensationRate::setContext ====" << std::endl;
  std::cout << "m_vapor_index = " << m_vapor_index << std::endl;
  std::cout << "m_solid_index = " << m_solid_index << std::endl;
  std::cout << "=======================================" << std::endl;
}

double CondensationRate::evalFromStruct(const CondensationData& shared_data)
{
  if ((shared_data.temperature <= m_T3 && m_ice) ||
      (shared_data.temperature >= m_T3 && !m_ice)) {
    return m_A;
  } else {
    return 0.;
  }
}

void CondensationRate::getParameters(AnyMap& node) const
{
  node["type"] = type();
  if (!valid()) {
    // object not fully set up
    return;
  }
  AnyMap rateNode;
  rateNode["A"].setQuantity(m_A, "s^{-1}");
  rateNode["P3"].setQuantity(m_P3, "pa");
  rateNode["T3"].setQuantity(m_T3, "K");
  rateNode["beta"].setQuantity(m_beta, "1");
  rateNode["gamma"].setQuantity(m_gamma, "1");
  rateNode["ice"].setQuantity(m_ice, "1");

  if (!rateNode.empty()) {
    node["rate-constants"] = std::move(rateNode);
  }
}

}
