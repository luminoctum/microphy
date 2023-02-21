#include <iostream>
#include <cantera/kinetics/EvaporationRate.h>
#include <cantera/kinetics/Reaction.h>
#include <cantera/kinetics/Kinetics.h>
#include <cantera/thermo/ThermoPhase.h>

namespace Cantera
{

bool EvaporationData::update(const ThermoPhase& phase, const Kinetics& kin)
{
  double T = phase.temperature();
  if (T != temperature) {
    update(T);
    return true;
  }
  return false;
}

// [conc] = kmol/m^3
void EvaporationRate::getActivityConcentration(double *actConc,
  double const *conc, EvaporationData const& shared_data)
{
  /*double svp = saturation_vapor_pressure(shared_data.temperature);
  if ((shared_data.temperature < m_T3 && m_ice) ||
      (shared_data.temperature > m_T3 && !m_ice))
  {
    double sup = conc[m_vapor_index]/(svp/(GasConstant*shared_data.temperature));
    if (sup < 1.) {
      actConc[m_solid_index] = conc[m_solid_index]/sup;
    } else {
      actConc[m_solid_index] = 0.;
    }
  }*/
}

void EvaporationRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
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
}

void EvaporationRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
  std::string vapor_name = rxn.reactants.begin()->first;
  std::string solid_name = rxn.products.begin()->first;
  m_vapor_index = kin.kineticsSpeciesIndex(vapor_name);
  m_solid_index = kin.kineticsSpeciesIndex(solid_name);
}

double EvaporationRate::evalFromStruct(const EvaporationData& shared_data)
{
  if ((shared_data.temperature <= m_T3 && m_ice) ||
      (shared_data.temperature >= m_T3 && !m_ice)) {
    return m_A;
  } else {
    return 0.;
  }
}

void EvaporationRate::getParameters(AnyMap& node) const
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
