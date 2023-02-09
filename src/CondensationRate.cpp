#include <iostream>
#include <cantera/kinetics/CondensationRate.h>
#include <cantera/kinetics/Reaction.h>
#include <cantera/kinetics/Kinetics.h>
#include <cantera/thermo/ThermoPhase.h>

namespace Cantera
{

bool CondensationData::update(const ThermoPhase& phase, const Kinetics& kin)
{
  double T = phase.temperature();
  double P = phase.pressure();

  if (P != pressure || T != temperature)
    update(T, P);

  phase.getMoleFractions(mole_fractions.data());

  return true;
}

void CondensationRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
  ReactionRate::setParameters(node, rate_units);
  if (!node.hasKey("rate-constant")) {
    return;
  }

  auto& rateNode = node["rate-constant"];
  m_A = rateNode["A"].asDouble();
  m_p = rateNode["p"].asDouble();
  m_a = rateNode["a"].asDouble();
  m_b = rateNode["b"].asDouble();
}

void CondensationRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
  std::string vapor_name = rxn.reactants.begin()->first;
  m_vapor_index = kin.kineticsSpeciesIndex(vapor_name);
}

double CondensationRate::evalFromStruct(const CondensationData& shared_data)
{
  double svp = saturation_vapor_pressure(shared_data.temperature);
  double mixr = shared_data.mole_fractions[m_vapor_index];
  if (mixr*shared_data.pressure > svp)
    return m_A;

  return 0.;
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
  rateNode["p"].setQuantity(m_p, "pa");
  rateNode["a"].setQuantity(m_a, "1");
  rateNode["b"].setQuantity(m_b, "K");

  if (!rateNode.empty()) {
    node["rate-constants"] = std::move(rateNode);
  }
}

}
