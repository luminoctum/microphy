#include <iostream>
#include <cantera/kinetics/FreezingRate.h>
#include <cantera/thermo/ThermoPhase.h>

namespace Cantera
{

bool FreezingData::update(const ThermoPhase& phase, const Kinetics& kin)
{
  double T = phase.temperature();
  if (T == temperature) {
    return false;
  }

  update(T);
  return true;
}

void FreezingRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
  ReactionRate::setParameters(node, rate_units);
  if (!node.hasKey("rate-constant")) {
    return;
  }

  auto& rateNode = node["rate-constant"];
  freezing_temperature = rateNode["T"].asDouble();
  freezing_rate = rateNode["A"].asDouble();
}

void FreezingRate::getParameters(AnyMap& node) const
{
  node["type"] = type();
  if (!valid()) {
    // object not fully set up
    return;
  }
  AnyMap rateNode;
  rateNode["T"].setQuantity(freezing_temperature, "K");
  rateNode["A"].setQuantity(freezing_rate, "kg^{-1} s^{-1}");

  if (!rateNode.empty()) {
    node["rate-constants"] = std::move(rateNode);
  }
}

void FreezingRate::updateFromStruct(const FreezingData& shared_data)
{
}

double FreezingRate::evalFromStruct(const FreezingData& shared_data)
{
  if (shared_data.temperature > freezing_temperature)
    return -freezing_rate;

  if (shared_data.temperature < freezing_temperature)
    return freezing_rate;

  return 0.;
}

double FreezingRate::ddTScaledFromStruct(const FreezingData& shared_data) const
{
  return 0.;
}

//! Check basic syntax and settings of reaction rate expression
void FreezingRate::check(const std::string& equation)
{}

void FreezingRate::validate(const std::string& equation, const Kinetics& kin)
{}

}
