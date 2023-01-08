
#include <cantera/kinetics/FreezingRate.h>
#include <cantera/thermo/ThermoPhase.h>

namespace Cantera
{

void FreezingData::update(double T)
{}

void FreezingData::update(double T, double P)
{}

bool FreezingData::update(const ThermoPhase& phase, const Kinetics& kin)
{}

FreezingRate::FreezingRate()
{}

void FreezingRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{}

void FreezingRate::getParameters(AnyMap& rateNode, const Units& rate_units) const
{}

void FreezingRate::updateFromStruct(const FreezingData& shared_data)
{}

double FreezingRate::evalFromStruct(const FreezingData& shared_data)
{
  return 0.;
}

void FreezingRate::validate(const std::string& equation, const Kinetics& kin)
{}

}
