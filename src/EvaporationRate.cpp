
#include <cantera/kinetics/EvaporationRate.h>

namespace Cantera
{

bool EvaporationData::update(const ThermoPhase& phase, const Kinetics& kin)
{
  return true;
}

void EvaporationRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
}

double EvaporationRate::evalFromStruct(const EvaporationData& shared_data)
{
  return 0.;
}

void EvaporationRate::getParameters(AnyMap& node) const
{}

}
