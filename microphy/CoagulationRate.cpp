
#include <cantera/kinetics/CoagulationRate.h>

namespace Cantera {

bool CoagulationData::update(const ThermoPhase& phase, const Kinetics& kin)
{
  return true;
}

void CoagulationRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
}

double CoagulationRate::evalFromStruct(const CoagulationData& shared_data)
{
  return 0.;
}

void CoagulationRate::getParameters(AnyMap& node) const
{}

}
