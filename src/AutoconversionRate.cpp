
#include <cantera/kinetics/AutoconversionRate.h>

namespace Cantera {

bool AutoconversionData::update(const ThermoPhase& phase, const Kinetics& kin)
{
  return true;
}

void AutoconversionRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
}

double AutoconversionRate::evalFromStruct(const AutoconversionData& shared_data)
{
  return 0.;
}

void AutoconversionRate::getParameters(AnyMap& node) const
{}

}
