#ifndef CONDENSATES_HPP
#define CONDENSATES_HPP

#include <athena_arrays.hpp>
#include "solution_group.hpp"

class MeshBlock;
class ParameterInput;

template<int NumTrackeringFields>
class Condensates : public SolutionGroup {
public:
// public access variables

  // conserved variables (extensive), kg/m^3
  AthenaArray<Real> c;

  // primitive variables (intensive), kg/kg
  AthenaArray<Real> q;

  // temperature and energy
  AthenaArray<Real> tem, energy;

  // differential velocity 1, 2, 3
  AthenaArray<Real> delv1, delv2, delv3;

  // flux in 3-dimensions
  AthenaArray<Real> flux[3];

// functions
  Condensates(MeshBlock *pmb, ParameterInput *pin);

  void calcSedimentaionFlux();

protected:
  PassiveScalars *pscalars_;
};

#endif
