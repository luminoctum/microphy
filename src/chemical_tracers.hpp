#ifndef CHEMICAL_TRACERS_HPP
#define CHEMICAL_TRACERS_HPP

#include <athena_arrays.hpp>
#include "solution_group.hpp"

class MeshBlock;
class ParameterInput;

class ChemicalTracers: public SolutionGroup {
public:
// public access variables

  // conserved variables (extensive), kg/m^3
  AthenaArray<Real> c;

  // primitive variables (intensive), kg/kg
  AthenaArray<Real> q;

  // flux in 3-dimensions
  AthenaArray<Real> flux[3];

// functions
  ChemicalTracers(MeshBlock *pmb, ParameterInput *pin);

protected:
  PassiveScalars *pscalars_;
};

#endif
