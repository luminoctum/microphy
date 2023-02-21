#include <iostream>
#include <mesh/mesh.hpp>
#include <parameter_input.hpp>
#include <scalars/scalars.hpp>
#include <debugger.hpp>
#include <configure.hpp>
#include "condensates.hpp"

int Globals::next_scalar_index = 0;

template<int NumTrackingFields>
Condensates<NumTrackingFields>::Condensates(MeshBlock *pmb, ParameterInput *pin):
  SolutionGroup(
    pin->GetString("chemistry", "kinetics_file"),
    pin->GetString("chemistry", "condensates")
    ),
  pscalars_(pmb->pscalars)
{
  std::stringstream msg;
  /*if (!pin->DoesParameterExist("chemistry", "condensates")) {
    pdebug->Leave();
    return;
  }*/

  // den, tem
  Globals::next_scalar_index += NumTrackingFields*nTotalSpecies();

  if (Globals::next_scalar_index > NumScalars) {
    msg << Globals::next_scalar_index << " scalars needed, but "
        << NumScalars << " configured. Increase the number of scalars";
    Debugger::Fatal("Condensates", msg.str());
  }

  q.InitWithShallowSlice(pscalars_->r, 4, scalarIndexFirst(), nTotalSpecies());
  c.InitWithShallowSlice(pscalars_->s, 4, scalarIndexFirst(), nTotalSpecies());

  if (NumTrackingFields > 1) {
    tem.InitWithShallowSlice(pscalars_->r, 4, scalarIndexFirst() + nTotalSpecies(), 
      nTotalSpecies());
    energy.InitWithShallowSlice(pscalars_->s, 4, scalarIndexFirst() + nTotalSpecies(),
      nTotalSpecies());
  }

  //pdebug->Leave();
}

// explicit instantiations
template class Condensates<1>;
