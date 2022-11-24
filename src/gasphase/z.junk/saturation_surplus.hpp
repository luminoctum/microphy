#ifndef SATURATION_SURPLUS_HPP
#define SATURATION_SURPLUD_HPP

// canoe headers
#include <configure.hpp>

/*! Saturation surplus for vapors can be both positive and negative
 * positive value represents supersaturation \n
 * negative value represents saturation deficit
 */
template<typename T>
inline void Thermodynamics::SaturationSurplus(Real dv[], T v, VariableType vtype) const {
  // mass to molar mixing ratio
  if (vtype == VariableType::prim) {
    PrimitiveToChemical(q1_, v);
  } else if (vtype == VariableType::cons) {
    ConservedToChemical(q1_, v);
  } else { // VariableType::chem
    for (int n = 0; n < NumHydros; ++n)
      q1_[n] = v[n];
  }
  // change molar density to molar mixing ratio
  Real mols = q1_[ipr]/(q1_[idn]*Rgas);
  for (int n = 1; n <= NumVapors; ++n) q1_[n] /= mols;

  for (int iv = 1; iv <= NumVapors; ++iv) { 
    int nc = q1_[idn] > t3_[iv] ? iv + NumVapors : iv + 2*NumVapors;
    int ic = NumHydros - NumVapors + nc - 1;
    Real rate = VaporCloudEquilibrium(q1_, iv, ic, t3_[iv], p3_[iv], 
        0., beta_[nc], delta_[nc], true);
    dv[iv] = rate/q1_[iv]*v[iv];
  }
}

#endif
