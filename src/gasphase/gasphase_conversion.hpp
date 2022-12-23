#ifndef GASPHASE_CONVERSION_HPP
#define GASPHASE_CONVERSION_HPP

#include <configure.hpp>
#include "gas_phase.hpp"
#include "../constants.hpp"

//! Change mass mixing ratio to molar mixing ratio
template<typename T1, typename T2>
inline void GasPhase::PrimitiveToChemical(T1 qmol, T2 const w) const
{
  // set molar mixing ratio
  Real sum = 1.;
  for (int n = 1; n <= NumVapors; ++n) {
    qmol[n] = w[n]/thermo_[n].mu_ratio;
    sum += w[n]*(1./thermo_[n].mu_ratio - 1.);
  }
  // set pressure, temperature, velocity
  qmol[ipr] = w[ipr];
  qmol[idn] = w[ipr]/(w[idn]*Rd_*sum);
  qmol[iv1] = w[iv1];
  qmol[iv2] = w[iv2];
  qmol[iv3] = w[iv3];

#pragma omp simd
  for (int n = 1; n <= NumVapors; ++n) {
    qmol[n] /= sum;
  }
}

//! Change molar mixing ratio to mass mixing ratio
template<typename T1, typename T2>
inline void GasPhase::ChemicalToPrimitive(T1 w, T2 const qmol) const {
  // set mass mixing ratio
  Real sum = 1.;
  for (int n = 1; n <= NumVapors; ++n) {
    w[n] = qmol[n]*thermo_[n].mu_ratio;
    sum += qmol[n]*(thermo_[n].mu_ratio - 1.); 
  }
#pragma omp simd
  for (int n = 1; n <= NumVapors; ++n)
    w[n] /= sum;

  // set pressure, density, velocity
  w[ipr] = qmol[ipr];
  w[idn] = sum*qmol[ipr]/(qmol[idn]*Rd_);
  w[iv1] = qmol[iv1];
  w[iv2] = qmol[iv2];
  w[iv3] = qmol[iv3];
}

//! Change density to molar mixing ratio
template<typename T1, typename T2>
inline void GasPhase::ConservedToChemical(T1 qmol, T2 const u) const {
  Real rho = 0., feps = 0., fsig = 0.;
  for (int n = 0; n <= NumVapors; ++n) {
    rho += u[n];
    qmol[n] = u[n]/thermo_[n].mu_ratio;
    feps += qmol[n];
    fsig += u[n]*thermo_[n].cv_ratio;
  }
  Real KE = 0.5*(u[iv1]*u[iv1] + u[iv2]*u[iv2] + u[iv3]*u[iv3])/rho;
  Real gm1 = peos_->getGamma() - 1.;
  qmol[ipr] = gm1*(u[ien] - KE)*feps/fsig;
  qmol[idn] = qmol[ipr]/(feps*Rd_);
  qmol[iv1] = u[iv1]/rho;
  qmol[iv2] = u[iv2]/rho;
  qmol[iv3] = u[iv3]/rho;

  Real mols = qmol[ipr]/(Constants::Rgas*qmol[idn]);
#pragma omp simd
  for (int n = 1; n <= NumVapors; ++n)
    qmol[n] *= mols/feps;
}

//! Change molar mixing ratio to density
template<typename T1, typename T2>
inline void GasPhase::ChemicalToConserved(T1 u, T2 const qmol) const {
  Real sum = 1., mols = qmol[ipr]/(Constants::Rgas*qmol[idn]);
  for (int n = 1; n <= NumVapors; ++n)
    sum += qmol[n]/mols*(thermo_[n].mu_ratio - 1.);
  Real rho = qmol[ipr]*sum/(Rd_*qmol[idn]);
  Real cvd = Rd_/(peos_->getGamma() - 1.);
  u[idn] = rho;
  u[ien] = 0.5*rho*(qmol[iv1]*qmol[iv1] 
    + qmol[iv2]*qmol[iv2] + qmol[iv3]*qmol[iv3]);
  for (int n = 1; n <= NumVapors; ++n) {
    u[n] = rho*qmol[n]/mols*thermo_[n].mu_ratio/sum;
    u[idn] -= u[n];
    u[ien] += u[n]*thermo_[n].cv_ratio*cvd*qmol[idn];
  }
  u[ien] += u[idn]*cvd*qmol[idn];
  u[iv1] = qmol[iv1]*rho;
  u[iv2] = qmol[iv2]*rho;
  u[iv3] = qmol[iv3]*rho;
}


#endif
