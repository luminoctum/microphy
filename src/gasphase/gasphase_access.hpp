#ifndef GASPHASE_ACCESS_HPP
#define GASPHASE_ACCESS_HPP

#include "gasphase.hpp"
#include "../constants.hpp"

template<typename T>
inline Real GasPhase::getGamma(T w) const {
  Real gamma = peos_->getGamma();
  Real fsig = 1., feps = 1.;
  for (int n = 1; n <= NumVapors; ++n) {
    fsig += w[n]*(thermo_[n].cv_ratio - 1.);
    feps += w[n]*(1./thermo_[n].mu_ratio - 1.);
  }
  return 1. + (gamma - 1.)*feps/fsig;
}

template<typename T>
inline Real GasPhase::RovRd(T w) const {
  Real feps = 1.;
  for (int n = 1; n <= NumVapors; ++n)
    feps += w[n]*(1./thermo_[n].mu_ratio - 1.);
  return feps;
}

//! Temperature
template<typename T>
inline Real GasPhase::getTemp(T w) const {
  return w[ipr]/(w[idn]*Rd_*RovRd(w));
}

template<typename T>
inline Real GasPhase::getPres(T u) const {
  Real gm1 = peos_->getGamma() - 1;
  Real rho = 0., fsig = 0., feps = 0.;
  for (int n = 0; n <= NumVapors; ++n) {
    rho += u[n];
    fsig += u[n]*thermo_[n].cv_ratio;
    feps += u[n]/thermo_[n].mu_ratio;
  }
  Real KE = 0.5*(u[iv1]*u[iv1] + u[iv2]*u[iv2] + u[iv3]*u[iv3])/rho;
  return gm1*(u[ien] - KE)*feps/fsig;
}

template<typename T>
inline Real GasPhase::getChi(T w) const {
  Real gamma = peos_->getGamma();
  Real tem[1] = {getTemp(w)};
  update_gamma(gamma, tem);
  Real qsig = 1., feps = 1.;
  for (int n = 1; n <= NumVapors; ++n) {
    qsig += w[n]*(thermo_[n].cp_ratio - 1.);
    feps += w[n]*(1./thermo_[n].mu_ratio - 1.);
  }
  return (gamma - 1.)/gamma*feps/qsig;
}

//! c_p = c_{pd}*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
//!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
template<typename T>
inline Real GasPhase::getSpecificCp(T w) const {
  Real gamma = peos_->getGamma();
  Real tem[1] = {getTemp(w)};
  update_gamma(gamma, tem);
  Real qsig = 1.;
  for (int n = 1; n <= NumVapors; ++n)
    qsig += w[n]*(thermo_[n].cp_ratio - 1.);
  return gamma/(gamma - 1.)*Rd_*qsig;
}

//! c_v = c_{vd}*(1 + \sum_i (q_i*(\hat{c}_{vi} - 1.)))
//!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{vi} - 1.)))
template<typename T>
inline Real GasPhase::getSpecificCv(T w) const {
  Real gamma = peos_->getGamma();
  Real tem[1] = {getTemp(w)};
  update_gamma(gamma, tem);
  Real qsig = 1.;
  for (int n = 1; n <= NumVapors; ++n)
    qsig += w[n]*(thermo_[n].cv_ratio - 1.);
  return 1./(gamma - 1.)*Rd_*qsig;
}

//! h = c_{pd}*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
//!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
template<typename T>
inline Real GasPhase::getSpecificEnthalpy(T w) const {
  Real gamma = peos_->getGamma();
  Real tem[1] = {getTemp(w)};
  update_gamma(gamma, tem);
  Real qsig = 1.;
  for (int n = 1; n <= NumVapors; ++n)
    qsig += w[n]*(thermo_[n].cp_ratio - 1.);
  return gamma/(gamma - 1.)*Rd_*qsig*tem[0];
}

template<typename T>
inline Real GasPhase::getMeanMolecularWeight(T w) const {
  Real feps = 1.;
  for (int n = 1; n <= NumVapors; ++n)
    feps += w[n]*(thermo_[n].mu_ratio - 1.);
  return Constants::Rgas/Rd_*feps;
}

#endif
