#ifndef VAPOR_HYDRO_HPP
#define VAPOR_HYDRO_HPP

#include "vapor.hpp"
#include <configure.hpp>

template<typename Base>
template<typename T>
inline Real Vapor<Base>::getGamma(T w) const {
  Real gammad = getDryGamma();
  Real fsig = 1., feps = 1.;
  for (int n = 1; n <= NumVapors; ++n) {
    fsig += w[n]*(m_cv_ratios_mass[n] - 1.);
    feps += w[n]*(1./m_molwt_ratios[n] - 1.);
  }
  return 1. + (gammad - 1.)*feps/fsig;
}

template<typename Base>
template<typename T>
inline Real Vapor<Base>::getMeanMolecularWeight(T w) const {
  Real feps = 1.;
  for (int n = 1; n <= NumVapors; ++n)
    feps += w[n]*(m_molwt_ratios[n] - 1.);
  return feps*m_dry_mmw;
}

template<typename Base>
template<typename T>
inline Real Vapor<Base>::getTemp(T w) const {
  return w[ipr]*getMeanMolecularWeight(w)/(w[idn]*Cantera::GasConstant);
}

template<typename Base>
template<typename T>
inline Real Vapor<Base>::getPres(T u) const {
  Real gm1 = getDryGamma() - 1.;
  Real rho = 0., fsig = 0., feps = 0.;
  for (int n = 0; n <= NumVapors; ++n) {
    rho += u[n];
    fsig += u[n]*m_cv_ratios_mass[n];
    feps += u[n]/m_molwt_ratios[n];
  }
  Real KE = 0.5*(u[iv1]*u[iv1] + u[iv2]*u[iv2] + u[iv3]*u[iv3])/rho;
  return gm1*(u[ien] - KE)*feps/fsig;
}

template<typename Base>
template<typename T>
inline Real Vapor<Base>::getRovCp(T w) const {
  Real gammad = getDryGamma();
  Real qsig = 1., feps = 1.;
  for (int n = 1; n <= NumVapors; ++n) {
    qsig += w[n]*(m_cp_ratios_mass[n] - 1.);
    feps += w[n]*(1./m_molwt_ratios[n] - 1.);
  }
  return (gammad - 1.)/gammad*feps/qsig;
}

//! c_p = c_{pd}*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
//!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
template<typename Base>
template<typename T>
inline Real Vapor<Base>::cp_JouleKilogramK(T w) const {
  Real gammad = getDryGamma();
  Real qsig = 1.;
  for (int n = 1; n <= NumVapors; ++n)
    qsig += w[n]*(m_cp_ratios_mass[n] - 1.);
  return gammad/(gammad - 1.)*qsig*Cantera::GasConstant/m_dry_mmw;
}

//! c_v = c_{vd}*(1 + \sum_i (q_i*(\hat{c}_{vi} - 1.)))
//!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{vi} - 1.)))
template<typename Base>
template<typename T>
inline Real Vapor<Base>::cv_JouleKilogramK(T w) const {
  Real gammad = getDryGamma();
  Real qsig = 1.;
  for (int n = 1; n <= NumVapors; ++n)
    qsig += w[n]*(m_cv_ratios_mass[n] - 1.);
  return 1./(gammad - 1.)*qsig*Cantera::GasConstant/m_dry_mmw;
}

/*! h = c_{pd}*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
//!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
template<typename T>
inline Real Vapor::getSpecificEnthalpy(T w) const {
  Real gamma = peos_->getGamma();
  Real tem[1] = {getTemp(w)};
  update_gamma(gamma, tem);
  Real qsig = 1.;
  for (int n = 1; n <= NumVapors; ++n)
    qsig += w[n]*(thermo_[n].cp_ratio - 1.);
  return gamma/(gamma - 1.)*Rd_*qsig*tem[0];
}*/

#endif
