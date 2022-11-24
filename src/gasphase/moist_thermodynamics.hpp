/** @file thermodynamic_funcs.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Tuesday Jun 01, 2021 15:42:38 PDT
 * @bug No known bugs.
 */

#ifndef THERMODYNAMICS_FUNCS_HPP
#define THERMODYNAMICS_FUNCS_HPP

#include "thermodynamics.hpp"
//#include "../particles/particles.hpp"

//! Potential temperature
template<typename T>
Real PotentialTemp(T w, Real p0, Thermodynamics const *pthermo) {
  Real chi = pthermo->getChi(w);
  Real temp = pthermo->getTemp(w);
  return temp*pow(p0/w[IPR], chi);
}

//! Moist static energy
template<typename T>
Real MoistStaticEnergy(T w, Real gz, Thermodynamics const *pthermo) {
  Real temp = pthermo->getTemp(w);
  Real IE = w[IDN]*pthermo->getSpecificCp(w)*temp;
  Real rho = w[IDN];
  /*if (ppart != nullptr) {
    for (int n = 0; n < NumVapors; ++n) {
      for (int t = 0; t < ppart->u.GetDim4(); ++t) {
        rho += ppart->u(t,k,j,i);
        IE -= ppart->u(t,k,j,i)*pthermo->GetLatent(1+NumVapors+n);
        IE += ppart->u(t,k,j,i)*ppart->GetCv(t)*temp;
      }
      ppart = ppart->next;
    }
  }*/
  return IE/rho + gz;
}

//! Relative humidity
template<typename T>
Real RelativeHumidity(T w, int iv, Thermodynamics const *pthermo) {
  Real dw[1+NumVapors];
  pthermo->SaturationSurplus(dw, w, VariableType::prim);
  return w[iv]/(w[iv] - dw[iv]);
}

//! Liquid water potential temperature
// approximate formula for liquid water potential temperature is
// theta_l = theta - theta/T*Lv/cpd*ql
template<typename T>
Real LiquidWaterPotentialTemp(T w, Real p0, Real lv_ov_cpd, Real qliq,
  Thermodynamics const *pthermo) {
  Real theta = PotentialTemp(w, p0, pthermo);
  Real temp = pthermo->getTemp(w);
  return theta - theta/temp*lv_ov_cpd*qliq;
}

//! Equivalent potential temperature
/*Real MoistEntropy(AthenaArray<Real> const& w, Thermodynamics *pthermo, Particles *ppart,
  int k, int j, int i) {
#if (NumVapors > 0)
  Real gamma = pthermo->pmy_block->peos->getGamma();
  Real tem[1] = {pthermo->getTemp(prim)};
  update_gamma(gamma, tem);
  Real cpd = Rd_*gamma/(gamma - 1.);
  Real temp = tem[0]
  Real pres = prim[IPR];

  Real qd = 1.;
  for (int n = 1; n <= NumVapors; ++n)
    qd -= prim[n];

  Real qc[1+NumVapors];
  std::fill(qc, qc + 1 + NumVapors, 0.);
  for (int n = 1 + NumVapors; n <= NumVapors; ++n)
    qc[1+(n-1)%NumVapors] += prim[n] + 1.0E-10;  // prevent devide by 0

  Real lv = 0.;
  for (int n = 1 + NumVapors; n < NMASS; ++n) {
    int ng = 1 + (n-1)%NumVapors;
    Real ratio = (prim[n] + 1.0E-10)/qc[ng];
    lv += GetLatent(n,temp)*prim[ng]*ratio;
  }

  Real st = 1.;
  for (int n = 1 + NumVapors; n < NMASS; ++n) {
    int ng = 1 + (n-1)%NumVapors;
    Real ratio = (prim[n] + 1.0E-10)/qc[ng];
    st += (prim[n] + prim[ng]*ratio)*(getCpRatio(n) - 1.);
  }
  Real lv_ov_cpt = lv/(cpd*st*temp);

  Real chi = Rd_/cpd*qd/st;

  Real xv = 1.;
  for (int n = 1; n <= NumVapors; ++n)
    xv += prim[n]/qd/mu_ratios_[n];

  Real pd = pres/xv;

  Real rh = 1.;
  for (int n = 1; n <= NumVapors; ++n) {
    Real eta = prim[n]/qd/mu_ratios_[n];
    Real pv = pres*eta/xv;
    int nc = n + NumVapors;
    Real esat;
    if (n == AMMONIA_VAPOR_ID)
      esat = sat_vapor_p_NH3_BriggsS(temp);
    else if (n == WATER_VAPOR_ID)
      esat = sat_vapor_p_H2O_BriggsS(temp);
    else
      esat = SatVaporPresIdeal(temp/t3_[nc], p3_[nc], beta_[nc], delta_[nc]);
    rh *= pow(pv/esat, -eta*Rd_/(cpd*st));
  }

  return temp*pow(p0/pd, chi)*exp(lv_ov_cpt)*rh;
#else
  return PotentialTemp(prim, p0, pthermo);
#endif
}*/

Real saha_ionization_electron_density(Real T, Real num, Real ion_ev);


#endif /* end of include guard THERMODYNAMICS_FUNCS_HPP */

