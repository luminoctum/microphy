#ifndef DLNT_DLNP_HPP
#define DLNT_DLNP_HPP

#include <configure.hpp>
#include "gas_thermo.hpp"

/*! Moist adiabatic lapse rate
 *
 * Calculate moist adiabatic lapse rate according to @cite li2018
 * @param q an array of [T,q1,q2,...,qn,P,q11,q12,q21,q22,...,qn1,qn2]
 * @param isat indicator array of saturated or not
 * @param rcp molar cp ratios with respect to dry air
 */
Real dlnTdlnP(Real const q[], int const isat[], GasThermo const thermo[], Real gamma)
{
  Real f_sig = 1.;
  // mass cp ratio to molar cp ratio
  for (int n = 1; n <= NumVapors; ++n)
    f_sig += q[n]*(thermo[n].cp_ratio*thermo[n].mu_ratio - 1.);
  Real cphat_ov_r = gamma/(gamma - 1.)*f_sig;

  // vapor
  Real xd = 1.;
  for (int n = 1; n <= NumVapors; ++n) xd -= q[n];

  Real c1 = 0., c2 = 0., c3 = 0.;
  for (int iv = 1; iv <= NumVapors; ++iv) {
    if (isat[iv] > 0) {
      Real latent = thermo[iv].beta*thermo[iv].t3/q[idn] - thermo[iv].delta;
      c1 += q[iv]/xd*latent;
      c2 += q[iv]/xd*latent*latent;
    }
    c3 += q[iv]/xd;
  }

  return (1. + c1)/(cphat_ov_r + (c2 + c1*c1)/(1. + c3));
}

#endif
