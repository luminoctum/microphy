/** @file constuct_atmosphere.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Wednesday May 26, 2021 18:45:55 UTC
 * @bug No known bugs.
 */

// C/C++ headers
#include <iostream>
#include <cstdlib>

// Athena++ headers
//#include <athena_arrays.hpp>
//#include <hydro/hydro.hpp>

// canoe
#include <configure.hpp>
#include "gasphase.hpp"
#include "construct_atmosphere.hpp"
#include "condense_vapor.hpp"
#include "../constants.hpp"

void GasPhase::ConstructAtmosphere(Real **w, Real Ts, Real Ps,
  Real grav, Real dzORdlnp, int len, Adiabat method, Real userp) const
{
  Real gamma = peos_->getGamma();

  Real q1[NumHydros];
  // hydro + vapor + condensate (0)
  w[0][idn] = w[0][ipr] = 1.;
  PrimitiveToChemical(q1, w[0]);

  // change molar concentration to mixing ratio
  Real mols = q1[ipr]/(q1[idn]*Constants::Rgas);
  for (int n = 1; n <= NumVapors; ++n) q1[n] /= mols;
  // reset TP
  q1[idn] = Ts;
  q1[ipr] = Ps;

  int isat[1+NumVapors];
  for (int n = 0; n <= NumVapors; ++n)
    isat[n] = 0;

  Real sum;
  do {
    sum = 1.;
    // condensation
    for (int iv = 1; iv <= NumVapors; ++iv) {
      Real rate = condense_vapor(q1, iv, moist_[iv]);
      if (rate > 0.) isat[iv] = 1;
      q1[iv] -= rate;
      sum -= rate;
    }

    // normalization
    for (int iv = 1; iv <= NumVapors; ++iv)
      q1[iv] /= sum;
  } while (std::abs(sum - 1.) < 1.E-10);

  // set vapor
  ChemicalToPrimitive(w[0], q1);

  for (int i = 1; i < len; ++i) {
    // RK4 integration 
    Real rdlnTdlnP, adTdz;
    int niter = 0;
#if HYDROSTATIC
    rdlnTdlnP = userp;
    rk4_integrate_lnp_adaptive(q1, isat, thermo_, gamma,
      dzORdlnp, ftol_, niter, (int)method, rdlnTdlnP);
#else
    adTdz = userp;
    rk4_integrate_z_adaptive(q1, isat, thermo_, gamma,
      grav/Rd_, dzORdlnp, ftol_, niter, (int)method, adTdz);
#endif

    // set vapor
    ChemicalToPrimitive(w[i], q1);
  }
}
