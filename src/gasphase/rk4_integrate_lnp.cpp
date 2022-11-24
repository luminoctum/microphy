#include <iostream>
#include <algorithm>

#include <molecule/molecule.hpp>
#include "construct_atmosphere.hpp"
#include "gasphase.hpp"
#include "dlnT_dlnP.hpp"
#include "condense_vapor.hpp"

int const max_iter_ = 10;

void rk4_integrate_lnp(Real q[], int isat[], GasThermo const thermo[], 
  Molecule const moist[], Real gamma, Real dlnp, int method, Real rdlnTdlnP)
{
  Real step[] = {0.5, 0.5, 1.};
  Real temp = q[idn];
  Real pres = q[ipr];
  Real chi[4];

  GasThermo thermo0[1+NumVapors];

  for (int rk = 0; rk < 4; ++rk) {
    // calculate gamma
    update_gamma(gamma, q);

    // calculate tendency
    if (method == 0 || method == 1)
      chi[rk] = dlnTdlnP(q, isat, thermo, gamma);
    else if (method == 2)
      chi[rk] = dlnTdlnP(q, isat, thermo0, gamma);
    else // isothermal
      chi[rk] = 0.;
    chi[rk] *= rdlnTdlnP;
    
    // integrate over dlnp
    if (rk < 3)
      q[idn] = temp*exp(chi[rk]*dlnp*step[rk]);
    else
      q[idn] = temp*exp(1./6.*(chi[0] + 2.*chi[1] + 2.*chi[2] + chi[3])*dlnp);
    q[ipr] = pres*exp(dlnp);
  }

  Real sum;
  do {
    // condensation
    sum = 1.;
    for (int iv = 1; iv <= NumVapors; ++iv) {
      Real rate = condense_vapor(q, iv, moist[iv]);
      if (rate > 0.) isat[iv] = 1;
      q[iv] -= rate;
      sum -= rate;
    }

    // normalization
    for (int iv = 1; iv <= NumVapors; ++iv)
      q[iv] /= sum;
  } while (std::abs(sum - 1.) < 1.E-10);
}

void rk4_integrate_lnp_adaptive(Real q[], int isat[], GasThermo const thermo[],
  Molecule const moist[], Real gamma, Real dlnp,
  Real ftol, int &niter, int method, Real rdlnTdlnP)
{
  Real q1[NumHydros], q2[NumHydros];
  int  isat1[1+NumVapors], isat2[1+NumVapors];

  for (int n = 0; n < NumHydros; ++n) {
    q1[n] = q[n];
    q2[n] = q[n];
  }
  for (int n = 0; n <=  NumVapors; ++n) {
    isat1[n] = isat[n];
    isat2[n] = isat[n];
  }

  // trail step
  rk4_integrate_lnp(q1, isat1, thermo, moist, gamma, dlnp, method, rdlnTdlnP);

  // refined step
  rk4_integrate_lnp(q2, isat2, thermo, moist, gamma, dlnp/2., method, rdlnTdlnP);
  rk4_integrate_lnp(q2, isat2, thermo, moist, gamma, dlnp/2., method, rdlnTdlnP);

  if ((fabs(q2[idn] - q1[idn]) > ftol) &&
      (niter++ < max_iter_))
  {
    rk4_integrate_lnp_adaptive(q, isat, thermo, moist, gamma, dlnp/2., ftol/2., 
      niter, method, rdlnTdlnP);
    rk4_integrate_lnp_adaptive(q, isat, thermo, moist, gamma, dlnp/2., ftol/2.,
      niter, method, rdlnTdlnP);
  } else {
    for (int n = 0; n < NumHydros; ++n)
      q[n] = q2[n];
    for (int n = 0; n <= NumVapors; ++n)
      isat[n] = isat2[n];
  }
}
