// C/C++
#include <iostream>
#include <algorithm>

#include <molecule/molecule.hpp>
#include "construct_atmosphere.hpp"
#include "gasphase.hpp"
#include "dlnT_dlnP.hpp"
#include "condense_vapor.hpp"

int const max_iter_ = 10;

void rk4_integrate_z(Real q[], int isat[], GasThermo const thermo[],
  Molecule const moist[], Real gamma, Real g_ov_Rd, Real dz, int method, Real adTdz)
{
  Real step[] = {0.5, 0.5, 1.};
  Real temp = q[idn];
  Real pres = q[ipr];
  Real dTdz[4], chi[4];

  GasThermo thermo0[1+NumVapors];

  for (int rk = 0; rk < 4; ++rk) {
    // calculate gamma
    update_gamma(gamma, q);

    Real q_eps = 1.;
    for (int iv = 1; iv <= NumVapors; ++iv)
      q_eps += q[iv]*(thermo[iv].mu_ratio - 1.);
    Real R_ov_Rd = 1./q_eps;

    if (method == 0 || method == 1)
      chi[rk] = dlnTdlnP(q, isat, thermo, gamma);
    else if (method == 2)
      chi[rk] = dlnTdlnP(q, isat, thermo, gamma);
    else  // isothermal
      chi[rk] = 0.;
    dTdz[rk] = - chi[rk]*g_ov_Rd/R_ov_Rd + adTdz;
    chi[rk] = - R_ov_Rd/g_ov_Rd*dTdz[rk];
    
    // integrate over dz
    Real chi_avg;
    if (rk < 3) {
      q[idn] = temp + dTdz[rk]*dz*step[rk];
      chi_avg = chi[rk];
    } else {
      q[idn] = temp + 1./6.*(dTdz[0] + 2.*dTdz[1] + 2.*dTdz[2] + dTdz[3])*dz;
      chi_avg = 1./6.*(chi[0] + 2.*chi[1] + 2.*chi[2] + chi[3]);
    }
    if (!(q[idn] > 0.)) q[idn] = temp;
    if (fabs(q[idn] - temp) > 0.1) // isothermal limit
      q[ipr] = pres*pow(q[idn]/temp, 1./chi_avg);
    else
      q[ipr] = pres*exp(-2.*g_ov_Rd*dz/(R_ov_Rd*(q[idn] + temp)));
  }

  Real sum;
  do {
    // condensation
    sum = 1.;
    for (int iv = 1; iv <= NumVapors; ++iv) {
      Real rate = condense_vapor(q, iv, moist[iv]);
      if (rate > 0.) isat[iv] = 1;
      q[iv] -= rate;
    }

    // normalization
    for (int iv = 1; iv <= NumVapors; ++iv)
      q[iv] /= sum;
  } while (std::abs(sum - 1.) < 1.E-10);
}

void rk4_integrate_z_adaptive(Real q[], int isat[], GasThermo const thermo[], 
  Molecule const moist[], Real gamma, Real g_ov_Rd, Real dz,
  Real ftol, int &niter, int method, Real adTdz)
{
  Real q1[NumHydros], q2[NumHydros];
  int  isat1[1+NumVapors], isat2[1+NumVapors];

  for (int n = 0; n < NumHydros; ++n) {
    q1[n] = q[n];
    q2[n] = q[n];
  }
  for (int n = 0; n <= NumVapors; ++n) {
    isat1[n] = isat[n];
    isat2[n] = isat[n];
  }

  // trail step
  rk4_integrate_z(q1, isat1, thermo, moist, gamma, g_ov_Rd, dz, method, adTdz);

  // refined step
  rk4_integrate_z(q2, isat2, thermo, moist, gamma, g_ov_Rd, dz/2., method, adTdz);
  rk4_integrate_z(q2, isat2, thermo, moist, gamma, g_ov_Rd, dz/2., method, adTdz);

  // abort if dz is less than 0.001 m
  //std::cout << dz << " ";
  //std::cout << fabs(q2[idn] - q1[idn]) << " ";
  //std::cout << ftol << std::endl;
  //assert(dz > 0.001);

  if ((fabs(q2[idn] - q1[idn]) > ftol) &&
      (niter++ < max_iter_))
  {
    rk4_integrate_z_adaptive(q, isat, thermo, moist, gamma,
                             g_ov_Rd, dz/2., ftol, niter, method, adTdz);
    rk4_integrate_z_adaptive(q, isat, thermo, moist, gamma,
                             g_ov_Rd, dz/2., ftol, niter, method, adTdz);
  } else {
    for (int n = 0; n < NumHydros; ++n)
      q[n] = q2[n];
    for (int n = 0; n <= NumVapors; ++n)
      isat[n] = isat2[n];
  }
}
