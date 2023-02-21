// C/C++ headers
#include <iostream>
#include <algorithm>

// debugger
#include <debugger.hpp>

#include "dlnT_dlnP.hpp"

//int const max_iter_ = 10;

template<typename VaporPhase, typename CloudPhase>
void rk4_integrate_z(
  AtmConstructionMethod method,
  Vapor<VaporPhase> *vapor,
  Suspension<VaporPhase, CloudPhase> *cloud,
  Cantera::Kinetics *kin)
{
  Real step[] = {0.5, 0.5, 1.};
  Real temp = vapor->temperature();
  Real pres = vapor->pressure();
  Real dTdz[4], chi[4];

  for (int rk = 0; rk < 4; ++rk) {
    vapor->updateThermo();
    // calculate gamma
    //update_gamma(gamma, q);

    Real mug_ov_R = vapor->meanMolecularWeight()*method.grav/Cantera::GasConstant;

    if (method.type == "reversible" || method.type == "pseudo") {
      chi[rk] = dlnTdlnP(vapor, cloud, kin);
    } else if (method.type == "dry") {
      chi[rk] = dlnTdlnP(vapor, cloud, kin);
    } else if (method.type == "isothermal") {// isothermal
      chi[rk] = 0.;
    } else {
      Debugger::Fatal("rk4_integrate_z", "Unrecognized method: " + method.type);
    }
    dTdz[rk] = - chi[rk]*mug_ov_R + method.adTdz;
    chi[rk] = - dTdz[rk]/mug_ov_R;
    
    // integrate over dz
    Real chi_avg, temp2, pres2;
    Real dz = method.dz;
    if (rk < 3) {
      temp2 = temp + dTdz[rk]*dz*step[rk];
      chi_avg = chi[rk];
    } else {
      temp2 = temp + 1./6.*(dTdz[0] + 2.*dTdz[1] + 2.*dTdz[2] + dTdz[3])*dz;
      chi_avg = 1./6.*(chi[0] + 2.*chi[1] + 2.*chi[2] + chi[3]);
    }
    if (!(temp2 > 0.)) temp2 = temp;
    if (fabs(temp2 - temp) > 0.1) // isothermal limit
      pres2 = pres*pow(temp2/temp, 1./chi_avg);
    else
      pres2 = pres*exp(-2.*mug_ov_R*dz/(temp2 + temp));
    vapor->setState_TP(temp2, pres2);
  }

  std::vector<Real> conc(vapor->nSpecies());
  std::vector<Real> act_conc(kin->nTotalSpecies());
  vapor->getConcentrations(conc.data());
  kin->getActivityConcentrations(act_conc.data());

  // remove supersaturated vapors
  for (int i = 1; i <= NumVapors; ++i) {
    std::string name = vapor->vaporName(i);
    conc[vapor->speciesIndex(name)] -= act_conc[kin->kineticsSpeciesIndex(name)];
  }
  vapor->setConcentrations(conc.data());
}

/*void rk4_integrate_z_adaptive(Real q[], int isat[], DimensionlessThermo const thermo[], 
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
}*/
