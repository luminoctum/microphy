// C/C++ headers
#include <iostream>
#include <algorithm>

// debugger
#include <debugger.hpp>

#include "dlnT_dlnP.hpp"

//int const max_iter_ = 10;

template<typename VaporPhase, typename CloudPhase>
void rk4_integrate_lnp(
  AtmConstructionMethod method,
  Vapor<VaporPhase> *vapor,
  Suspension<VaporPhase, CloudPhase> *cloud,
  Cantera::Kinetics *kin)
{
  Real step[] = {0.5, 0.5, 1.};
  Real temp = vapor->temperature();
  Real pres = vapor->pressure();
  Real chi[4];

  //DimensionlessThermo thermo0[1+NumVapors];

  for (int rk = 0; rk < 4; ++rk) {
    vapor->updateThermo();
    // calculate gamma
    // update_gamma(gamma, q);

    // calculate tendency
    if (method.type == "reversible" || method.type == "pseudo") {
      chi[rk] = dlnTdlnP(vapor, cloud, kin);
    } else if (method.type == "dry") {
      chi[rk] = dlnTdlnP(vapor, cloud, kin);
    } else if (method.type == "isothermal") {
      chi[rk] = 0.;
    } else {
      Debugger::Fatal("rk4_integrate_lnp", "Unrecognized method: " + method.type);
    }
    chi[rk] *= method.rdlnTdlnP;
    
    // integrate over dlnp
    Real dlnp = method.dlnp;
    if (rk < 3) {
      vapor->setState_TP(temp*exp(chi[rk]*dlnp*step[rk]), pres*exp(dlnp));
    } else {
      vapor->setState_TP(temp*exp(1./6.*(chi[0] + 2.*chi[1] + 2.*chi[2] + chi[3])*dlnp),
        pres*exp(dlnp));
    }
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

/*void rk4_integrate_lnp_adaptive(Real q[], int isat[], DimensionlessThermo const thermo[],
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
}*/
