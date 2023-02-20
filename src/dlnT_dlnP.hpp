#ifndef DLNT_DLNP_HPP
#define DLNT_DLNP_HPP

#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/thermo/IdealSolidSolnPhase.h>
#include <cantera/kinetics.h>

#include "vapor.hpp"

struct AtmConstructionMethod {
  std::string type = "pseudo";
  Real grav = 0.;
  Real dlnp = 0.;
  Real dz = 0.;
  Real rdlnTdlnP = 1.;
  Real adTdz = 0.;
};

template<typename VaporPhase, typename CloudPhase>
Real dlnTdlnP(
  Vapor<VaporPhase> *vapor,
  Suspension<VaporPhase, CloudPhase> *cloud,
  Cantera::Kinetics *kin
  );

/** Moist adiabatic lapse rate
 *
 * Calculate moist adiabatic lapse rate according to @cite li2018
 * @param q an array of [T,q1,q2,...,qn,P,q11,q12,q21,q22,...,qn1,qn2]
 * @param rcp molar cp ratios with respect to dry air
 */

// speciealization
template<>
Real dlnTdlnP(
  Vapor<Cantera::IdealGasPhase> *vapor,
  Suspension<Cantera::IdealGasPhase, Cantera::IdealSolidSolnPhase> *cloud,
  Cantera::Kinetics *kin
  )
{
  Real cp_gas = vapor->cp_mole();
  Real cp_cloud = cloud->cp_mole()*cloud->getRelativeMole();
  Real cphat_ov_r = (cp_gas + cp_cloud)/Cantera::GasConstant;

  // get mole fractions
  std::vector<Real> mole_fractions;
  vapor->getMoleFractions(mole_fractions.data());

  Real xdry = 0.;
  for (int i = 0; i < vapor->nDryMixtures(); ++i) xdry += mole_fractions[i];

  Real c1 = 0., c2 = 0., c3 = 0.;
  vapor->setLatentHeats(kin);

  for (int i = 1; i <= NumVapors; ++i) {
    Real temperature = vapor->temperature();
    Real eta = mole_fractions[vapor->nDryMixtures() + i - 1]/xdry;
    Real beta = vapor->getLatentHeat_JouleKilomole(i)/(Cantera::GasConstant*temperature);

    c1 += eta*beta;
    c2 += eta*beta*beta;
    c3 += eta;

    /*if (isat[iv] > 0) {
      Real latent = thermo[iv].beta*thermo[iv].t3/q[idn] - thermo[iv].delta;
      c1 += q[iv]/xd*latent;
      c2 += q[iv]/xd*latent*latent;
    }
    c3 += q[iv]/xd;*/
  }

  return (1. + c1)/(cphat_ov_r + (c2 + c1*c1)/(1. + c3));
}

#endif
