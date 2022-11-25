// C/C++ header
#include <vector>  // fill
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <stdexcept>

// Athena++ 
#include <parameter_input.hpp>
#include <mesh/mesh.hpp>

#include <debugger.hpp>
#include <molecule/molecules.hpp>
#include "gasphase.hpp"
#include "../vectorize.hpp"

std::ostream& operator<<(std::ostream& os, GasPhase const& my)
{
  std::cout << "Dry component = " << std::endl;
  std::cout << my.dry_ << std::endl;
  std::cout << "Rd = " << my.Rd_ << std::endl;
  std::cout << "Tref = " << my.Tref_ << std::endl;
  std::cout << "GasThermo data = "<< std::endl;
  for (int i = 0; i <= NumVapors; ++i)
    std::cout << "#" << i << " " << my.thermo_[i] << std::endl;
  return os;
}

GasPhase::GasPhase(MeshBlock *pmb, ParameterInput *pin):
  peos_(pmb->peos), Tref_(pin->GetOrAddReal("hydro", "Tref", 300.))
{
  pdebug->Enter("GasPhase");
  std::string str = pin->GetOrAddString("hydro", "vapors", "Dummy:1.");
  std::vector<std::string> names = Vectorize<std::string>(str.c_str(), ", ");

  if (names.size() <= NumVapors) {
    Debugger::Fatal("GasPhase", "not enough vapors");
  }

  // dry component
  dry_ = GasMixture(names[0]);
  moist_[0] = dry_;

  // temporary clouds
  Molecule clouds[1+NumVapors];

  // vapors
  for (int i = 1; i <= NumVapors; ++i) {
    if (names[i] == "H2O") {
      moist_[i] = Water<0>();
      clouds[i] = Water<ReferenceCloudPhase>();
    } else if (names[i] == "NH3") {
      moist_[i] = Ammonia<0>();
      clouds[i] = Ammonia<ReferenceCloudPhase>();
    } else if (names[i] == "CH4") {
      moist_[i] = Methane<0>();
      clouds[i] = Methane<ReferenceCloudPhase>();
    } else if (names[i] == "H2S") {
      moist_[i] = HydrogenSulfide<0>();
      clouds[i] = HydrogenSulfide<ReferenceCloudPhase>();
    } else {
      Debugger::Fatal("GasPhase", names[i], "unrecognized");
    }
  }

  // set up thermodynamic data
  Real gamma = peos_->getGamma();
  for (int i = 1; i <= NumVapors; ++i) {
    thermo_[i].mu_ratio = moist_[i].mu/dry_.mu; 
    thermo_[i].cp_ratio = moist_[i].cp(Tref_)/dry_.cp(Tref_)/thermo_[i].mu_ratio;
    // cv_ratio = $\gamma\sigma_i + (1. - \gamma)/\epsilon_i$
    thermo_[i].cv_ratio = gamma*thermo_[i].cp_ratio 
      + (1. - gamma)/thermo_[i].mu_ratio;
    thermo_[i].t3 = moist_[i].tripleT; 
    thermo_[i].p3 = moist_[i].tripleP; 
    thermo_[i].beta = -clouds[i].getEnthalpyOffset()/(Constants::Rgas*moist_[i].tripleT);
    // delta = $(\sigma_j - \sigma_i)*\epsilon_i*\gamma/(\gamma - 1)$
    Real cloud_cp_ratio = clouds[i].cp(Tref_)/dry_.cp(Tref_)/thermo_[i].mu_ratio;
    thermo_[i].delta = (cloud_cp_ratio - thermo_[i].cp_ratio)*
      thermo_[i].mu_ratio/(1. - 1./gamma);
  }

  if (pin->DoesParameterExist("hydro", "Rd")) {
    Rd_ = pin->GetReal("hydro", "Rd");
  } else {
    Rd_ = Constants::Rgas/dry_.mu;
    pin->GetOrAddReal("hydro", "Rd", Rd_);
  }

  ftol_ = pin->GetOrAddReal("hydro", "ftol", 1.0E-4);

  pdebug->Leave();
}

void __attribute__((weak)) update_gamma(Real& gamma, Real const q[]) {}
