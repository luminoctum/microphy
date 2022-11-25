#include <iostream>
#include <mesh/mesh.hpp>
#include <parameter_input.hpp>
#include <scalars/scalars.hpp>
#include <debugger.hpp>
#include <molecule/molecules.hpp>
#include "cloud_phase.hpp"
#include "../vectorize.hpp"

namespace Globals {
  extern int next_scalar_index;
}

std::ostream& operator<<(std::ostream& os, CloudPhase const& my)
{
  for (int i = 0; i < my.clouds_.size(); ++i)
    std::cout << "#" << i << " " << my.clouds_[i] << std::endl;
  return os;
}

CloudPhase::CloudPhase(MeshBlock *pmb, ParameterInput *pin):
  ns_(Globals::next_scalar_index), active_(false)
{
  pdebug->Enter("CloudPhase");
  if (!pin->DoesParameterExist("hydro", "clouds")) {
    pdebug->Leave();
    return;
  }

  active_ = true;

  std::string str = pin->GetString("hydro", "clouds");
  std::vector<std::string> names = Vectorize<std::string>(str.c_str(), ", ");

  Globals::next_scalar_index += names.size();
  ne_ = Globals::next_scalar_index - 1;

  if (ne_ >= NumScalars) {
    Debugger::Fatal("CloudPhase", "Number of scalars less than",
      std::to_string(ne_));
  }

  // clouds
  for (int i = 0; i < names.size(); ++i) {
    char buf[80];
    strcpy(buf, names[i].c_str());

    char *p = strchr(buf, ':');
    if (p != NULL) {
      *p = '\0';
    }

    std::string name = buf;
    if (name == "H2O(l)") {
      clouds_[i] = Water<1>();
    } else if (name == "H2O(s)") {
      clouds_[i] = Water<2>();
    } else if (name == "NH3(l)") {
      clouds_[i] = Ammonia<1>();
    } else if (name == "NH3(s)") {
      clouds_[i] = Ammonia<2>();
    } else if (name == "H2S(l)") {
      clouds_[i] = HydrogenSulfide<1>();
    } else if (name == "H2S(s)") {
      clouds_[i] = HydrogenSulfide<2>();
    } else if (name == "NH4SH(s)") {
      clouds_[i] = AmmoniumHydrosulfide<2>();
    } else {
      Debugger::Fatal("GasPhase", name, "unrecognized");
    }

    if (p != NULL) {
      clouds_[i].id = p+1;
    }
  }

  c.InitWithShallowSlice(pscalars_->s, 4, ns_, ne_ - ns_ + 1);
  q.InitWithShallowSlice(pscalars_->r, 4, ns_, ne_ - ns_ + 1);

  pdebug->Leave();
}
