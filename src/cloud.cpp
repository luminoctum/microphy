#include <iostream>
#include <mesh/mesh.hpp>
#include <parameter_input.hpp>
#include <scalars/scalars.hpp>
#include <debugger.hpp>
#include <molecule/molecules.hpp>
#include <cliutils/vectorize.hpp>
#include "cloud.hpp"

namespace Globals {
  extern int next_scalar_index;
}

std::ostream& operator<<(std::ostream& os, Cloud const& my)
{
  for (int i = 0; i < my.molecules_.size(); ++i)
    std::cout << "#" << i << " " << my.molecules_[i] << std::endl;
  return os;
}

Cloud::Cloud(MeshBlock *pmb, ParameterInput *pin)
{
  if (!pin->DoesParameterExist("hydro", "clouds")) {
    pdebug->Leave();
    return;
  }

  std::string str = pin->GetString("hydro", "clouds");
  std::vector<std::string> names = Vectorize<std::string>(str.c_str(), ", ");

  int ns = Globals::next_scalar_index;
  // den, tem
  Globals::next_scalar_index += 2*names.size();

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
      molecules_[i] = Water<1>();
    } else if (name == "H2O(s)") {
      molecules_[i] = Water<2>();
    } else if (name == "NH3(l)") {
      molecules_[i] = Ammonia<1>();
    } else if (name == "NH3(s)") {
      molecules_[i] = Ammonia<2>();
    } else if (name == "H2S(l)") {
      molecules_[i] = HydrogenSulfide<1>();
    } else if (name == "H2S(s)") {
      molecules_[i] = HydrogenSulfide<2>();
    } else if (name == "NH4SH(s)") {
      molecules_[i] = AmmoniumHydrosulfide<2>();
    } else {
      Debugger::Fatal("GasPhase", name, "unrecognized");
    }

    if (p != NULL) {
      molecules_[i].id = p+1;
    }
  }

  q.InitWithShallowSlice(pscalars_->r, 4, ns, names.size());
  c.InitWithShallowSlice(pscalars_->s, 4, ns, names.size());
  tem.InitWithShallowSlice(pscalars_->r, 4, ns + names.size(), names.size());
  energy.InitWithShallowSlice(pscalars_->s, 4, ns + names.size(), names.size());

  pdebug->Leave();
}
