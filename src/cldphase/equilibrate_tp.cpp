#include <debugger.hpp>
#include "cloud_phase.hpp"

void CloudPhase::EquilibrateTP(Real qmol[], Real temp, Real pres) const
{
  qmol[idn] = temp;
  qmol[ipr] = pres;

  Real rate_max;
  Real ftol = 1.E-4;
  int max_iter = 10;
  int iter = 0;
  do {
    rate_max = 0;
    /*for (int i = 0; i < reactions_.size(); ++i) {
      Real rate = reactions_[i].fptr(qmol, reactions_[i].indx, pgas_, this);
      for (int j = 0; j < MaxReactants; ++j) {
        int k = reactions_[i].indx[j];
        if (k == 0) break;
        if (k < 0) {
          qmol[-k] -= rate;
        } else {
          qmol[k] += rate;
        }
      }
      rate_max = std::max(std::abs(rate), rate_max);
    }*/
  } while ((rate_max > ftol) && (iter++ < max_iter));

  if (iter >= max_iter) {
    Debugger::Fatal("CloudPhase::EquilibrateTP", "fails to converge");
  }
}
