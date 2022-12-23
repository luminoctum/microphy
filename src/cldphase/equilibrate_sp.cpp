#include <sstream>
#include <debugger.hpp>

extern "C" {
  #include <root.h>
}
#include "cloud_phase.hpp"

double CloudPhase::solve_entropy_(double T, void *aux) {
  EntropySolver *psolver = static_cast<EntropySolver*>(aux);

  psolver->pcld->EquilibrateTP(psolver->qmol, T, psolver->pres);
  return psolver->pcld->getMolarEntropy(psolver->qmol) - psolver->entropy;
}

void CloudPhase::EquilibrateSP(Real qmol[], Real entropy, Real pres) const
{
  Real dt = 1.5, Tmin = qmol[idn], Tmax = qmol[idn];

  EquilibrateTP(qmol, Tmin, pres);
  while (getMolarEntropy(qmol) > entropy) {
    Tmin /= dt;
    EquilibrateTP(qmol, Tmin, pres);
  }

  EquilibrateTP(qmol, Tmax, pres);
  while (getMolarEntropy(qmol) < entropy) {
    Tmax *= dt;
    EquilibrateTP(qmol, Tmax, pres);
  }

  EntropySolver  solver;

  solver.qmol = qmol;
  solver.pres = pres;
  solver.entropy = entropy;
  solver.pcld = this;

  Real Tc;
  Real ftol = 1.E-4;
  int error = root(Tmin, Tmax, ftol, &Tc, solve_entropy_, &solver);

  if (error) {
    std::stringstream msg;
    msg << "Tmin = " << Tmin << std::endl;
    msg << "Tmax = " << Tmax << std::endl;
    msg << "Required Entropy = " << entropy << std::endl;
    EquilibrateTP(qmol, Tmin, pres);
    msg << "Entropy at Tmin = " << getMolarEntropy(qmol) << std::endl;
    EquilibrateTP(qmol, Tmax, pres);
    msg << "Entropy at Tmax = " << getMolarEntropy(qmol) << std::endl;
    Debugger::Fatal("Clouds::EquilibrateSP", "Fails to converge", msg.str());
  }

  // take care of mixed phase
  Real *qmol1 = new Real [NumHydros + clouds_.size()];
  Real *qmol2 = new Real [NumHydros + clouds_.size()];

  for (int i = 0; i < NumHydros + clouds_.size(); ++i) {
    qmol1[i] = qmol[i];
    qmol2[i] = qmol[i];
  }

  EquilibrateTP(qmol1, Tc - 1.E-4, pres);
  Real s1 = getMolarEntropy(qmol1);

  EquilibrateTP(qmol2, Tc + 1.E-4, pres);
  Real s2 = getMolarEntropy(qmol2);

  Real xfrac = (s2 - entropy)/(s2 - s1);

  for (int i = 1; i < NumHydros + clouds_.size(); ++i)
    qmol[i] = xfrac*qmol1[i] + (1. - xfrac)*qmol2[i];

  qmol[idn] = Tc;
  qmol[ipr] = pres;
}
