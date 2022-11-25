#ifndef CLOUD_PHASE_HPP
#define CLOUD_PHASE_HPP

#include <iosfwd>
#include <vector>
#include <athena.hpp>
#include <molecule/molecules.hpp>
#include <configure.hpp>
#include "../reaction/cloud_reaction.hpp"

class MeshBlock;
class ParameterInput;
class GasPhase;
class CloudPhase;
class PassiveScalars;

struct EntropySolver {
  Real                  *qmol;
  Real                  pres;
  Real                  entropy;
  CloudPhase const      *pcld;
};

class CloudPhase {
  friend std::ostream& operator<<(std::ostream& os, CloudPhase const& my);
public:
  CloudPhase(MeshBlock *pmb, ParameterInput *pin);

  // conserved variables (extensive), kg/m^3
  AthenaArray<Real> c;

  // primitive variables (intensive), kg/kg
  AthenaArray<Real> q;

  void EquilibrateTP(Real qmol[], Real temp, Real pres) const;

  void EquilibrateSP(Real qmol[], Real entropy, Real pres) const;

  template<typename T> Real getMolarEntropy(T qmol) const;

protected:
  PassiveScalars        *pscalars_;

  GasPhase const        *pgas_;

  bool                  active_;

  int                   ns_, ne_;

  std::vector<Molecule> clouds_;

  std::vector<CloudReaction>
                        reactions_;

  int                   max_iter_;

  Real                  ftol_;

  static double solve_entropy_(double, void*);
};


#endif
