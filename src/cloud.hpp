#ifndef CLOUD_HPP
#define CLOUD_HPP

#include <map>
#include <iosfwd>
#include <vector>
#include <athena.hpp>
#include <molecule/molecules.hpp>
#include <configure.hpp>

class MeshBlock;
class ParameterInput;
class Cloud;
class PassiveScalars;
class Microphysics;

struct EntropySolver {
  Real             *qmol;
  Real             pres;
  Real             entropy;
  Cloud const      *pcld;
};

class Cloud {
  friend std::ostream& operator<<(std::ostream& os, Cloud const& my);
public:
// public access variables

  // conserved variables (extensive), kg/m^3
  AthenaArray<Real> c;

  // primitive variables (intensive), kg/kg
  AthenaArray<Real> q;

  // temperature and energy
  AthenaArray<Real> tem, energy;

  // differential velocity 1, 2, 3
  AthenaArray<Real> delv1, delv2, delv3;

  // flux in 3-dimensions
  AthenaArray<Real>     flux[3];

  Cloud(MeshBlock *pmb, ParameterInput *pin);

  Molecule const& getMolecule(int i) const {
    return molecules_[i];
  }

  void calcSedimentaionFlux();

  void EquilibrateTP(Real qmol[], Real temp, Real pres) const;

  void EquilibrateSP(Real qmol[], Real entropy, Real pres) const;

  void EquilibrateUV(Real qmol[], Real energy, Real volume) const;

  template<typename T> Real getMolarEntropy(T qmol) const;

protected:
  PassiveScalars        *pscalars_;

  std::vector<Molecule> molecules_;

  static double solve_entropy_(double, void*);
};

#endif
