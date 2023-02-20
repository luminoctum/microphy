#ifndef SUSPENSION_HPP
#define SUSPENSION_HPP

#include <cantera/thermo.h>

template<typename VaporPhase, typename SolidPhase>
class Suspension : public SolidPhase
{
public:
  Suspension(ParameterInput *pin, VaporPhase const *vapor, std::string name):
    SolidPhase(pin->GetString("chemistry", "kinetics_file"), name),
    m_vapor(vapor), m_mole_ratio(0), m_mass_ratio(0)
  {}

  void setRelativeMole(Real mf) {
    m_mole_ratio = mf;
    m_mass_ratio = this->meanMolecularWeight()*mf/
      (m_vapor->meanMolecularWeight() + this->meanMolecularWeight()*mf);
  }

  void setRelativeMass(Real mf) {
    m_mass_ratio = mf;
    m_mole_ratio = mf/this->meanMolecularWeight()/
      (1./m_vapor->meanMolecularWeight() + mf/this->meanMolecularWeight());
  }

  Real getRelativeMole() const {
    return m_mole_ratio;
  }

  Real getRelativeMass() const {
    return m_mass_ratio;
  }

  void getConcentrations(double* const c) const override {
    Real mols = m_vapor->molarDensity()*m_mole_ratio;
    for (size_t i = 0; i < this->nSpecies(); ++i) {
      c[i] = mols*this->moleFraction(i);
    }
  }

  void getActivityConcentrations(doublereal* c) const override {
    SolidPhase::getActivityConcentrations(c);
    for (size_t i = 0; i < this->nSpecies(); ++i) {
      c[i] *= m_vapor->density()/this->density();
    }
  }

  double concentration(const size_t k) const override {
    return m_vapor->molarDensity()*m_mole_ratio*this->moleFraction(k);
  }

protected:
  VaporPhase const *m_vapor;

  //! molar fraction of total solution
  Real m_mole_ratio;

  //! mass fraction of total solution
  Real m_mass_ratio;
};

#endif
