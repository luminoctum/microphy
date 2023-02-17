#ifndef VAPOR_HPP
#define VAPOR_HPP

#include <configure.hpp>
#include <cliutils/vectorize.hpp>
#include <parameter_input.hpp>
#include <cantera/base/stringUtils.h>
#include <cantera/thermo.h>

template<typename VaporPhase>
class Vapor : public VaporPhase
{
public:
  //Vapor(MeshBlock *pmb, ParamterInput *pin);
  Vapor(ParameterInput *pin):
    VaporPhase(pin->GetString("chemistry", "kinetics_file"), "vapor")
  {
    std::string str = pin->GetString("hydro", "vapors");
    std::vector<std::string> vapor_comps = Vectorize<std::string>(str.c_str(), ",");
    for (size_t i = 0; i < vapor_comps.size(); ++i)
      this->setMoleFractionsByName(vapor_comps[i]);
    m_num_dry_mixtures = Cantera::parseCompString(vapor_comps[0]).size();
    m_molwt_ratios[0] = 1.;
    m_cv_mass_ratios[0] = 1.;
  }

  ~Vapor() {}

  size_t nDryMixtures() const {
    return m_num_dry_mixtures;
  }

  Real getMolecularWeightRatio(int n) const {
    return m_molwt_ratios[n];
  }

  Real getCvMassRatio(int n) const {
    return m_cv_mass_ratios[n];
  }

  Real getGamma() const {
    return this->cp_mole()/this->cv_mole();
  }

  Real getGammaDry() const {
    Real cp = 0., cv = 0.;
    for (int i = 0; i < nDryMixtures(); ++i) {
      cp += this->m_cp0_R[i]*this->moleFraction(i);
      cv += (this->m_cp0_R[i] - 1.)*this->moleFraction(i);
    }
    return cp/cv;
  }

  // only works for ideal gas
  template<typename T>
  Real getTemp(T w) const {
    return w[ipr]*this->meanMolecularWeight()/(w[idn]*Cantera::GasConstant);
  }

  //! Update the species reference state thermodynamic functions
  /*!
   *  This method is called each time a thermodynamic property is requested,
   *  to check whether the internal species properties within the object
   *  need to be updated. Currently, this updates the species thermo
   *  polynomial values for the current value of the temperature. A check is
   *  made to see if the temperature has changed since the last evaluation.
   *  This object does not contain any persistent data that depends on the
   *  concentration, that needs to be updated. The state object modifies its
   *  concentration dependent information at the time the setMoleFractions()
   *  (or equivalent) call is made.
   */
  void updateThermodynamics(Real temperature, Real pressure) {
    this->setTemperature(temperature);
    this->setPressure(pressure);
    VaporPhase::updateThermo();

    Real cp_dry = 0., ym = 0.;
    Real m_mmw_dry = 0.;

    for (int i = 0; i < nDryMixtures(); ++i) {
      cp_dry += this->m_cp0_R[i]*this->moleFraction(i);
      m_mmw_dry += this->molecularWeight(i)*this->moleFraction(i);
      ym += this->moleFraction(i);
    }

    cp_dry /= ym;
    m_mmw_dry /= ym;

    for (int i = 1; i <= NumVapors; ++i)
      m_cv_mass_ratios[i] = (this->m_cp0_R[i] - 1.)*m_mmw_dry/(cp_dry - 1.)
        /this->molecularWeight(i);
  }

protected:
  size_t                        m_num_dry_mixtures;
  Real                          m_mmw_dry;
  std::array<Real,1+NumVapors>  m_molwt_ratios;
  std::array<Real,1+NumVapors>  m_cv_mass_ratios;
};

#endif
