#ifndef VAPOR_HPP
#define VAPOR_HPP

#include <sstream>
#include <configure.hpp>
#include <debugger.hpp>
#include <cliutils/vectorize.hpp>
#include <parameter_input.hpp>
#include <cantera/base/stringUtils.h>
#include <cantera/thermo.h>
#include <cantera/kinetics.h>

template<typename BasePhase>
class Vapor : public BasePhase
{
public:
  //Vapor(MeshBlock *pmb, ParamterInput *pin);
  Vapor(ParameterInput *pin):
    BasePhase(pin->GetString("chemistry", "kinetics_file"), "vapor"),
    m_dry_mmw(0.), m_dry_cp_mole(0.)
  {
    std::stringstream msg;
    std::string compstr = pin->GetString("hydro", "vapors");
    std::vector<std::string> dry_comps = Vectorize<std::string>(
      compstr.substr(0, compstr.find(',')).c_str(), ";");
    std::replace(compstr.begin(), compstr.end(), ',', ';');

    this->setMoleFractionsByName(compstr);
    m_num_dry_mixtures = dry_comps.size();

    if (this->nSpecies() - nDryMixtures() > NumVapors) {
      msg << this->nSpecies() - nDryMixtures() << " vapors needed, but "
          << NumVapors << " configured. Increase the number of vapors";
      Debugger::Fatal("Vapor", msg.str());
    }

    Real ym = 0.;
    for (int i = 0; i < nDryMixtures(); ++i) {
      m_dry_mmw += this->molecularWeight(i)*this->moleFraction(i);
      m_dry_cp_mole += this->m_cp0_R[i]*this->moleFraction(i);
      ym += this->moleFraction(i);
    }
    m_dry_mmw /= ym;
    m_dry_cp_mole /= ym;

    m_molwt_ratios[0] = 1.;
    m_cv_ratios_mass[0] = 1.;
    for (int i = 1; i <= NumVapors; ++i) {
      m_molwt_ratios[i] = this->molecularWeight(nDryMixtures() + i - 1)/m_dry_mmw;
    }

    updateThermo();
    /*std::cout << this->m_cp0_R[0] << std::endl;
    std::cout << m_dry_mmw << std::endl;
    std::cout << m_dry_cp_mole << std::endl;

    for (int i = 0; i <= NumVapors; ++i)
      std::cout << m_molwt_ratios[i] << " " << m_cv_ratios_mass[i] << std::endl;
    */
  }

  ~Vapor() {}

  size_t nDryMixtures() const {
    return m_num_dry_mixtures;
  }

  /*Real getMolecularWeightRatio(int n) const {
    return m_molwt_ratios[n];
  }

  Real getCvRatio_mass(int n) const {
    return m_cv_ratios_mass[n];
  }*/

  Real getGamma() const {
    return this->cp_mole()/this->cv_mole();
  }

  Real dry_cp_mass() const {
    return m_dry_cp_mole*Cantera::GasConstant/m_dry_mmw;
  }

  Real dry_cv_mass() const {
    return (m_dry_cp_mole - 1.)*Cantera::GasConstant/m_dry_mmw;
  }

  Real getDryGamma() const {
    return m_dry_cp_mole/(m_dry_cp_mole - 1.);
  }

  // to test whether a vapor is saturated
  // we test saturation by adding small amount of condensates
  void setSaturationIndicator(Cantera::Kinetics &kin) {
    size_t nphases = kin.nPhases();

    if (m_wdot.size() < kin.nTotalSpecies()) {
      m_wdot.resize(kin.nTotalSpecies());
    }

    kin.getNetProductionRates(m_wdot.data());
    for (int i = 1; i <= NumVapors; ++i) {
      std::string name = this->speciesName(i);
      size_t j = kin.kineticsSpeciesIndex(name);
      if (m_wdot[j] > 0.) {
        m_isat[i] = false;
      } else {
        m_isat[i] = true;
      }
    }
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
  void updateThermo() const override {
    BasePhase::updateThermo();

    Real ym = 0.;
    m_dry_cp_mole = 0.;

    for (int i = 0; i < nDryMixtures(); ++i) {
      m_dry_cp_mole += this->m_cp0_R[i]*this->moleFraction(i);
      ym += this->moleFraction(i);
    }

    m_dry_cp_mole /= ym;

    for (int i = 1; i <= NumVapors; ++i) {
      m_cv_ratios_mass[i] = (this->m_cp0_R[nDryMixtures() + i - 1] - 1.)
        /(m_dry_cp_mole - 1.);
      m_cv_ratios_mass[i] /= m_molwt_ratios[i];
    }
  }

protected:
  //! mean molecular weight of dry gas
  Real                          m_dry_mmw;

  //! dimensionless molar heat capacity at constant pressure of dry gas
  mutable Real                  m_dry_cp_mole;

  //! number of dry species whose mixing ratios never change
  size_t                        m_num_dry_mixtures;

  //! molecular weight ratios between vapors and dry gas
  std::array<Real,1+NumVapors>  m_molwt_ratios;

  //! mass heat capacity ratios at constant volume
  mutable std::array<Real,1+NumVapors>  m_cv_ratios_mass;

  //! scract arrays, net Production rate
  std::vector<Real>             m_wdot;

  //! latent heat of condensation in J/(kmol)
  Real                          m_molar_latent_heat;

  //! m_isat is an indicator array of whether vapor i is saturated or not
  std::array<bool,1+NumVapors>  m_isat;
};

#endif
