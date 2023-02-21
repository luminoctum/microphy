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

    std::fill(m_latent_heats.begin(), m_latent_heats.end(), 0.);
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
  }

  Real getGamma() const {
    return this->cp_mole()/this->cv_mole();
  }

  Real dry_cp_mass() const {
    return m_dry_cp_mole*Cantera::GasConstant/m_dry_mmw;
  }

  Real dry_cv_mass() const {
    return (m_dry_cp_mole - 1.)*Cantera::GasConstant/m_dry_mmw;
  }*/

  Real getDryGamma() const {
    return m_dry_cp_mole/(m_dry_cp_mole - 1.);
  }

  size_t vaporIndex(std::string const& name) const {
    return this->speciesIndex(name) - nDryMixtures() + 1;
  }

  std::string vaporName(size_t k) const {
    return this->speciesName(k + nDryMixtures() - 1);
  }

  bool testSaturation(size_t k) const {
    return m_latent_heats[k] > 0;
  }

  Real getLatentHeat_JouleKilomole(size_t k) const {
    return m_latent_heats[k];
  }

  // to test whether a vapor is saturated
  // we test saturation by adding small amount of condensates
  // this function sets the latent heat variable
  // if the vapor is not saturated, latent heat is 0. Otherwise, it is the enthalpy
  // difference between vapor and cloud
  void setLatentHeats(Cantera::Kinetics *kin) {
    size_t nphases = kin->nPhases();

    if (m_wdot.size() < kin->nTotalSpecies()) {
      m_wdot.resize(kin->nTotalSpecies());
    }

    if (m_deltaH.size() < kin->nReactions()) {
      m_deltaH.resize(kin->nReactions());
    }

    std::vector<Real> hbar(this->nSpecies());
    this->getPartialMolarEnthalpies(hbar.data());
    std::cout << "enthalpies = ";
    for (int i = 0; i < this->nSpecies(); ++i)
      std::cout << hbar[i] << " ";
    std::cout << std::endl;

    kin->getNetProductionRates(m_wdot.data());
    kin->getDeltaEnthalpy(m_deltaH.data());
    Real deltaH = 0;
    for (auto it = m_deltaH.begin(); it != m_deltaH.end(); ++it)
      deltaH += *it;
    for (int i = 1; i <= NumVapors; ++i) {
      std::string name = this->speciesName(nDryMixtures() + i - 1);
      size_t j = kin->kineticsSpeciesIndex(name);
      if (m_wdot[j] > 0.) {
        m_latent_heats[i] = 0.;
        //m_isat[i] = false;
      } else {
        m_latent_heats[i] = -deltaH;
        std::cout << ">>> I'm setting latent heat = " << -deltaH << " <<<" << std::endl;
        //m_isat[i] = true;
      }
    }
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
  void updateThermo() {
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

  //! the following functions are for interaction with hydrodynamics
  template<typename T> Real getGamma(T w) const;

  template<typename T> Real getMeanMolecularWeight(T w) const;

  //! only works for ideal gas
  template<typename T> Real getTemp(T w) const;

  template<typename T> Real getPres(T u) const;

  template<typename T> Real getRovCp(T w) const;

  template<typename T> Real cp_JouleKilogramK(T w) const;

  template<typename T> Real cv_JouleKilogramK(T w) const;

protected:
  //! mean molecular weight of dry gas
  Real                          m_dry_mmw;

  //! number of dry species whose mixing ratios never change
  size_t                        m_num_dry_mixtures;

  //! molecular weight ratios between vapors and dry gas
  std::array<Real,1+NumVapors>  m_molwt_ratios;

  //! dimensionless molar heat capacity at constant pressure of dry gas
  Real                          m_dry_cp_mole;

  //! mass heat capacity ratios at constant volume
  std::array<Real,1+NumVapors>  m_cv_ratios_mass;

  //! mass heat capacity ratios at constant pressure
  std::array<Real,1+NumVapors>  m_cp_ratios_mass;

  //! scract arrays, production rate
  std::vector<Real>             m_wdot, m_deltaH;

  //! m_latent_heats stores the latent heat of a vapor if it condenses
  std::array<Real,1+NumVapors>  m_latent_heats;
};

#include "vapor_hydro.hpp"

#endif
