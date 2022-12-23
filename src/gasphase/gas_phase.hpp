/** @file gas_phase.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Tuesday May 25, 2021 18:04:45 UTC
 * @bug No known bugs.
 */

#ifndef GAS_PHASE_HPP
#define GAS_PHASE_HPP

// Athena++
#include <eos/eos.hpp>

// C/C++ headers
#include <cfloat>
#include <iosfwd>
#include <map>
#include <vector>

#include <molecule/gas_mixture.hpp>
#include <configure.hpp>
#include "gas_thermo.hpp"

class Molecule;
class MeshBlock;
class ParameterInput;

// dynamic variables are ordered in the array as the following
// 0: dry air (non-condensible)
// iphase = 0 - moist air (condensible)
// iphase = 1 - primary condensible species
// iphase = 2..(N-1) - all other condensible species

enum class Adiabat {reversible = 0, pseudo = 1, dry = 2};

class GasPhase {
  friend std::ostream& operator<<(std::ostream& os, GasPhase const& my);
public:
  // member functions
  GasPhase(MeshBlock *pmb, ParameterInput *pin);
  ~GasPhase() {}

  GasMixture const& getDry() const {
    return dry_;
  }

  Molecule const& getVapor(int i) const {
    return moist_[i];
  }

  /*! Ideal gas constant of dry air in [J/(kg K)]
   * @return $R_d=\hat{R}/m_d$
   */
  Real getRd() const {
    return Rd_;
  }

  /*! Ratio of specific heat capacity at constant volume
   * @param i index of the vapor
   * @return $c_{v,i}/c_{v,d}$
   */
  Real getCvRatio(int i) const {
    return thermo_[i].cv_ratio;
  }

  /*! Specific heat capacity at constant volume
   *
   * $c_{v,d} = \frac{R_d}{\gamma_d - 1}$ \n
   * $c_{v,i} = \frac{c_{v,i}}{c_{v,d}}\times c_{v,d}$
   * @param i index of the vapor
   * @return $c_v$ [J/(kg K)]
   */
  Real getCv(int i) const {
    Real cvd = Rd_/(peos_->getGamma() - 1.);
    return getCvRatio(i)*cvd;
  }

  /*! Ratio of specific heat capacity at constant pressure
   * @param i index of the vapor
   * @return $c_{p,i}/c_{p,d}$
   */
  Real getCpRatio(int i) const {
    return thermo_[i].cp_ratio;
  }

  /*! Specific heat capacity at constant pressure
   *
   * $c_{p,d} = \frac{\gamma_d}{\gamma_d - 1}R_d$ \n
   * $c_{p,i} = \frac{c_{p,i}}{c_{p,d}}\times c_{p,d}$
   * @param i index of the vapor
   * @return $c_p$ [J/(kg K)]
   * */
  Real getCp(int i) const {
    Real gamma = peos_->getGamma();
    Real cpd = Rd_*gamma /(gamma - 1.);
    return getCpRatio(i)*cpd;
  }

  /*/ Temperature dependent specific latent heat of condensates at constant volume
   *
   * $L_{ij}(T) = L_{ij}^r - (c_{ij} - c_{p,i})\times(T - T^r)$
   * $= L_{ij}^r - \delta_{ij}R_i(T - T^r)$
   * @param j index of condensate
   * @param temp temperature (default to 0)
   * @return $L_{ij}(T)$ [J/kg]
  Real getLatent(int j, Real temp = 0.) const {
    return thermo_[j].latent - 
      thermo_[j].delta*Rd_/thermo_[j].mu_ratio*temp;
  }*/

  /*! Ratio of molecular weights
   * @param i index of the vapor or the condensate
   * @return $\epsilon_i=m_i/m_d$
   */
  Real getMassRatio(int i) const {
    return thermo_[i].mu_ratio;
  }

  /*! $\beta$ parameter of a condensate
   *
   * $\beta_{ij} = \frac{\Delta U_{ij}}{R_i T^r}$ \n
   * $\Delta U_{ij}$ is the difference in internal energy between the vapor $i$ 
   * and the condensate $j$ \n
   * $R_i=\hat{R}/m_i=R_d/\epsilon_i$ \n
   * $T^r$ is the triple point temperature
   * @param j index of the condensate
   * @return $\beta_{ij}$
   */
  Real getBeta(int j) const {
    return thermo_[j].beta;
  }

  /*! $\delta$ parameter of a condensate
   *
   * $\delta_{ij} = \frac{\Delta c_{ij}}{R_i}$ \n
   * $\Delta c_{ij} = c_{ij} - c_{p,j}$ is the difference in specific heat capacity 
   * at constant temperature.\n
   * $c_{ij}$ is the heat capacity of condensate $j$ of vapor $i$
   * @param j index of the condensate
   * @return $\delta_{ij}$
   */
  Real getDelta(int j) const {
    return thermo_[j].delta;
  }

  // conversion functions
  //! Change mass mixing ratio to molar mixing ratio
  template<typename T1, typename T2>
  void PrimitiveToChemical(T1 qmol, T2 const w) const;

  //! Change molar mixing ratio to mass mixing ratio
  template<typename T1, typename T2>
  void ChemicalToPrimitive(T1 w, T2 const qmol) const;

  //! Change density to molar mixing ratio
  template<typename T1, typename T2>
  void ConservedToChemical(T1 qmol, T2 const u) const;

  //! Change molar mixing ratio to density
  template<typename T1, typename T2>
  void ChemicalToConserved(T1 u, T2 const qmol) const;

  // access functions
  //! polytropic index $\gamma=c_p/c_v$
  template<typename T> Real getGamma(T w) const;

  template<typename T> Real RovRd(T w) const;

  //! Temperature
  template<typename T> Real getTemp(T w) const;

  template<typename T> Real getPres(T u) const;

  template<typename T> Real getChi(T w) const;

  //! c_p = c_{pd}*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
  template<typename T> Real getSpecificCp(T w) const;

  //! c_v = c_{vd}*(1 + \sum_i (q_i*(\hat{c}_{vi} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{vi} - 1.)))
  template<typename T> Real getSpecificCv(T w) const;

  //! h = c_{pd}*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
  //!   = \gamma_d/(\gamma_d - 1.)*R_d*T*(1 + \sum_i (q_i*(\hat{c}_{pi} - 1.)))
  template<typename T> Real getSpecificEnthalpy(T w) const;

  template<typename T> Real getMeanMolecularWeight(T w) const;

  /*! Construct an 1d atmosphere
   * @param method choose from [reversible, pseudo, dry, isothermal]
   */
  void ConstructAtmosphere(Real **w, Real Ts, Real Ps,
    Real grav, Real dzORdlnp, int len, Adiabat method, Real userp) const;

protected:
  EquationOfState const       *peos_;

  GasMixture                  dry_;
  Molecule                    moist_[1+NumVapors];

  //! ideal gas constant of dry air in J/kg
  Real                        Rd_;

  //! reference temperature
  Real                        Tref_;

  //! thermodynamic data
  GasThermo                   thermo_[1+NumVapors];

  Real ftol_;
};

void update_gamma(Real& gamma, Real const qmol[]);

#include "gasphase_conversion.hpp"
#include "gasphase_access.hpp"

#endif
