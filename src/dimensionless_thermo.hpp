#ifndef DIMENSIONLESS_THERMO_HPP
#define DIMENSIONLESS_THERMO_HPP

// C/C++
#include <iosfwd>

// canoe
#include <configure.hpp>

struct DimensionlessThermo {
  friend std::ostream& operator<<(std::ostream& os, DimensionlessThermo const& my);
  DimensionlessThermo() :
    mu_ratio(1), cp_ratio(1), cv_ratio(1),
    t3(0), p3(0), beta(0)
  {}

  //! ratio of mean molecular weights
  Real                        mu_ratio;

  //! ratio of specific heat capacities at constant pressure
  Real                        cp_ratio;

  //! ratio of specific heat capacities at constant volume
  Real                        cv_ratio;

  //! triple point temperature [K]
  Real                        t3;

  //! triple point pressure [pa]
  Real                        p3;

  /*! dimensionless latent heat
   *$\beta_{ij} == \frac{L_{ij}^r + \Delta c_{ij}T^r}{R_i T^r} $
   */
  Real                        beta;

  /*! dimensionless differences in specific heat capacity
   *$(c_{ij} - c_{p,i})/R_i$
   */
  Real                        delta;
};

#endif
