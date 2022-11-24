#ifndef CONDENSE_VAPOR_HPP
#define CONDENSE_VAPOR_HPP

#include <cmath>

#include <molecule/molecule.hpp>
#include <configure.hpp>

//! solving for
//! \frac{x_v - r}{1 - r} = s
//! r = x_v - \frac{s}{1-s}(1 - x_v)
inline Real condense_vapor(Real const xmol[], int iv, Molecule const& moist)
{
  Real xv = xmol[iv];
  Real s = moist.svp(xmol[idn])/xmol[ipr];

  Real g = 1. - xv;
  Real s1 = s/(1. - s);
  Real rate = xv - s1*g;

  // condensate at most xv vapor
  if (rate > 0.) rate = std::min(rate, xv);

  return rate;
}

#endif
