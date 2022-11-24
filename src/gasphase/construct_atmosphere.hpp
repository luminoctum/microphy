#ifndef CONSTRUCT_ATMOSPHERE_HPP
#define CONSTRUCT_ATMOSPHERE_HPP

// C/C++ header
#include <vector>
#include <cmath>

// canoe header
#include <configure.hpp>

class GasThermo;

void rk4_integrate_z(Real q[], int isat[],
  GasThermo const thermo[], Molecule const moist[],
  Real gamma, Real g_ov_Rd, Real dz, int method, Real adTdz = 0.);

void rk4_integrate_z_adaptive(Real q[], int isat[],
  GasThermo const thermo[], Molecule const moist[],
  Real gamma, Real g_ov_Rd, Real dz, Real ftol, int &niter, int method, Real adTdz = 0.);

void rk4_integrate_lnp(Real q[], int isat[],
  GasThermo const thermo[], Molecule const moist[],
  Real gamma, Real dlnp, int method, Real rdlnTdlnP = 1.);

void rk4_integrate_lnp_adaptive(Real q[], int isat[],
  GasThermo const thermo[], Molecule const moist[],
  Real gamma, Real dlnp, Real ftol, int &niter, int method, Real rdlnTdlnP = 1.);

#endif
