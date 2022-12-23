// C/C++
#include <iostream>

// canoe
#include "gas_thermo.hpp"

std::ostream& operator<<(std::ostream& os, GasThermo const& my)
{
  std::cout << "mu_ratio = " << my.mu_ratio << ", ";
  std::cout << "cp_ratio = " << my.cp_ratio << ", ";
  std::cout << "cv_ratio = " << my.cv_ratio << ", ";
  std::cout << "beta = " << my.beta << ", ";
  std::cout << "delta = " << my.delta << ", ";
  std::cout << "t3 = " << my.t3 << ", ";
  std::cout << "p3 = " << my.p3;
  return os;
}
