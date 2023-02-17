#include <configure.hpp>
#include <cliutils/vectorize.hpp>
#include "solution_group.hpp"

SolutionGroup::SolutionGroup(std::string kinetics_file, std::string solstrs):
  m_kk(0), m_iscalar(Globals::next_scalar_index)
{
  std::vector<std::string> names = Vectorize<std::string>(solstrs.c_str(), ", ");

  for (size_t i = 0; i < names.size(); ++i) {
    m_solutions.push_back(Cantera::newSolution(kinetics_file, names[i]));
    m_solutionsByName[names[i]] = m_solutions[i];
    m_solutionIndex[names[i]] = i;
    m_start.push_back(m_kk);
    m_kk += m_solutions[i]->thermo()->nSpecies();
  }
}
