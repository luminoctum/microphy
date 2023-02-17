#ifndef SOLUTION_GROUP_HPP
#define SOLUTION_GROUP_HPP

#include <memory>
#include <map>
#include <cantera/base/Solution.h>
#include <cantera/thermo/ThermoPhase.h>

class SolutionGroup {
public:
  SolutionGroup(std::string kinetics_file, std::string solstrs);

  //! get the Solution object for a solution in a group by index
  std::shared_ptr<Cantera::Solution> solution(int i) {
    return m_solutions[i];
  }

  std::shared_ptr<Cantera::Solution> solution(std::string const& name) {
    return m_solutionsByName[name];
  }

  size_t nSolutions() const {
    return m_solutions.size();
  }

  size_t nTotalSpecies() const {
    return m_kk;
  }

  size_t solutionIndex(std::string const& name) const {
    if (m_solutionIndex.find(name) == m_solutionIndex.end()) {
      return Cantera::npos;
    } else {
      return m_solutionIndex.at(name);
    }
  }

  size_t scalarIndexFirst() const {
    return m_iscalar;
  }

  size_t scalarIndexLast() const {
    return m_iscalar + nTotalSpecies() - 1;
  }

  size_t speciesScalarIndexStart(std::string const& phase) const
  {
    size_t i = solutionIndex(phase);
    if (i != Cantera::npos)
      return m_iscalar + m_start[i];
    else
      return Cantera::npos;
  }

  size_t speciesScalarIndex(std::string const& name) const
  {
    for (size_t i = 0; i < m_solutions.size(); ++i) {
      size_t j = m_solutions[i]->thermo()->speciesIndex(name);
      if (j != Cantera::npos) {
        return m_iscalar + m_start[i] + j;
      }
    }
    return Cantera::npos;
  }

  size_t speciesScalarIndex(std::string const& phase, std::string const& name) const
  {
    size_t i = solutionIndex(phase);
    if (i != Cantera::npos)
      return m_iscalar + m_start[i] + m_solutions[i]->thermo()->speciesIndex(name);
    else
      return Cantera::npos;
  }

protected:
  //! total number of species in the solution group
  size_t m_kk;

  //! index mapping to athena passive scalar fields
  size_t m_iscalar;

  /**
   * m_start is a vector of integers specifying the beginning position for the
   * species vector for the n'th solution group in the condensates class.
   */
  std::vector<size_t> m_start;

  std::vector<std::shared_ptr<Cantera::Solution>> m_solutions;

  //! solutions, for access by name
  std::map<std::string, std::shared_ptr<Cantera::Solution>> m_solutionsByName;

  //! solution indices, for access by name
  std::map<std::string, size_t> m_solutionIndex;
};

#endif
