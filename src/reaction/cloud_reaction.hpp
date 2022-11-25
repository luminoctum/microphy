#ifndef CLOUD_REACTION_HPP
#define CLOUD_REACTION_HPP

#include <string>
#include <configure.hpp>

class GasPhase;
class CloudPhase;

using ReactionFunc = Real (*)(Real qmol[], int const indx[],
  GasPhase const *pgas, CloudPhase const *pcld);

struct CloudReaction {
  CloudReaction(std::string rstr,
    GasPhase const *pgas, CloudPhase const *pcld);

  int           indx[MaxReactants];
  ReactionFunc  fptr;
};

#endif
