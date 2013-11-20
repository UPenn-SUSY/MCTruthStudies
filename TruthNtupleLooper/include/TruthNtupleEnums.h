#ifndef TRUTHNTUPLEENUMS_H
#define TRUTHNTUPLEENUMS_H

#include <string>

// =============================================================================
namespace TruthNtuple
{
  enum FLAVOR_CHANNEL { FLAVOR_ALL
                      , FLAVOR_EE
                      , FLAVOR_MM
                      , FLAVOR_EM
                      , FLAVOR_ME
                      , FLAVOR_NONE
                      , FLAVOR_N
                      };
  const std::string FlavorChannelStrings[] = { "flavor_all"
                                             , "flavor_ee"
                                             , "flavor_mm"
                                             , "flavor_em"
                                             , "flavor_me"
                                             , "flavor_none"
                                             , "flavor_n"
                                             };
}

#endif
