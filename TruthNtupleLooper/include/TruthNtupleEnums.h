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
  const std::string FlavorChannelStrings[] = { "fc_all"
                                             , "fc_ee"
                                             , "fc_mm"
                                             , "fc_em"
                                             , "fc_me"
                                             , "fc_none"
                                             , "fc_n"
                                             };
}

#endif
