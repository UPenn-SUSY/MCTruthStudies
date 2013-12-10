#ifndef BMINUSLCUTFLOW_CALCULATORS_H
#define BMINUSLCUTFLOW_CALCULATORS_H

#include <vector>
#include "TruthNtupleLooper/include/ObjectDefs.h"

namespace BMinusL
{
  // do pairing by finding combiniations of b and l which provide the minumum
  // difference in mbl
  std::vector< std::pair< TruthNtuple::Particle*
                        , TruthNtuple::Particle*
                        >
             > doMblDiffPairing( const std::vector<TruthNtuple::Particle*>& b_list
                               , const std::vector<TruthNtuple::Particle*>& l_list
                               );

  // do pairing by finding combiniations of b and l which provide the maximum
  // ratio mbl(subleading)/mbl(leading)
  std::vector< std::pair< TruthNtuple::Particle*
                        , TruthNtuple::Particle*
                        >
             > doMblRatioPairing( const std::vector<TruthNtuple::Particle*>& b_list
                                , const std::vector<TruthNtuple::Particle*>& l_list
                                );

  // do pairing by finding combiniations of b and l which provide the maximum
  // quadratic sum of the mbl values
  std::vector< std::pair< TruthNtuple::Particle*
                        , TruthNtuple::Particle*
                        >
             > doMblSqSumPairing( const std::vector<TruthNtuple::Particle*>& b_list
                                , const std::vector<TruthNtuple::Particle*>& l_list
                                );

}

#endif
