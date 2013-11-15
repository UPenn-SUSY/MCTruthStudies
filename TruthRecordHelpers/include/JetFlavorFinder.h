#ifndef JET_FLAVOR_FINDER_H
#define JET_FLAVOR_FINDER_H

#include <vector>

// =============================================================================
namespace TruthNtuple
{
  class Jet;
}

// =============================================================================
namespace TruthRecordHelpers
{
  // -----------------------------------------------------------------------------
  bool isBJet( const float my_jet_eta
             , const float my_jet_phi
             , const std::vector<int>* mc_pdgId
             , const std::vector<int>* /*mc_status*/
             , const std::vector<int>* /*mc_barcode*/
             , const std::vector<float>* mc_pt
             , const std::vector<float>* mc_eta
             , const std::vector<float>* mc_phi
             , bool verbose = false
             );
}

#endif
