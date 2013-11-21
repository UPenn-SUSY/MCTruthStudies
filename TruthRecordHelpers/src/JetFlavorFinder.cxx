#include <iostream>
#include <vector>
#include <math.h>
#include "include/JetFlavorFinder.h"
#include "TruthNtupleLooper/include/Calculators.h"

// =============================================================================
namespace TruthRecordHelpers
{
  // -----------------------------------------------------------------------------
  int isBJet( const float my_jet_eta
            , const float my_jet_phi
            , const std::vector<int>* mc_pdgId
            , const std::vector<int>* /*mc_status*/
            , const std::vector<int>* /*mc_barcode*/
            , const std::vector<float>* mc_pt
            , const std::vector<float>* mc_eta
            , const std::vector<float>* mc_phi
            , bool /*verbose*/
            )
  {
    double delta_r_b   = 999.;
    // double delta_r_c   = 999.;
    // double delta_r_tau = 999.;
    double delta_r     = 999.;

    int index_b = -1;
    // int barcode_b = 0;
    // int barcode_c = 0;

    const double pt_cut_value = 5.e3;
    const float delta_r_cut = 0.30;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // loop over mc record
    size_t num_mc_records = mc_pdgId->size();
    for (size_t mc_it = 0; mc_it != num_mc_records; ++mc_it) {
      // check that mc particle passes pt cut
      if (  mc_pt->at(mc_it) < pt_cut_value
          && abs(mc_pdgId->at(mc_it)) != 5
         ) {
        continue;
      }
      // find dR between jet and mc particle
      delta_r = TruthNtuple::deltaR( my_jet_eta, my_jet_phi
                                   , mc_eta->at(mc_it), mc_phi->at(mc_it)
                                   );

      // if this particle is a b-quark, and it is closer than previous b quarks,
      // make this the selected b-quark
      if (abs(mc_pdgId->at(mc_it)) == 5 && delta_r < delta_r_b) {
        delta_r_b = delta_r;
        index_b = mc_it;
      }

      /*
      // leave this here for when we do full flavor finding
      if (abs(mc_pdgId->at(mc_it)) == 5 && delta_r < delta_r_b) {
      // delta_r_b = delta_r;
      // barcode_b = mc_barcode->at(mc_it);
      }
      */
    }

    // if there is a b-quark close to our jet, return b-quark index
    if (delta_r_b < delta_r_cut)
      return index_b;

    return -1;
  }
}
