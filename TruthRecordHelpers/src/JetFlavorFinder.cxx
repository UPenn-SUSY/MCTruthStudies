#include <iostream>
#include <vector>
#include <math.h>
#include "TruthRecordHelpers/include/JetFlavorFinder.h"

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
            , bool verbose
            )
  {
    // double delta_r_b   = 999.;
    // double delta_r_c   = 999.;
    // double delta_r_tau = 999.;
    // double delta_r     = 999.;

    // int barcode_b = 0;
    // int barcode_c = 0;

    const double pt_cut_value = 10.e3;
    const float delta_r_cut = 0.30;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    size_t num_mc_records = mc_pdgId->size();
    for (size_t mc_it = 0; mc_it != num_mc_records; ++mc_it) {
      if (  mc_pt->at(mc_it) < pt_cut_value
          && abs(mc_pdgId->at(mc_it)) != 5
         ) {
        continue;
      }
      double delta_eta = fabs(fabs(mc_eta->at(mc_it)) - fabs(my_jet_eta));
      double delta_phi = fabs(fabs(mc_phi->at(mc_it)) - fabs(my_jet_phi));
      while (delta_phi > 3.14159) {
        delta_phi -= 3.14159;
      }
      double delta_r = sqrt(delta_eta*delta_eta + delta_phi*delta_phi);

      if (abs(mc_pdgId->at(mc_it)) == 5 && delta_r < delta_r_cut) {
        return mc_it;
      }

      /*
      // leave this here for when we do full flavor finding
      if (abs(mc_pdgId->at(mc_it)) == 5 && delta_r < delta_r_b) {
      // delta_r_b = delta_r;
      // barcode_b = mc_barcode->at(mc_it);
      }
      */
    }

    /*
       if (delta_r_b < delta_r_cut)
       return true;
       */

    return -1;
  }
}
