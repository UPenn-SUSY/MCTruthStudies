#include <iostream>
#include <vector>
#include "include/ParentFinder.h"

#include "TPython.h"

// =============================================================================
namespace TruthRecordHelpers
{
  // -----------------------------------------------------------------------------
  int getParentPdgIdFromBarcode( int barcode
                               , const std::vector<int>* mc_barcode
                               , const std::vector<int>* mc_pdg_id
                               , const std::vector<std::vector<int> >* mc_parent_index
                               , bool verbose
                               )
  {
    size_t term = mc_barcode->size();
    for (size_t mc_index= 0; mc_index != term; ++mc_index) {
      if (mc_barcode->at(mc_index) == barcode) {
        if (verbose) {
          std::cout << "found a barcode match at index " << mc_index << " of " << mc_barcode->size() << ":\n";
          std::cout << "\tbarcodes: " << mc_barcode->at(mc_index) << " = " << barcode << "\n";
        }
        return getParentPdgId( mc_index
                             , mc_pdg_id
                             , mc_parent_index
                             );
      }
    }
    return 0;
  }

  // -----------------------------------------------------------------------------
  int getParentPdgId( int mc_index
                    , const std::vector<int>* mc_pdg_id
                    , const std::vector<std::vector<int> >* mc_parent_index
                    , bool verbose
                    )
  {
    int original_pdgid = mc_pdg_id->at(mc_index);
    int current_index = mc_index;
    int mother = 0;

    for (unsigned int itr = 0; itr != 100 && mother == 0; ++itr) {
      size_t num_parents = mc_parent_index->at(current_index).size();

      for (size_t it_parent = 0; it_parent != num_parents; ++it_parent) {
        int parent_index = mc_parent_index->at(current_index).at(it_parent);
        int parent_pdgid = mc_pdg_id->at(parent_index);

        if (parent_pdgid == original_pdgid) {
          current_index = parent_index;
          break;
        }
        else {
          mother = parent_pdgid;
        }
      }
    }
    return mother;
  }
}
