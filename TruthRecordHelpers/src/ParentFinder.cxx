#include <iostream>
#include <vector>
#include <math.h>
#include "include/ParentFinder.h"

#include "TPython.h"

// =============================================================================
namespace TruthRecordHelpers
{
  // -----------------------------------------------------------------------------
  int getIndexFromBarcode( int barcode
                         , const std::vector<int>* mc_barcode
                         , bool verbose
                         )
  {
    // loop over mc truth block
    size_t term = mc_barcode->size();
    for (size_t mc_index= 0; mc_index != term; ++mc_index) {
      if (mc_barcode->at(mc_index) == barcode) {
        if (verbose) {
          std::cout << "found a barcode match at index " << mc_index
                    << " of " << mc_barcode->size() << ":\n";
          std::cout << "\tbarcodes: " << mc_barcode->at(mc_index)
                    << " = " << barcode << "\n";
        }
        return mc_index;
      }
    }

    return -1;
  }

  // -----------------------------------------------------------------------------
  int getParentIndexFromBarcode( int barcode
                               , const std::vector<int>* mc_barcode
                               , const std::vector<int>* mc_pdg_id
                               , const std::vector<std::vector<int> >* mc_parent_index
                               , bool verbose
                               )
  {
    // first, get mc index
    int mc_index = getIndexFromBarcode(barcode, mc_barcode, verbose);

    // if valid index, find the parent pdg id
    if (mc_index >= 0)
      return getParentIndex(mc_index, mc_pdg_id, mc_parent_index);

    // else, return 0
    return 0;
  }

  // -----------------------------------------------------------------------------
  int getParentPdgIdFromBarcode( int barcode
                               , const std::vector<int>* mc_barcode
                               , const std::vector<int>* mc_pdg_id
                               , const std::vector<std::vector<int> >* mc_parent_index
                               , bool verbose
                               )
  {
    // first, get mc index
    int mc_index = getIndexFromBarcode(barcode, mc_barcode, verbose);

    // if valid index, find the parent pdg id
    if (mc_index >= 0)
      return getParentPdgId(mc_index, mc_pdg_id, mc_parent_index);

    // else, return 0
    return 0;
  }

  // -----------------------------------------------------------------------------
  int getParentBarcodeFromBarcode( int barcode
                                 , const std::vector<int>* mc_barcode
                                 , const std::vector<int>* mc_pdg_id
                                 , const std::vector<std::vector<int> >* mc_parent_index
                                 , bool verbose
                                 )
  {
    // first, get mc index
    int mc_index = getIndexFromBarcode(barcode, mc_barcode, verbose);

    // if valid index, find the parent pdg id
    if (mc_index >= 0)
      return getParentBarcode(mc_index, mc_barcode, mc_pdg_id, mc_parent_index);

    // else, return 0
    return 0;
  }

  // -----------------------------------------------------------------------------
  int getInitialIndex( int mc_index
                     , const std::vector<int>* mc_pdg_id
                     , const std::vector<std::vector<int> >* mc_parent_index
                     , bool follow_tau_parent
                     , bool verbose
                     )
  {
    if (mc_index < 0) return -1;

    int original_pdgid = mc_pdg_id->at(mc_index);
    int current_index = mc_index;
    int mother = 0;

    if (verbose) {
      std::cout << "looking for initial index: "
                << "\n\tfinal index: " << mc_index
                << "\n\toriginal pdgid: " << original_pdgid
                << "\n";
    }

    // this counter is just to prevent infinite loops. mother == 0 check should
    // get us out of this loop!
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
          if (verbose) {
            std::cout << "found initial index! -- original_pdgid" << original_pdgid
                      << " -- mother pdg id: " << parent_pdgid
                      << " -- initial index: " << current_index
                      << " -- mother index: " << parent_index
                      << "\n";
          }
          if (follow_tau_parent && fabs(parent_pdgid) == 15) {
            if (verbose) {
              std::cout << "The parent is a tau: Stepping back a level to find the tau's parent\n";
            }
            return getInitialIndex(parent_index, mc_pdg_id, mc_parent_index, follow_tau_parent, verbose);
          }
          return current_index;
        }
      }
    }
    if (verbose) {
      std::cout << "did not find initial index -- original_pdgid: "
                << original_pdgid << "\n";
    }
    return -1;
  }

  // -----------------------------------------------------------------------------
  int getParentIndex( int mc_index
                    , const std::vector<int>* mc_pdg_id
                    , const std::vector<std::vector<int> >* mc_parent_index
                    , bool follow_tau_parent
                    , bool verbose
                    )
  {
    if (mc_index < 0) return -1;

    int original_pdgid = mc_pdg_id->at(mc_index);
    int current_index = mc_index;
    int mother = 0;

    // this counter is just to prevent infinite loops. mother == 0 check should
    // get us out of this loop!
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
          if (verbose) {
            std::cout << "found mother! -- original_pdgid " << original_pdgid
                      << " -- mother pdg id: " << parent_pdgid
                      << " -- mother index: " << parent_index
                      << " -- follow tau parent: " << follow_tau_parent
                      << "\n";
          }
          if (follow_tau_parent && fabs(parent_pdgid) == 15) {
            if (verbose) {
              std::cout << "The parent is a tau: Stepping back a level to find the tau's parent\n";
            }
            return getParentIndex(parent_index, mc_pdg_id, mc_parent_index, follow_tau_parent, verbose);
          }
          // mother = parent_pdgid;
          return parent_index;
        }
      }
    }
    if (verbose) {
      std::cout << "did not find mother -- original_pdgid: "
                << original_pdgid << "\n";
    }
    return -1;
  }

  // -----------------------------------------------------------------------------
  int getParentPdgId( int mc_index
                    , const std::vector<int>* mc_pdg_id
                    , const std::vector<std::vector<int> >* mc_parent_index
                    , bool verbose
                    )
  {
    int parent_index = getParentIndex( mc_index
                                     , mc_pdg_id
                                     , mc_parent_index
                                     , verbose
                                     );
    if (parent_index >= 0) {
      return mc_pdg_id->at(parent_index);
    }
    return 0;
  }

  // -----------------------------------------------------------------------------
  int getParentBarcode( int mc_index
                      , const std::vector<int>* mc_barcode
                      , const std::vector<int>* mc_pdg_id
                      , const std::vector<std::vector<int> >* mc_parent_index
                      , bool verbose
                      )
  {
    int parent_index = getParentIndex( mc_index
                                     , mc_pdg_id
                                     , mc_parent_index
                                     , verbose
                                     );
    return mc_barcode->at(parent_index);
  }
}
