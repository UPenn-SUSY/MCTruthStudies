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
    if (mc_index < 0) {
      if (verbose)
        std::cout << "mc_index = " << mc_index << " -- assigning default value of -1\n";
      return -1;
    }

    int original_pdgid = mc_pdg_id->at(mc_index);
    int current_index = mc_index;
    int mother = 0;

    if (verbose) {
      std::cout << "\noriginal pdgid: " << original_pdgid
                << " -- current index: " << current_index
                << " -- mother: " << mother
                << " -- size of mc_parent_index: " << mc_parent_index->size()
                << "\n";
      for (size_t parent_it = 0; parent_it != mc_parent_index->size() ; ++parent_it) {
        std::cout << "\tparent_it: " << parent_it
                  << " num parents: " << mc_parent_index->at(parent_it).size()
                  << "\n";
      }
    }

    // this counter is just to prevent infinite loops. mother == 0 check should
    // get us out of this loop!
    for (unsigned int itr = 0; itr != 100 && mother == 0; ++itr) {
      size_t num_parents = mc_parent_index->at(current_index).size();
      if (verbose) {
        std::cout << "\titr: " << itr
                  << " -- current index: " << current_index
                  << " -- num_parents: " << num_parents
                  << "\n";
      }

      for (size_t it_parent = 0; it_parent != num_parents; ++it_parent) {
        int parent_index = mc_parent_index->at(current_index).at(it_parent);
        int parent_pdgid = mc_pdg_id->at(parent_index);

        if (verbose) {
          std::cout << "\t\tit parent: " << it_parent
                    << " -- parent index: " << parent_index
                    << " -- parent_pdgid: " << parent_pdgid
                    << "\n";
        }

        if (parent_pdgid == original_pdgid) {
          if (verbose) {
            std::cout << "\t\t\tparent pdg id == original_pdgid -- update current index\n";
          }
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

  // -----------------------------------------------------------------------------
  int getImmediateParentIndex( int mc_index
                             , const std::vector<int>* mc_pdg_id
                             , const std::vector<std::vector<int> >* mc_parent_index
                             , bool verbose
                             )
  {
    if (mc_index < 0) {
      if (verbose)
        std::cout << "mc_index = " << mc_index << " -- assigning default value of -1\n";
      return -1;
    }

    if (mc_parent_index->at(mc_index).size() > 0) {
      return mc_parent_index->at(mc_index).at(0);
    }
    return -1;
  }

  // -----------------------------------------------------------------------------
  int doDrMatchForParent( int mc_index
                        , const std::vector<int>* mc_pdg_id
                        , const std::vector<int>* mc_status_code
                        , const std::vector<float>* mc_eta
                        , const std::vector<float>* mc_phi
                        , int limit_status_code
                        , float dr_threshold
                        )
  {
    int start_pdgid = mc_pdg_id->at(mc_index);
    float start_eta = mc_eta->at(mc_index);
    float start_phi = mc_phi->at(mc_index);
    float dr_threshold_2 = dr_threshold*dr_threshold;
    int closest_mc_index = -1;

    size_t mc_n = mc_pdg_id->size();
    for (size_t mc_it = 0; mc_it != mc_n; ++mc_it) {
      // Don't match to same mc particle!
      if (mc_it == mc_index) continue;

      // if the pdg id of the start particle and the particle in the loop don't agree, skip
      if (start_pdgid != mc_pdg_id->at(mc_it)) continue;

      // if we want to limit the search to certain status codes, and this does not match, skip
      if (  limit_status_code > 0
         && mc_status_code->at(mc_it) != limit_status_code
         )
        continue;

      // calculate deta/dphi
      float deta = (start_eta - mc_eta->at(mc_it));
      float dphi = fabs(start_phi - mc_phi->at(mc_it));
      while (dphi > +3.14159) dphi -= 2*3.14159;
      while (dphi < -3.14159) dphi += 2*3.14159;

      // Find dR^2 and check if it is below threshold - if yes, update threshold and the closest mc index
      float dr2 = (deta*deta + dphi*dphi);
      if (dr2 < dr_threshold_2) {
        dr_threshold_2 = dr2;
        closest_mc_index = mc_it;
      }
    }

    return closest_mc_index;
  }
}
