#ifndef PARENT_FINDER_H
#define PARENT_FINDER_H

#include <vector>

namespace TruthRecordHelpers
{
  int getParentPdgIdFromBarcode( int barcode
                               , const std::vector<int>* mc_barcode
                               , const std::vector<int>* mc_pdg_id
                               , const std::vector<std::vector<int> >* mc_parent_index
                               , bool verbose = false
                               );

  int getParentPdgId( int mc_index
                    , const std::vector<int>* mc_mc_pdg_id
                    , const std::vector<std::vector<int> >* mc_parent_index
                    , bool verbose = false
                    );
}

#endif
