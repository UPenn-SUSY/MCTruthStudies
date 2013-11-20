#ifndef PARENT_FINDER_H
#define PARENT_FINDER_H

#include <vector>

namespace TruthRecordHelpers
{
  // given a barcode, find the corresponding index within the mc truth block
  int getIndexFromBarcode( int barcode
                         , const std::vector<int>* mc_barcode
                         , bool verbose
                         );

  // given the particle barcode, get the parent index within the mc truth block
  int getParentIndexFromBarcode( int barcode
                               , const std::vector<int>* mc_barcode
                               , const std::vector<int>* mc_pdg_id
                               , const std::vector<std::vector<int> >* mc_parent_index
                               , bool verbose = false
                               );

  // from the particle barcode, get the parent pdg
  int getParentPdgIdFromBarcode( int barcode
                               , const std::vector<int>* mc_barcode
                               , const std::vector<int>* mc_pdg_id
                               , const std::vector<std::vector<int> >* mc_parent_index
                               , bool verbose = false
                               );

  // From the particle barcode, get the parent barcode
  int getParentBarcodeFromBarcode( int barcode
                                 , const std::vector<int>* mc_barcode
                                 , const std::vector<int>* mc_pdg_id
                                 , const std::vector<std::vector<int> >* mc_parent_index
                                 , bool verbose = false
                                 );

  // from the particle's mc index, get the parent index
  int getParentIndex( int mc_index
                    , const std::vector<int>* mc_mc_pdg_id
                    , const std::vector<std::vector<int> >* mc_parent_index
                    , bool verbose = false
                    );
  // from the particle's mc index, get the parent pdgid
  int getParentPdgId( int mc_index
                    , const std::vector<int>* mc_mc_pdg_id
                    , const std::vector<std::vector<int> >* mc_parent_index
                    , bool verbose = false
                    );
  // from the particle's mc index, get the parent barcode
  int getParentBarcode( int mc_index
                      , const std::vector<int>* mc_barcode
                      , const std::vector<int>* mc_mc_pdg_id
                      , const std::vector<std::vector<int> >* mc_parent_index
                      , bool verbose = false
                      );
}

#endif
