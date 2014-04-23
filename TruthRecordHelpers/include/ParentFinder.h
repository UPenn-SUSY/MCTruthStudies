#ifndef PARENT_FINDER_H
#define PARENT_FINDER_H

#include <vector>

// ==============================================================================
namespace TruthRecordHelpers
{
  // given a barcode, find the corresponding index within the mc truth block
  int getIndexFromBarcode( int barcode
                         , const std::vector<int>* mc_barcode
                         , bool verbose = false
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

  // from the index of this particle when it is created (before any scattering)
  int getInitialIndex( int mc_index
                     , const std::vector<int>* mc_mc_pdg_id
                     , const std::vector<std::vector<int> >* mc_parent_index
                     , bool follow_tau_parent = false
                     , bool verbose = false
                     );
  // from the particle's mc index, get the parent index
  int getParentIndex( int mc_index
                    , const std::vector<int>* mc_mc_pdg_id
                    , const std::vector<std::vector<int> >* mc_parent_index
                    , bool follow_tau_parent = false
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

  // from the particle's mc index, get the immediate parent index
  int getImmediateParentIndex( int mc_index
                             , const std::vector<int>* mc_mc_pdg_id
                             , const std::vector<std::vector<int> >* mc_parent_index
                             , bool verbose = false
                             );

  // search through the truth record for the particle with the same pdgID which is closest in dR
  // return the mc index of the found particle
  int doDrMatchForParent( int mc_index
                        , const std::vector<int>* mc_pdg_id
                        , const std::vector<int>* mc_status_code
                        , const std::vector<float>* mc_eta
                        , const std::vector<float>* mc_phi
                        , int limit_status_code = -1
                        , float dr_threshold = 0.1
                        );


}

#endif
