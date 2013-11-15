#ifndef TruthNtupleLooper_h
#define TruthNtupleLooper_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

#include "include/ObjectDefs.h"

// =============================================================================
class TBranch;

// =============================================================================
namespace TruthNtuple
{
  // ===========================================================================
  class TruthNtupleLooper {
    public :
      TruthNtupleLooper(TTree *tree=0);
      virtual ~TruthNtupleLooper();
      virtual Int_t    GetEntry(Long64_t entry);
      virtual Long64_t LoadTree(Long64_t entry);
      virtual void     Init(TTree *tree);
      virtual void     Loop();
      virtual Bool_t   Notify();

      virtual void clearObjects();
      virtual void constructObjects();
      virtual void processEvent();

    protected:
      std::vector<TruthNtuple::Electron> m_el_list;
      std::vector<TruthNtuple::Muon>     m_mu_list;
      std::vector<TruthNtuple::Jet>      m_jet_list;

      TTree          *fChain;   //!pointer to the analyzed TTree or TChain
      Int_t           fCurrent; //!current Tree number in a TChain

      // Declaration of leaf types
      UInt_t          RunNumber;
      UInt_t          EventNumber;
      UInt_t          timestamp;
      UInt_t          timestamp_ns;
      UInt_t          lbn;
      UInt_t          bcid;
      UInt_t          detmask0;
      UInt_t          detmask1;
      Float_t         actualIntPerXing;
      Float_t         averageIntPerXing;
      UInt_t          mc_channel_number;
      UInt_t          mc_event_number;
      Float_t         mc_event_weight;
      UInt_t          pixelFlags;
      UInt_t          sctFlags;
      UInt_t          trtFlags;
      UInt_t          larFlags;
      UInt_t          tileFlags;
      UInt_t          muonFlags;
      UInt_t          fwdFlags;
      UInt_t          coreFlags;
      UInt_t          pixelError;
      UInt_t          sctError;
      UInt_t          trtError;
      UInt_t          larError;
      UInt_t          tileError;
      UInt_t          muonError;
      UInt_t          fwdError;
      UInt_t          coreError;
      Bool_t          isSimulation;
      Bool_t          isCalibration;
      Bool_t          isTestBeam;
      Int_t           mcevt_n;
      std::vector<int>     *mcevt_signal_process_id;
      std::vector<int>     *mcevt_event_number;
      std::vector<double>  *mcevt_event_scale;
      std::vector<double>  *mcevt_alphaQCD;
      std::vector<double>  *mcevt_alphaQED;
      std::vector<int>     *mcevt_pdf_id1;
      std::vector<int>     *mcevt_pdf_id2;
      std::vector<double>  *mcevt_pdf_x1;
      std::vector<double>  *mcevt_pdf_x2;
      std::vector<double>  *mcevt_pdf_scale;
      std::vector<double>  *mcevt_pdf1;
      std::vector<double>  *mcevt_pdf2;
      std::vector<std::vector<double> > *mcevt_weight;
      std::vector<int>     *mcevt_nparticle;
      std::vector<short>   *mcevt_pileUpType;
      Int_t           mcVxn;
      std::vector<float>   *mcVxx;
      std::vector<float>   *mcVxy;
      std::vector<float>   *mcVxz;
      Int_t           vxn;
      std::vector<float>   *vxx;
      std::vector<float>   *vxy;
      std::vector<float>   *vxz;
      Int_t           jet_AntiKt4TruthJets_n;
      std::vector<float>   *jet_AntiKt4TruthJets_E;
      std::vector<float>   *jet_AntiKt4TruthJets_pt;
      std::vector<float>   *jet_AntiKt4TruthJets_m;
      std::vector<float>   *jet_AntiKt4TruthJets_eta;
      std::vector<float>   *jet_AntiKt4TruthJets_phi;
      std::vector<float>   *jet_AntiKt4TruthJets_EtaOrigin;
      std::vector<float>   *jet_AntiKt4TruthJets_PhiOrigin;
      std::vector<float>   *jet_AntiKt4TruthJets_MOrigin;
      std::vector<float>   *jet_AntiKt4TruthJets_EtaOriginEM;
      std::vector<float>   *jet_AntiKt4TruthJets_PhiOriginEM;
      std::vector<float>   *jet_AntiKt4TruthJets_MOriginEM;
      std::vector<float>   *jet_AntiKt4TruthJets_WIDTH;
      std::vector<float>   *jet_AntiKt4TruthJets_n90;
      std::vector<float>   *jet_AntiKt4TruthJets_Timing;
      std::vector<float>   *jet_AntiKt4TruthJets_LArQuality;
      std::vector<float>   *jet_AntiKt4TruthJets_nTrk;
      std::vector<float>   *jet_AntiKt4TruthJets_sumPtTrk;
      std::vector<float>   *jet_AntiKt4TruthJets_OriginIndex;
      std::vector<float>   *jet_AntiKt4TruthJets_HECQuality;
      std::vector<float>   *jet_AntiKt4TruthJets_NegativeE;
      std::vector<float>   *jet_AntiKt4TruthJets_AverageLArQF;
      std::vector<float>   *jet_AntiKt4TruthJets_YFlip12;
      std::vector<float>   *jet_AntiKt4TruthJets_YFlip23;
      std::vector<float>   *jet_AntiKt4TruthJets_BCH_CORR_CELL;
      std::vector<float>   *jet_AntiKt4TruthJets_BCH_CORR_DOTX;
      std::vector<float>   *jet_AntiKt4TruthJets_BCH_CORR_JET;
      std::vector<float>   *jet_AntiKt4TruthJets_BCH_CORR_JET_FORCELL;
      std::vector<float>   *jet_AntiKt4TruthJets_ENG_BAD_CELLS;
      std::vector<float>   *jet_AntiKt4TruthJets_N_BAD_CELLS;
      std::vector<float>   *jet_AntiKt4TruthJets_N_BAD_CELLS_CORR;
      std::vector<float>   *jet_AntiKt4TruthJets_BAD_CELLS_CORR_E;
      std::vector<float>   *jet_AntiKt4TruthJets_NumTowers;
      std::vector<float>   *jet_AntiKt4TruthJets_ootFracCells5;
      std::vector<float>   *jet_AntiKt4TruthJets_ootFracCells10;
      std::vector<float>   *jet_AntiKt4TruthJets_ootFracClusters5;
      std::vector<float>   *jet_AntiKt4TruthJets_ootFracClusters10;
      std::vector<int>     *jet_AntiKt4TruthJets_SamplingMax;
      std::vector<float>   *jet_AntiKt4TruthJets_fracSamplingMax;
      std::vector<float>   *jet_AntiKt4TruthJets_hecf;
      std::vector<float>   *jet_AntiKt4TruthJets_tgap3f;
      std::vector<int>     *jet_AntiKt4TruthJets_isUgly;
      std::vector<int>     *jet_AntiKt4TruthJets_isBadLooseMinus;
      std::vector<int>     *jet_AntiKt4TruthJets_isBadLoose;
      std::vector<int>     *jet_AntiKt4TruthJets_isBadMedium;
      std::vector<int>     *jet_AntiKt4TruthJets_isBadTight;
      std::vector<float>   *jet_AntiKt4TruthJets_emfrac;
      std::vector<float>   *jet_AntiKt4TruthJets_Offset;
      std::vector<float>   *jet_AntiKt4TruthJets_EMJES;
      std::vector<float>   *jet_AntiKt4TruthJets_EMJES_EtaCorr;
      std::vector<float>   *jet_AntiKt4TruthJets_EMJESnooffset;
      std::vector<float>   *jet_AntiKt4TruthJets_GCWJES;
      std::vector<float>   *jet_AntiKt4TruthJets_GCWJES_EtaCorr;
      std::vector<float>   *jet_AntiKt4TruthJets_CB;
      std::vector<float>   *jet_AntiKt4TruthJets_LCJES;
      std::vector<float>   *jet_AntiKt4TruthJets_emscale_E;
      std::vector<float>   *jet_AntiKt4TruthJets_emscale_pt;
      std::vector<float>   *jet_AntiKt4TruthJets_emscale_m;
      std::vector<float>   *jet_AntiKt4TruthJets_emscale_eta;
      std::vector<float>   *jet_AntiKt4TruthJets_emscale_phi;
      Int_t           jet_AntiKt4TopoNewEM_n;
      std::vector<float>   *jet_AntiKt4TopoNewEM_E;
      std::vector<float>   *jet_AntiKt4TopoNewEM_pt;
      std::vector<float>   *jet_AntiKt4TopoNewEM_m;
      std::vector<float>   *jet_AntiKt4TopoNewEM_eta;
      std::vector<float>   *jet_AntiKt4TopoNewEM_phi;
      std::vector<float>   *jet_AntiKt4TopoNewEM_EtaOrigin;
      std::vector<float>   *jet_AntiKt4TopoNewEM_PhiOrigin;
      std::vector<float>   *jet_AntiKt4TopoNewEM_MOrigin;
      std::vector<float>   *jet_AntiKt4TopoNewEM_EtaOriginEM;
      std::vector<float>   *jet_AntiKt4TopoNewEM_PhiOriginEM;
      std::vector<float>   *jet_AntiKt4TopoNewEM_MOriginEM;
      std::vector<float>   *jet_AntiKt4TopoNewEM_WIDTH;
      std::vector<float>   *jet_AntiKt4TopoNewEM_n90;
      std::vector<float>   *jet_AntiKt4TopoNewEM_Timing;
      std::vector<float>   *jet_AntiKt4TopoNewEM_LArQuality;
      std::vector<float>   *jet_AntiKt4TopoNewEM_nTrk;
      std::vector<float>   *jet_AntiKt4TopoNewEM_sumPtTrk;
      std::vector<float>   *jet_AntiKt4TopoNewEM_OriginIndex;
      std::vector<float>   *jet_AntiKt4TopoNewEM_HECQuality;
      std::vector<float>   *jet_AntiKt4TopoNewEM_NegativeE;
      std::vector<float>   *jet_AntiKt4TopoNewEM_AverageLArQF;
      std::vector<float>   *jet_AntiKt4TopoNewEM_YFlip12;
      std::vector<float>   *jet_AntiKt4TopoNewEM_YFlip23;
      std::vector<float>   *jet_AntiKt4TopoNewEM_BCH_CORR_CELL;
      std::vector<float>   *jet_AntiKt4TopoNewEM_BCH_CORR_DOTX;
      std::vector<float>   *jet_AntiKt4TopoNewEM_BCH_CORR_JET;
      std::vector<float>   *jet_AntiKt4TopoNewEM_BCH_CORR_JET_FORCELL;
      std::vector<float>   *jet_AntiKt4TopoNewEM_ENG_BAD_CELLS;
      std::vector<float>   *jet_AntiKt4TopoNewEM_N_BAD_CELLS;
      std::vector<float>   *jet_AntiKt4TopoNewEM_N_BAD_CELLS_CORR;
      std::vector<float>   *jet_AntiKt4TopoNewEM_BAD_CELLS_CORR_E;
      std::vector<float>   *jet_AntiKt4TopoNewEM_NumTowers;
      std::vector<float>   *jet_AntiKt4TopoNewEM_ootFracCells5;
      std::vector<float>   *jet_AntiKt4TopoNewEM_ootFracCells10;
      std::vector<float>   *jet_AntiKt4TopoNewEM_ootFracClusters5;
      std::vector<float>   *jet_AntiKt4TopoNewEM_ootFracClusters10;
      std::vector<int>     *jet_AntiKt4TopoNewEM_SamplingMax;
      std::vector<float>   *jet_AntiKt4TopoNewEM_fracSamplingMax;
      std::vector<float>   *jet_AntiKt4TopoNewEM_hecf;
      std::vector<float>   *jet_AntiKt4TopoNewEM_tgap3f;
      std::vector<int>     *jet_AntiKt4TopoNewEM_isUgly;
      std::vector<int>     *jet_AntiKt4TopoNewEM_isBadLooseMinus;
      std::vector<int>     *jet_AntiKt4TopoNewEM_isBadLoose;
      std::vector<int>     *jet_AntiKt4TopoNewEM_isBadMedium;
      std::vector<int>     *jet_AntiKt4TopoNewEM_isBadTight;
      std::vector<float>   *jet_AntiKt4TopoNewEM_emfrac;
      std::vector<float>   *jet_AntiKt4TopoNewEM_Offset;
      std::vector<float>   *jet_AntiKt4TopoNewEM_EMJES;
      std::vector<float>   *jet_AntiKt4TopoNewEM_EMJES_EtaCorr;
      std::vector<float>   *jet_AntiKt4TopoNewEM_EMJESnooffset;
      std::vector<float>   *jet_AntiKt4TopoNewEM_GCWJES;
      std::vector<float>   *jet_AntiKt4TopoNewEM_GCWJES_EtaCorr;
      std::vector<float>   *jet_AntiKt4TopoNewEM_CB;
      std::vector<float>   *jet_AntiKt4TopoNewEM_LCJES;
      std::vector<float>   *jet_AntiKt4TopoNewEM_emscale_E;
      std::vector<float>   *jet_AntiKt4TopoNewEM_emscale_pt;
      std::vector<float>   *jet_AntiKt4TopoNewEM_emscale_m;
      std::vector<float>   *jet_AntiKt4TopoNewEM_emscale_eta;
      std::vector<float>   *jet_AntiKt4TopoNewEM_emscale_phi;
      Int_t           mc_n;
      std::vector<float>   *mc_E;
      std::vector<float>   *mc_pt;
      std::vector<float>   *mc_m;
      std::vector<float>   *mc_eta;
      std::vector<float>   *mc_phi;
      std::vector<float>   *mc_px;
      std::vector<float>   *mc_py;
      std::vector<float>   *mc_pz;
      std::vector<int>     *mc_status;
      std::vector<int>     *mc_barcode;
      std::vector<int>     *mc_pdgId;
      std::vector<float>   *mc_charge;
      std::vector<float>   *mc_vx_x;
      std::vector<float>   *mc_vx_y;
      std::vector<float>   *mc_vx_z;
      std::vector<int>     *mc_vx_barcode;
      std::vector<std::vector<int> > *mc_child_index;
      std::vector<std::vector<int> > *mc_parent_index;
      Float_t         MET_Truth_NonInt_etx;
      Float_t         MET_Truth_NonInt_ety;
      Float_t         MET_Truth_Int_etx;
      Float_t         MET_Truth_Int_ety;
      Float_t         MET_Truth_IntCentral_etx;
      Float_t         MET_Truth_IntCentral_ety;
      Float_t         MET_Truth_IntFwd_etx;
      Float_t         MET_Truth_IntFwd_ety;
      Float_t         MET_Truth_IntOutCover_etx;
      Float_t         MET_Truth_IntOutCover_ety;
      Float_t         MET_Truth_IntMuons_etx;
      Float_t         MET_Truth_IntMuons_ety;
      Int_t           el_n;
      std::vector<float>   *el_E;
      std::vector<float>   *el_pt;
      std::vector<float>   *el_m;
      std::vector<float>   *el_px;
      std::vector<float>   *el_py;
      std::vector<float>   *el_pz;
      std::vector<float>   *el_eta;
      std::vector<float>   *el_phi;
      std::vector<int>     *el_status;
      std::vector<int>     *el_barcode;
      std::vector<int>     *el_charge;
      Int_t           mu_muid_n;
      std::vector<float>   *mu_muid_E;
      std::vector<float>   *mu_muid_pt;
      std::vector<float>   *mu_muid_m;
      std::vector<float>   *mu_muid_px;
      std::vector<float>   *mu_muid_py;
      std::vector<float>   *mu_muid_pz;
      std::vector<float>   *mu_muid_eta;
      std::vector<float>   *mu_muid_phi;
      std::vector<int>     *mu_muid_status;
      std::vector<int>     *mu_muid_barcode;
      std::vector<int>     *mu_muid_charge;
      Int_t           mu_staco_n;
      std::vector<float>   *mu_staco_E;
      std::vector<float>   *mu_staco_pt;
      std::vector<float>   *mu_staco_m;
      std::vector<float>   *mu_staco_px;
      std::vector<float>   *mu_staco_py;
      std::vector<float>   *mu_staco_pz;
      std::vector<float>   *mu_staco_eta;
      std::vector<float>   *mu_staco_phi;
      std::vector<int>     *mu_staco_status;
      std::vector<int>     *mu_staco_barcode;
      std::vector<int>     *mu_staco_charge;
      Int_t           tau_n;
      std::vector<float>   *tau_pt;
      std::vector<float>   *tau_m;
      std::vector<float>   *tau_eta;
      std::vector<float>   *tau_phi;
      std::vector<int>     *tau_status;
      std::vector<int>     *tau_barcode;
      std::vector<int>     *tau_charge;
      Int_t           ph_n;
      std::vector<float>   *ph_E;
      std::vector<float>   *ph_pt;
      std::vector<float>   *ph_m;
      std::vector<float>   *ph_px;
      std::vector<float>   *ph_py;
      std::vector<float>   *ph_pz;
      std::vector<float>   *ph_eta;
      std::vector<float>   *ph_phi;
      std::vector<int>     *ph_status;
      std::vector<int>     *ph_barcode;
      Int_t           el_truth_n;
      std::vector<float>   *el_truth_E;
      std::vector<float>   *el_truth_pt;
      std::vector<float>   *el_truth_m;
      std::vector<float>   *el_truth_px;
      std::vector<float>   *el_truth_py;
      std::vector<float>   *el_truth_pz;
      std::vector<float>   *el_truth_eta;
      std::vector<float>   *el_truth_phi;
      std::vector<int>     *el_truth_status;
      std::vector<int>     *el_truth_barcode;
      std::vector<int>     *el_truth_charge;
      Int_t           mu_muid_truth_n;
      std::vector<float>   *mu_muid_truth_E;
      std::vector<float>   *mu_muid_truth_pt;
      std::vector<float>   *mu_muid_truth_m;
      std::vector<float>   *mu_muid_truth_px;
      std::vector<float>   *mu_muid_truth_py;
      std::vector<float>   *mu_muid_truth_pz;
      std::vector<float>   *mu_muid_truth_eta;
      std::vector<float>   *mu_muid_truth_phi;
      std::vector<int>     *mu_muid_truth_status;
      std::vector<int>     *mu_muid_truth_barcode;
      std::vector<int>     *mu_muid_truth_charge;
      Int_t           mu_staco_truth_n;
      std::vector<float>   *mu_staco_truth_E;
      std::vector<float>   *mu_staco_truth_pt;
      std::vector<float>   *mu_staco_truth_m;
      std::vector<float>   *mu_staco_truth_px;
      std::vector<float>   *mu_staco_truth_py;
      std::vector<float>   *mu_staco_truth_pz;
      std::vector<float>   *mu_staco_truth_eta;
      std::vector<float>   *mu_staco_truth_phi;
      std::vector<int>     *mu_staco_truth_status;
      std::vector<int>     *mu_staco_truth_barcode;
      std::vector<int>     *mu_staco_truth_charge;
      Int_t           trueTau_n;
      std::vector<float>   *trueTau_pt;
      std::vector<float>   *trueTau_m;
      std::vector<float>   *trueTau_eta;
      std::vector<float>   *trueTau_phi;
      std::vector<int>     *trueTau_status;
      std::vector<int>     *trueTau_barcode;
      std::vector<int>     *trueTau_charge;
      Int_t           ph_truth_n;
      std::vector<float>   *ph_truth_E;
      std::vector<float>   *ph_truth_pt;
      std::vector<float>   *ph_truth_m;
      std::vector<float>   *ph_truth_px;
      std::vector<float>   *ph_truth_py;
      std::vector<float>   *ph_truth_pz;
      std::vector<float>   *ph_truth_eta;
      std::vector<float>   *ph_truth_phi;
      std::vector<int>     *ph_truth_status;
      std::vector<int>     *ph_truth_barcode;

      // List of branches
      TBranch        *b_RunNumber;   //!
      TBranch        *b_EventNumber;   //!
      TBranch        *b_timestamp;   //!
      TBranch        *b_timestamp_ns;   //!
      TBranch        *b_lbn;   //!
      TBranch        *b_bcid;   //!
      TBranch        *b_detmask0;   //!
      TBranch        *b_detmask1;   //!
      TBranch        *b_actualIntPerXing;   //!
      TBranch        *b_averageIntPerXing;   //!
      TBranch        *b_mc_channel_number;   //!
      TBranch        *b_mc_event_number;   //!
      TBranch        *b_mc_event_weight;   //!
      TBranch        *b_pixelFlags;   //!
      TBranch        *b_sctFlags;   //!
      TBranch        *b_trtFlags;   //!
      TBranch        *b_larFlags;   //!
      TBranch        *b_tileFlags;   //!
      TBranch        *b_muonFlags;   //!
      TBranch        *b_fwdFlags;   //!
      TBranch        *b_coreFlags;   //!
      TBranch        *b_pixelError;   //!
      TBranch        *b_sctError;   //!
      TBranch        *b_trtError;   //!
      TBranch        *b_larError;   //!
      TBranch        *b_tileError;   //!
      TBranch        *b_muonError;   //!
      TBranch        *b_fwdError;   //!
      TBranch        *b_coreError;   //!
      TBranch        *b_isSimulation;   //!
      TBranch        *b_isCalibration;   //!
      TBranch        *b_isTestBeam;   //!
      TBranch        *b_mcevt_n;   //!
      TBranch        *b_mcevt_signal_process_id;   //!
      TBranch        *b_mcevt_event_number;   //!
      TBranch        *b_mcevt_event_scale;   //!
      TBranch        *b_mcevt_alphaQCD;   //!
      TBranch        *b_mcevt_alphaQED;   //!
      TBranch        *b_mcevt_pdf_id1;   //!
      TBranch        *b_mcevt_pdf_id2;   //!
      TBranch        *b_mcevt_pdf_x1;   //!
      TBranch        *b_mcevt_pdf_x2;   //!
      TBranch        *b_mcevt_pdf_scale;   //!
      TBranch        *b_mcevt_pdf1;   //!
      TBranch        *b_mcevt_pdf2;   //!
      TBranch        *b_mcevt_weight;   //!
      TBranch        *b_mcevt_nparticle;   //!
      TBranch        *b_mcevt_pileUpType;   //!
      TBranch        *b_mcVxn;   //!
      TBranch        *b_mcVxx;   //!
      TBranch        *b_mcVxy;   //!
      TBranch        *b_mcVxz;   //!
      TBranch        *b_vxn;   //!
      TBranch        *b_vxx;   //!
      TBranch        *b_vxy;   //!
      TBranch        *b_vxz;   //!
      TBranch        *b_jet_AntiKt4TruthJets_n;   //!
      TBranch        *b_jet_AntiKt4TruthJets_E;   //!
      TBranch        *b_jet_AntiKt4TruthJets_pt;   //!
      TBranch        *b_jet_AntiKt4TruthJets_m;   //!
      TBranch        *b_jet_AntiKt4TruthJets_eta;   //!
      TBranch        *b_jet_AntiKt4TruthJets_phi;   //!
      TBranch        *b_jet_AntiKt4TruthJets_EtaOrigin;   //!
      TBranch        *b_jet_AntiKt4TruthJets_PhiOrigin;   //!
      TBranch        *b_jet_AntiKt4TruthJets_MOrigin;   //!
      TBranch        *b_jet_AntiKt4TruthJets_EtaOriginEM;   //!
      TBranch        *b_jet_AntiKt4TruthJets_PhiOriginEM;   //!
      TBranch        *b_jet_AntiKt4TruthJets_MOriginEM;   //!
      TBranch        *b_jet_AntiKt4TruthJets_WIDTH;   //!
      TBranch        *b_jet_AntiKt4TruthJets_n90;   //!
      TBranch        *b_jet_AntiKt4TruthJets_Timing;   //!
      TBranch        *b_jet_AntiKt4TruthJets_LArQuality;   //!
      TBranch        *b_jet_AntiKt4TruthJets_nTrk;   //!
      TBranch        *b_jet_AntiKt4TruthJets_sumPtTrk;   //!
      TBranch        *b_jet_AntiKt4TruthJets_OriginIndex;   //!
      TBranch        *b_jet_AntiKt4TruthJets_HECQuality;   //!
      TBranch        *b_jet_AntiKt4TruthJets_NegativeE;   //!
      TBranch        *b_jet_AntiKt4TruthJets_AverageLArQF;   //!
      TBranch        *b_jet_AntiKt4TruthJets_YFlip12;   //!
      TBranch        *b_jet_AntiKt4TruthJets_YFlip23;   //!
      TBranch        *b_jet_AntiKt4TruthJets_BCH_CORR_CELL;   //!
      TBranch        *b_jet_AntiKt4TruthJets_BCH_CORR_DOTX;   //!
      TBranch        *b_jet_AntiKt4TruthJets_BCH_CORR_JET;   //!
      TBranch        *b_jet_AntiKt4TruthJets_BCH_CORR_JET_FORCELL;   //!
      TBranch        *b_jet_AntiKt4TruthJets_ENG_BAD_CELLS;   //!
      TBranch        *b_jet_AntiKt4TruthJets_N_BAD_CELLS;   //!
      TBranch        *b_jet_AntiKt4TruthJets_N_BAD_CELLS_CORR;   //!
      TBranch        *b_jet_AntiKt4TruthJets_BAD_CELLS_CORR_E;   //!
      TBranch        *b_jet_AntiKt4TruthJets_NumTowers;   //!
      TBranch        *b_jet_AntiKt4TruthJets_ootFracCells5;   //!
      TBranch        *b_jet_AntiKt4TruthJets_ootFracCells10;   //!
      TBranch        *b_jet_AntiKt4TruthJets_ootFracClusters5;   //!
      TBranch        *b_jet_AntiKt4TruthJets_ootFracClusters10;   //!
      TBranch        *b_jet_AntiKt4TruthJets_SamplingMax;   //!
      TBranch        *b_jet_AntiKt4TruthJets_fracSamplingMax;   //!
      TBranch        *b_jet_AntiKt4TruthJets_hecf;   //!
      TBranch        *b_jet_AntiKt4TruthJets_tgap3f;   //!
      TBranch        *b_jet_AntiKt4TruthJets_isUgly;   //!
      TBranch        *b_jet_AntiKt4TruthJets_isBadLooseMinus;   //!
      TBranch        *b_jet_AntiKt4TruthJets_isBadLoose;   //!
      TBranch        *b_jet_AntiKt4TruthJets_isBadMedium;   //!
      TBranch        *b_jet_AntiKt4TruthJets_isBadTight;   //!
      TBranch        *b_jet_AntiKt4TruthJets_emfrac;   //!
      TBranch        *b_jet_AntiKt4TruthJets_Offset;   //!
      TBranch        *b_jet_AntiKt4TruthJets_EMJES;   //!
      TBranch        *b_jet_AntiKt4TruthJets_EMJES_EtaCorr;   //!
      TBranch        *b_jet_AntiKt4TruthJets_EMJESnooffset;   //!
      TBranch        *b_jet_AntiKt4TruthJets_GCWJES;   //!
      TBranch        *b_jet_AntiKt4TruthJets_GCWJES_EtaCorr;   //!
      TBranch        *b_jet_AntiKt4TruthJets_CB;   //!
      TBranch        *b_jet_AntiKt4TruthJets_LCJES;   //!
      TBranch        *b_jet_AntiKt4TruthJets_emscale_E;   //!
      TBranch        *b_jet_AntiKt4TruthJets_emscale_pt;   //!
      TBranch        *b_jet_AntiKt4TruthJets_emscale_m;   //!
      TBranch        *b_jet_AntiKt4TruthJets_emscale_eta;   //!
      TBranch        *b_jet_AntiKt4TruthJets_emscale_phi;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_n;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_E;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_pt;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_m;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_eta;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_phi;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_EtaOrigin;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_PhiOrigin;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_MOrigin;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_EtaOriginEM;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_PhiOriginEM;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_MOriginEM;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_WIDTH;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_n90;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_Timing;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_LArQuality;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_nTrk;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_sumPtTrk;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_OriginIndex;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_HECQuality;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_NegativeE;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_AverageLArQF;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_YFlip12;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_YFlip23;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_BCH_CORR_CELL;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_BCH_CORR_DOTX;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_BCH_CORR_JET;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_BCH_CORR_JET_FORCELL;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_ENG_BAD_CELLS;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_N_BAD_CELLS;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_N_BAD_CELLS_CORR;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_BAD_CELLS_CORR_E;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_NumTowers;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_ootFracCells5;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_ootFracCells10;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_ootFracClusters5;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_ootFracClusters10;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_SamplingMax;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_fracSamplingMax;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_hecf;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_tgap3f;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_isUgly;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_isBadLooseMinus;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_isBadLoose;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_isBadMedium;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_isBadTight;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_emfrac;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_Offset;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_EMJES;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_EMJES_EtaCorr;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_EMJESnooffset;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_GCWJES;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_GCWJES_EtaCorr;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_CB;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_LCJES;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_emscale_E;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_emscale_pt;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_emscale_m;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_emscale_eta;   //!
      TBranch        *b_jet_AntiKt4TopoNewEM_emscale_phi;   //!
      TBranch        *b_mc_n;   //!
      TBranch        *b_mc_E;   //!
      TBranch        *b_mc_pt;   //!
      TBranch        *b_mc_m;   //!
      TBranch        *b_mc_eta;   //!
      TBranch        *b_mc_phi;   //!
      TBranch        *b_mc_px;   //!
      TBranch        *b_mc_py;   //!
      TBranch        *b_mc_pz;   //!
      TBranch        *b_mc_status;   //!
      TBranch        *b_mc_barcode;   //!
      TBranch        *b_mc_pdgId;   //!
      TBranch        *b_mc_charge;   //!
      TBranch        *b_mc_vx_x;   //!
      TBranch        *b_mc_vx_y;   //!
      TBranch        *b_mc_vx_z;   //!
      TBranch        *b_mc_vx_barcode;   //!
      TBranch        *b_mc_child_index;   //!
      TBranch        *b_mc_parent_index;   //!
      TBranch        *b_MET_Truth_NonInt_etx;   //!
      TBranch        *b_MET_Truth_NonInt_ety;   //!
      TBranch        *b_MET_Truth_Int_etx;   //!
      TBranch        *b_MET_Truth_Int_ety;   //!
      TBranch        *b_MET_Truth_IntCentral_etx;   //!
      TBranch        *b_MET_Truth_IntCentral_ety;   //!
      TBranch        *b_MET_Truth_IntFwd_etx;   //!
      TBranch        *b_MET_Truth_IntFwd_ety;   //!
      TBranch        *b_MET_Truth_IntOutCover_etx;   //!
      TBranch        *b_MET_Truth_IntOutCover_ety;   //!
      TBranch        *b_MET_Truth_IntMuons_etx;   //!
      TBranch        *b_MET_Truth_IntMuons_ety;   //!
      TBranch        *b_el_n;   //!
      TBranch        *b_el_E;   //!
      TBranch        *b_el_pt;   //!
      TBranch        *b_el_m;   //!
      TBranch        *b_el_px;   //!
      TBranch        *b_el_py;   //!
      TBranch        *b_el_pz;   //!
      TBranch        *b_el_eta;   //!
      TBranch        *b_el_phi;   //!
      TBranch        *b_el_status;   //!
      TBranch        *b_el_barcode;   //!
      TBranch        *b_el_charge;   //!
      TBranch        *b_mu_muid_n;   //!
      TBranch        *b_mu_muid_E;   //!
      TBranch        *b_mu_muid_pt;   //!
      TBranch        *b_mu_muid_m;   //!
      TBranch        *b_mu_muid_px;   //!
      TBranch        *b_mu_muid_py;   //!
      TBranch        *b_mu_muid_pz;   //!
      TBranch        *b_mu_muid_eta;   //!
      TBranch        *b_mu_muid_phi;   //!
      TBranch        *b_mu_muid_status;   //!
      TBranch        *b_mu_muid_barcode;   //!
      TBranch        *b_mu_muid_charge;   //!
      TBranch        *b_mu_staco_n;   //!
      TBranch        *b_mu_staco_E;   //!
      TBranch        *b_mu_staco_pt;   //!
      TBranch        *b_mu_staco_m;   //!
      TBranch        *b_mu_staco_px;   //!
      TBranch        *b_mu_staco_py;   //!
      TBranch        *b_mu_staco_pz;   //!
      TBranch        *b_mu_staco_eta;   //!
      TBranch        *b_mu_staco_phi;   //!
      TBranch        *b_mu_staco_status;   //!
      TBranch        *b_mu_staco_barcode;   //!
      TBranch        *b_mu_staco_charge;   //!
      TBranch        *b_tau_n;   //!
      TBranch        *b_tau_pt;   //!
      TBranch        *b_tau_m;   //!
      TBranch        *b_tau_eta;   //!
      TBranch        *b_tau_phi;   //!
      TBranch        *b_tau_status;   //!
      TBranch        *b_tau_barcode;   //!
      TBranch        *b_tau_charge;   //!
      TBranch        *b_ph_n;   //!
      TBranch        *b_ph_E;   //!
      TBranch        *b_ph_pt;   //!
      TBranch        *b_ph_m;   //!
      TBranch        *b_ph_px;   //!
      TBranch        *b_ph_py;   //!
      TBranch        *b_ph_pz;   //!
      TBranch        *b_ph_eta;   //!
      TBranch        *b_ph_phi;   //!
      TBranch        *b_ph_status;   //!
      TBranch        *b_ph_barcode;   //!
      TBranch        *b_el_truth_n;   //!
      TBranch        *b_el_truth_E;   //!
      TBranch        *b_el_truth_pt;   //!
      TBranch        *b_el_truth_m;   //!
      TBranch        *b_el_truth_px;   //!
      TBranch        *b_el_truth_py;   //!
      TBranch        *b_el_truth_pz;   //!
      TBranch        *b_el_truth_eta;   //!
      TBranch        *b_el_truth_phi;   //!
      TBranch        *b_el_truth_status;   //!
      TBranch        *b_el_truth_barcode;   //!
      TBranch        *b_el_truth_charge;   //!
      TBranch        *b_mu_muid_truth_n;   //!
      TBranch        *b_mu_muid_truth_E;   //!
      TBranch        *b_mu_muid_truth_pt;   //!
      TBranch        *b_mu_muid_truth_m;   //!
      TBranch        *b_mu_muid_truth_px;   //!
      TBranch        *b_mu_muid_truth_py;   //!
      TBranch        *b_mu_muid_truth_pz;   //!
      TBranch        *b_mu_muid_truth_eta;   //!
      TBranch        *b_mu_muid_truth_phi;   //!
      TBranch        *b_mu_muid_truth_status;   //!
      TBranch        *b_mu_muid_truth_barcode;   //!
      TBranch        *b_mu_muid_truth_charge;   //!
      TBranch        *b_mu_staco_truth_n;   //!
      TBranch        *b_mu_staco_truth_E;   //!
      TBranch        *b_mu_staco_truth_pt;   //!
      TBranch        *b_mu_staco_truth_m;   //!
      TBranch        *b_mu_staco_truth_px;   //!
      TBranch        *b_mu_staco_truth_py;   //!
      TBranch        *b_mu_staco_truth_pz;   //!
      TBranch        *b_mu_staco_truth_eta;   //!
      TBranch        *b_mu_staco_truth_phi;   //!
      TBranch        *b_mu_staco_truth_status;   //!
      TBranch        *b_mu_staco_truth_barcode;   //!
      TBranch        *b_mu_staco_truth_charge;   //!
      TBranch        *b_trueTau_n;   //!
      TBranch        *b_trueTau_pt;   //!
      TBranch        *b_trueTau_m;   //!
      TBranch        *b_trueTau_eta;   //!
      TBranch        *b_trueTau_phi;   //!
      TBranch        *b_trueTau_status;   //!
      TBranch        *b_trueTau_barcode;   //!
      TBranch        *b_trueTau_charge;   //!
      TBranch        *b_ph_truth_n;   //!
      TBranch        *b_ph_truth_E;   //!
      TBranch        *b_ph_truth_pt;   //!
      TBranch        *b_ph_truth_m;   //!
      TBranch        *b_ph_truth_px;   //!
      TBranch        *b_ph_truth_py;   //!
      TBranch        *b_ph_truth_pz;   //!
      TBranch        *b_ph_truth_eta;   //!
      TBranch        *b_ph_truth_phi;   //!
      TBranch        *b_ph_truth_status;   //!
      TBranch        *b_ph_truth_barcode;   //!

    private:

      friend class Particle;
      friend class Electron;
      friend class Muon;
      friend class Jet;
  };
}

#endif
