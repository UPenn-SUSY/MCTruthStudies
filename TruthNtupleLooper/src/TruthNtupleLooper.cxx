#include "include/TruthNtupleLooper.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <dirent.h>
#include <math.h>
#include <algorithm>

// #include "TH2.h"

#include "include/ObjectDefs.h"
#include "ProgressBar/include/ProgressBar.h"

// -----------------------------------------------------------------------------
TruthNtuple::TruthNtupleLooper::TruthNtupleLooper(TTree *tree) : fChain(0)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../evgen.SampleProcess.TRUTH.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("../../evgen.SampleProcess.TRUTH.root");
    }
    f->GetObject("truth",tree);
  }
  Init(tree);
}

// -----------------------------------------------------------------------------
TruthNtuple::TruthNtupleLooper::~TruthNtupleLooper()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

// -----------------------------------------------------------------------------
Int_t TruthNtuple::TruthNtupleLooper::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

// -----------------------------------------------------------------------------
Long64_t TruthNtuple::TruthNtupleLooper::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

// -----------------------------------------------------------------------------
void TruthNtuple::TruthNtupleLooper::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  mcevt_signal_process_id = 0;
  mcevt_event_number = 0;
  mcevt_event_scale = 0;
  mcevt_alphaQCD = 0;
  mcevt_alphaQED = 0;
  mcevt_pdf_id1 = 0;
  mcevt_pdf_id2 = 0;
  mcevt_pdf_x1 = 0;
  mcevt_pdf_x2 = 0;
  mcevt_pdf_scale = 0;
  mcevt_pdf1 = 0;
  mcevt_pdf2 = 0;
  mcevt_weight = 0;
  mcevt_nparticle = 0;
  mcevt_pileUpType = 0;
  mcVxx = 0;
  mcVxy = 0;
  mcVxz = 0;
  vxx = 0;
  vxy = 0;
  vxz = 0;
  jet_AntiKt4TruthJets_E = 0;
  jet_AntiKt4TruthJets_pt = 0;
  jet_AntiKt4TruthJets_m = 0;
  jet_AntiKt4TruthJets_eta = 0;
  jet_AntiKt4TruthJets_phi = 0;
  jet_AntiKt4TruthJets_EtaOrigin = 0;
  jet_AntiKt4TruthJets_PhiOrigin = 0;
  jet_AntiKt4TruthJets_MOrigin = 0;
  jet_AntiKt4TruthJets_EtaOriginEM = 0;
  jet_AntiKt4TruthJets_PhiOriginEM = 0;
  jet_AntiKt4TruthJets_MOriginEM = 0;
  jet_AntiKt4TruthJets_WIDTH = 0;
  jet_AntiKt4TruthJets_n90 = 0;
  jet_AntiKt4TruthJets_Timing = 0;
  jet_AntiKt4TruthJets_LArQuality = 0;
  jet_AntiKt4TruthJets_nTrk = 0;
  jet_AntiKt4TruthJets_sumPtTrk = 0;
  jet_AntiKt4TruthJets_OriginIndex = 0;
  jet_AntiKt4TruthJets_HECQuality = 0;
  jet_AntiKt4TruthJets_NegativeE = 0;
  jet_AntiKt4TruthJets_AverageLArQF = 0;
  jet_AntiKt4TruthJets_YFlip12 = 0;
  jet_AntiKt4TruthJets_YFlip23 = 0;
  jet_AntiKt4TruthJets_BCH_CORR_CELL = 0;
  jet_AntiKt4TruthJets_BCH_CORR_DOTX = 0;
  jet_AntiKt4TruthJets_BCH_CORR_JET = 0;
  jet_AntiKt4TruthJets_BCH_CORR_JET_FORCELL = 0;
  jet_AntiKt4TruthJets_ENG_BAD_CELLS = 0;
  jet_AntiKt4TruthJets_N_BAD_CELLS = 0;
  jet_AntiKt4TruthJets_N_BAD_CELLS_CORR = 0;
  jet_AntiKt4TruthJets_BAD_CELLS_CORR_E = 0;
  jet_AntiKt4TruthJets_NumTowers = 0;
  jet_AntiKt4TruthJets_ootFracCells5 = 0;
  jet_AntiKt4TruthJets_ootFracCells10 = 0;
  jet_AntiKt4TruthJets_ootFracClusters5 = 0;
  jet_AntiKt4TruthJets_ootFracClusters10 = 0;
  jet_AntiKt4TruthJets_SamplingMax = 0;
  jet_AntiKt4TruthJets_fracSamplingMax = 0;
  jet_AntiKt4TruthJets_hecf = 0;
  jet_AntiKt4TruthJets_tgap3f = 0;
  jet_AntiKt4TruthJets_isUgly = 0;
  jet_AntiKt4TruthJets_isBadLooseMinus = 0;
  jet_AntiKt4TruthJets_isBadLoose = 0;
  jet_AntiKt4TruthJets_isBadMedium = 0;
  jet_AntiKt4TruthJets_isBadTight = 0;
  jet_AntiKt4TruthJets_emfrac = 0;
  jet_AntiKt4TruthJets_Offset = 0;
  jet_AntiKt4TruthJets_EMJES = 0;
  jet_AntiKt4TruthJets_EMJES_EtaCorr = 0;
  jet_AntiKt4TruthJets_EMJESnooffset = 0;
  jet_AntiKt4TruthJets_GCWJES = 0;
  jet_AntiKt4TruthJets_GCWJES_EtaCorr = 0;
  jet_AntiKt4TruthJets_CB = 0;
  jet_AntiKt4TruthJets_LCJES = 0;
  jet_AntiKt4TruthJets_emscale_E = 0;
  jet_AntiKt4TruthJets_emscale_pt = 0;
  jet_AntiKt4TruthJets_emscale_m = 0;
  jet_AntiKt4TruthJets_emscale_eta = 0;
  jet_AntiKt4TruthJets_emscale_phi = 0;
  jet_AntiKt4TopoNewEM_E = 0;
  jet_AntiKt4TopoNewEM_pt = 0;
  jet_AntiKt4TopoNewEM_m = 0;
  jet_AntiKt4TopoNewEM_eta = 0;
  jet_AntiKt4TopoNewEM_phi = 0;
  jet_AntiKt4TopoNewEM_EtaOrigin = 0;
  jet_AntiKt4TopoNewEM_PhiOrigin = 0;
  jet_AntiKt4TopoNewEM_MOrigin = 0;
  jet_AntiKt4TopoNewEM_EtaOriginEM = 0;
  jet_AntiKt4TopoNewEM_PhiOriginEM = 0;
  jet_AntiKt4TopoNewEM_MOriginEM = 0;
  jet_AntiKt4TopoNewEM_WIDTH = 0;
  jet_AntiKt4TopoNewEM_n90 = 0;
  jet_AntiKt4TopoNewEM_Timing = 0;
  jet_AntiKt4TopoNewEM_LArQuality = 0;
  jet_AntiKt4TopoNewEM_nTrk = 0;
  jet_AntiKt4TopoNewEM_sumPtTrk = 0;
  jet_AntiKt4TopoNewEM_OriginIndex = 0;
  jet_AntiKt4TopoNewEM_HECQuality = 0;
  jet_AntiKt4TopoNewEM_NegativeE = 0;
  jet_AntiKt4TopoNewEM_AverageLArQF = 0;
  jet_AntiKt4TopoNewEM_YFlip12 = 0;
  jet_AntiKt4TopoNewEM_YFlip23 = 0;
  jet_AntiKt4TopoNewEM_BCH_CORR_CELL = 0;
  jet_AntiKt4TopoNewEM_BCH_CORR_DOTX = 0;
  jet_AntiKt4TopoNewEM_BCH_CORR_JET = 0;
  jet_AntiKt4TopoNewEM_BCH_CORR_JET_FORCELL = 0;
  jet_AntiKt4TopoNewEM_ENG_BAD_CELLS = 0;
  jet_AntiKt4TopoNewEM_N_BAD_CELLS = 0;
  jet_AntiKt4TopoNewEM_N_BAD_CELLS_CORR = 0;
  jet_AntiKt4TopoNewEM_BAD_CELLS_CORR_E = 0;
  jet_AntiKt4TopoNewEM_NumTowers = 0;
  jet_AntiKt4TopoNewEM_ootFracCells5 = 0;
  jet_AntiKt4TopoNewEM_ootFracCells10 = 0;
  jet_AntiKt4TopoNewEM_ootFracClusters5 = 0;
  jet_AntiKt4TopoNewEM_ootFracClusters10 = 0;
  jet_AntiKt4TopoNewEM_SamplingMax = 0;
  jet_AntiKt4TopoNewEM_fracSamplingMax = 0;
  jet_AntiKt4TopoNewEM_hecf = 0;
  jet_AntiKt4TopoNewEM_tgap3f = 0;
  jet_AntiKt4TopoNewEM_isUgly = 0;
  jet_AntiKt4TopoNewEM_isBadLooseMinus = 0;
  jet_AntiKt4TopoNewEM_isBadLoose = 0;
  jet_AntiKt4TopoNewEM_isBadMedium = 0;
  jet_AntiKt4TopoNewEM_isBadTight = 0;
  jet_AntiKt4TopoNewEM_emfrac = 0;
  jet_AntiKt4TopoNewEM_Offset = 0;
  jet_AntiKt4TopoNewEM_EMJES = 0;
  jet_AntiKt4TopoNewEM_EMJES_EtaCorr = 0;
  jet_AntiKt4TopoNewEM_EMJESnooffset = 0;
  jet_AntiKt4TopoNewEM_GCWJES = 0;
  jet_AntiKt4TopoNewEM_GCWJES_EtaCorr = 0;
  jet_AntiKt4TopoNewEM_CB = 0;
  jet_AntiKt4TopoNewEM_LCJES = 0;
  jet_AntiKt4TopoNewEM_emscale_E = 0;
  jet_AntiKt4TopoNewEM_emscale_pt = 0;
  jet_AntiKt4TopoNewEM_emscale_m = 0;
  jet_AntiKt4TopoNewEM_emscale_eta = 0;
  jet_AntiKt4TopoNewEM_emscale_phi = 0;
  mc_E = 0;
  mc_pt = 0;
  mc_m = 0;
  mc_eta = 0;
  mc_phi = 0;
  mc_px = 0;
  mc_py = 0;
  mc_pz = 0;
  mc_status = 0;
  mc_barcode = 0;
  mc_pdgId = 0;
  mc_charge = 0;
  mc_vx_x = 0;
  mc_vx_y = 0;
  mc_vx_z = 0;
  mc_vx_barcode = 0;
  mc_child_index = 0;
  mc_parent_index = 0;
  el_E = 0;
  el_pt = 0;
  el_m = 0;
  el_px = 0;
  el_py = 0;
  el_pz = 0;
  el_eta = 0;
  el_phi = 0;
  el_status = 0;
  el_barcode = 0;
  el_charge = 0;
  mu_muid_E = 0;
  mu_muid_pt = 0;
  mu_muid_m = 0;
  mu_muid_px = 0;
  mu_muid_py = 0;
  mu_muid_pz = 0;
  mu_muid_eta = 0;
  mu_muid_phi = 0;
  mu_muid_status = 0;
  mu_muid_barcode = 0;
  mu_muid_charge = 0;
  mu_staco_E = 0;
  mu_staco_pt = 0;
  mu_staco_m = 0;
  mu_staco_px = 0;
  mu_staco_py = 0;
  mu_staco_pz = 0;
  mu_staco_eta = 0;
  mu_staco_phi = 0;
  mu_staco_status = 0;
  mu_staco_barcode = 0;
  mu_staco_charge = 0;
  tau_pt = 0;
  tau_m = 0;
  tau_eta = 0;
  tau_phi = 0;
  tau_status = 0;
  tau_barcode = 0;
  tau_charge = 0;
  ph_E = 0;
  ph_pt = 0;
  ph_m = 0;
  ph_px = 0;
  ph_py = 0;
  ph_pz = 0;
  ph_eta = 0;
  ph_phi = 0;
  ph_status = 0;
  ph_barcode = 0;
  el_truth_E = 0;
  el_truth_pt = 0;
  el_truth_m = 0;
  el_truth_px = 0;
  el_truth_py = 0;
  el_truth_pz = 0;
  el_truth_eta = 0;
  el_truth_phi = 0;
  el_truth_status = 0;
  el_truth_barcode = 0;
  el_truth_charge = 0;
  mu_muid_truth_E = 0;
  mu_muid_truth_pt = 0;
  mu_muid_truth_m = 0;
  mu_muid_truth_px = 0;
  mu_muid_truth_py = 0;
  mu_muid_truth_pz = 0;
  mu_muid_truth_eta = 0;
  mu_muid_truth_phi = 0;
  mu_muid_truth_status = 0;
  mu_muid_truth_barcode = 0;
  mu_muid_truth_charge = 0;
  mu_staco_truth_E = 0;
  mu_staco_truth_pt = 0;
  mu_staco_truth_m = 0;
  mu_staco_truth_px = 0;
  mu_staco_truth_py = 0;
  mu_staco_truth_pz = 0;
  mu_staco_truth_eta = 0;
  mu_staco_truth_phi = 0;
  mu_staco_truth_status = 0;
  mu_staco_truth_barcode = 0;
  mu_staco_truth_charge = 0;
  trueTau_pt = 0;
  trueTau_m = 0;
  trueTau_eta = 0;
  trueTau_phi = 0;
  trueTau_status = 0;
  trueTau_barcode = 0;
  trueTau_charge = 0;
  ph_truth_E = 0;
  ph_truth_pt = 0;
  ph_truth_m = 0;
  ph_truth_px = 0;
  ph_truth_py = 0;
  ph_truth_pz = 0;
  ph_truth_eta = 0;
  ph_truth_phi = 0;
  ph_truth_status = 0;
  ph_truth_barcode = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
  fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
  fChain->SetBranchAddress("timestamp", &timestamp, &b_timestamp);
  fChain->SetBranchAddress("timestamp_ns", &timestamp_ns, &b_timestamp_ns);
  fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
  fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
  fChain->SetBranchAddress("detmask0", &detmask0, &b_detmask0);
  fChain->SetBranchAddress("detmask1", &detmask1, &b_detmask1);
  fChain->SetBranchAddress("actualIntPerXing", &actualIntPerXing, &b_actualIntPerXing);
  fChain->SetBranchAddress("averageIntPerXing", &averageIntPerXing, &b_averageIntPerXing);
  fChain->SetBranchAddress("mc_channel_number", &mc_channel_number, &b_mc_channel_number);
  fChain->SetBranchAddress("mc_event_number", &mc_event_number, &b_mc_event_number);
  fChain->SetBranchAddress("mc_event_weight", &mc_event_weight, &b_mc_event_weight);
  fChain->SetBranchAddress("pixelFlags", &pixelFlags, &b_pixelFlags);
  fChain->SetBranchAddress("sctFlags", &sctFlags, &b_sctFlags);
  fChain->SetBranchAddress("trtFlags", &trtFlags, &b_trtFlags);
  fChain->SetBranchAddress("larFlags", &larFlags, &b_larFlags);
  fChain->SetBranchAddress("tileFlags", &tileFlags, &b_tileFlags);
  fChain->SetBranchAddress("muonFlags", &muonFlags, &b_muonFlags);
  fChain->SetBranchAddress("fwdFlags", &fwdFlags, &b_fwdFlags);
  fChain->SetBranchAddress("coreFlags", &coreFlags, &b_coreFlags);
  fChain->SetBranchAddress("pixelError", &pixelError, &b_pixelError);
  fChain->SetBranchAddress("sctError", &sctError, &b_sctError);
  fChain->SetBranchAddress("trtError", &trtError, &b_trtError);
  fChain->SetBranchAddress("larError", &larError, &b_larError);
  fChain->SetBranchAddress("tileError", &tileError, &b_tileError);
  fChain->SetBranchAddress("muonError", &muonError, &b_muonError);
  fChain->SetBranchAddress("fwdError", &fwdError, &b_fwdError);
  fChain->SetBranchAddress("coreError", &coreError, &b_coreError);
  fChain->SetBranchAddress("isSimulation", &isSimulation, &b_isSimulation);
  fChain->SetBranchAddress("isCalibration", &isCalibration, &b_isCalibration);
  fChain->SetBranchAddress("isTestBeam", &isTestBeam, &b_isTestBeam);
  fChain->SetBranchAddress("mcevt_n", &mcevt_n, &b_mcevt_n);
  fChain->SetBranchAddress("mcevt_signal_process_id", &mcevt_signal_process_id, &b_mcevt_signal_process_id);
  fChain->SetBranchAddress("mcevt_event_number", &mcevt_event_number, &b_mcevt_event_number);
  fChain->SetBranchAddress("mcevt_event_scale", &mcevt_event_scale, &b_mcevt_event_scale);
  fChain->SetBranchAddress("mcevt_alphaQCD", &mcevt_alphaQCD, &b_mcevt_alphaQCD);
  fChain->SetBranchAddress("mcevt_alphaQED", &mcevt_alphaQED, &b_mcevt_alphaQED);
  fChain->SetBranchAddress("mcevt_pdf_id1", &mcevt_pdf_id1, &b_mcevt_pdf_id1);
  fChain->SetBranchAddress("mcevt_pdf_id2", &mcevt_pdf_id2, &b_mcevt_pdf_id2);
  fChain->SetBranchAddress("mcevt_pdf_x1", &mcevt_pdf_x1, &b_mcevt_pdf_x1);
  fChain->SetBranchAddress("mcevt_pdf_x2", &mcevt_pdf_x2, &b_mcevt_pdf_x2);
  fChain->SetBranchAddress("mcevt_pdf_scale", &mcevt_pdf_scale, &b_mcevt_pdf_scale);
  fChain->SetBranchAddress("mcevt_pdf1", &mcevt_pdf1, &b_mcevt_pdf1);
  fChain->SetBranchAddress("mcevt_pdf2", &mcevt_pdf2, &b_mcevt_pdf2);
  fChain->SetBranchAddress("mcevt_weight", &mcevt_weight, &b_mcevt_weight);
  fChain->SetBranchAddress("mcevt_nparticle", &mcevt_nparticle, &b_mcevt_nparticle);
  fChain->SetBranchAddress("mcevt_pileUpType", &mcevt_pileUpType, &b_mcevt_pileUpType);
  fChain->SetBranchAddress("mcVxn", &mcVxn, &b_mcVxn);
  fChain->SetBranchAddress("mcVxx", &mcVxx, &b_mcVxx);
  fChain->SetBranchAddress("mcVxy", &mcVxy, &b_mcVxy);
  fChain->SetBranchAddress("mcVxz", &mcVxz, &b_mcVxz);
  fChain->SetBranchAddress("vxn", &vxn, &b_vxn);
  fChain->SetBranchAddress("vxx", &vxx, &b_vxx);
  fChain->SetBranchAddress("vxy", &vxy, &b_vxy);
  fChain->SetBranchAddress("vxz", &vxz, &b_vxz);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_n", &jet_AntiKt4TruthJets_n, &b_jet_AntiKt4TruthJets_n);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_E", &jet_AntiKt4TruthJets_E, &b_jet_AntiKt4TruthJets_E);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_pt", &jet_AntiKt4TruthJets_pt, &b_jet_AntiKt4TruthJets_pt);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_m", &jet_AntiKt4TruthJets_m, &b_jet_AntiKt4TruthJets_m);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_eta", &jet_AntiKt4TruthJets_eta, &b_jet_AntiKt4TruthJets_eta);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_phi", &jet_AntiKt4TruthJets_phi, &b_jet_AntiKt4TruthJets_phi);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_EtaOrigin", &jet_AntiKt4TruthJets_EtaOrigin, &b_jet_AntiKt4TruthJets_EtaOrigin);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_PhiOrigin", &jet_AntiKt4TruthJets_PhiOrigin, &b_jet_AntiKt4TruthJets_PhiOrigin);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_MOrigin", &jet_AntiKt4TruthJets_MOrigin, &b_jet_AntiKt4TruthJets_MOrigin);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_EtaOriginEM", &jet_AntiKt4TruthJets_EtaOriginEM, &b_jet_AntiKt4TruthJets_EtaOriginEM);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_PhiOriginEM", &jet_AntiKt4TruthJets_PhiOriginEM, &b_jet_AntiKt4TruthJets_PhiOriginEM);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_MOriginEM", &jet_AntiKt4TruthJets_MOriginEM, &b_jet_AntiKt4TruthJets_MOriginEM);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_WIDTH", &jet_AntiKt4TruthJets_WIDTH, &b_jet_AntiKt4TruthJets_WIDTH);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_n90", &jet_AntiKt4TruthJets_n90, &b_jet_AntiKt4TruthJets_n90);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_Timing", &jet_AntiKt4TruthJets_Timing, &b_jet_AntiKt4TruthJets_Timing);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_LArQuality", &jet_AntiKt4TruthJets_LArQuality, &b_jet_AntiKt4TruthJets_LArQuality);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_nTrk", &jet_AntiKt4TruthJets_nTrk, &b_jet_AntiKt4TruthJets_nTrk);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_sumPtTrk", &jet_AntiKt4TruthJets_sumPtTrk, &b_jet_AntiKt4TruthJets_sumPtTrk);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_OriginIndex", &jet_AntiKt4TruthJets_OriginIndex, &b_jet_AntiKt4TruthJets_OriginIndex);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_HECQuality", &jet_AntiKt4TruthJets_HECQuality, &b_jet_AntiKt4TruthJets_HECQuality);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_NegativeE", &jet_AntiKt4TruthJets_NegativeE, &b_jet_AntiKt4TruthJets_NegativeE);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_AverageLArQF", &jet_AntiKt4TruthJets_AverageLArQF, &b_jet_AntiKt4TruthJets_AverageLArQF);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_YFlip12", &jet_AntiKt4TruthJets_YFlip12, &b_jet_AntiKt4TruthJets_YFlip12);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_YFlip23", &jet_AntiKt4TruthJets_YFlip23, &b_jet_AntiKt4TruthJets_YFlip23);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_BCH_CORR_CELL", &jet_AntiKt4TruthJets_BCH_CORR_CELL, &b_jet_AntiKt4TruthJets_BCH_CORR_CELL);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_BCH_CORR_DOTX", &jet_AntiKt4TruthJets_BCH_CORR_DOTX, &b_jet_AntiKt4TruthJets_BCH_CORR_DOTX);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_BCH_CORR_JET", &jet_AntiKt4TruthJets_BCH_CORR_JET, &b_jet_AntiKt4TruthJets_BCH_CORR_JET);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_BCH_CORR_JET_FORCELL", &jet_AntiKt4TruthJets_BCH_CORR_JET_FORCELL, &b_jet_AntiKt4TruthJets_BCH_CORR_JET_FORCELL);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_ENG_BAD_CELLS", &jet_AntiKt4TruthJets_ENG_BAD_CELLS, &b_jet_AntiKt4TruthJets_ENG_BAD_CELLS);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_N_BAD_CELLS", &jet_AntiKt4TruthJets_N_BAD_CELLS, &b_jet_AntiKt4TruthJets_N_BAD_CELLS);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_N_BAD_CELLS_CORR", &jet_AntiKt4TruthJets_N_BAD_CELLS_CORR, &b_jet_AntiKt4TruthJets_N_BAD_CELLS_CORR);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_BAD_CELLS_CORR_E", &jet_AntiKt4TruthJets_BAD_CELLS_CORR_E, &b_jet_AntiKt4TruthJets_BAD_CELLS_CORR_E);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_NumTowers", &jet_AntiKt4TruthJets_NumTowers, &b_jet_AntiKt4TruthJets_NumTowers);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_ootFracCells5", &jet_AntiKt4TruthJets_ootFracCells5, &b_jet_AntiKt4TruthJets_ootFracCells5);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_ootFracCells10", &jet_AntiKt4TruthJets_ootFracCells10, &b_jet_AntiKt4TruthJets_ootFracCells10);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_ootFracClusters5", &jet_AntiKt4TruthJets_ootFracClusters5, &b_jet_AntiKt4TruthJets_ootFracClusters5);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_ootFracClusters10", &jet_AntiKt4TruthJets_ootFracClusters10, &b_jet_AntiKt4TruthJets_ootFracClusters10);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_SamplingMax", &jet_AntiKt4TruthJets_SamplingMax, &b_jet_AntiKt4TruthJets_SamplingMax);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_fracSamplingMax", &jet_AntiKt4TruthJets_fracSamplingMax, &b_jet_AntiKt4TruthJets_fracSamplingMax);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_hecf", &jet_AntiKt4TruthJets_hecf, &b_jet_AntiKt4TruthJets_hecf);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_tgap3f", &jet_AntiKt4TruthJets_tgap3f, &b_jet_AntiKt4TruthJets_tgap3f);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_isUgly", &jet_AntiKt4TruthJets_isUgly, &b_jet_AntiKt4TruthJets_isUgly);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_isBadLooseMinus", &jet_AntiKt4TruthJets_isBadLooseMinus, &b_jet_AntiKt4TruthJets_isBadLooseMinus);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_isBadLoose", &jet_AntiKt4TruthJets_isBadLoose, &b_jet_AntiKt4TruthJets_isBadLoose);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_isBadMedium", &jet_AntiKt4TruthJets_isBadMedium, &b_jet_AntiKt4TruthJets_isBadMedium);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_isBadTight", &jet_AntiKt4TruthJets_isBadTight, &b_jet_AntiKt4TruthJets_isBadTight);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_emfrac", &jet_AntiKt4TruthJets_emfrac, &b_jet_AntiKt4TruthJets_emfrac);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_Offset", &jet_AntiKt4TruthJets_Offset, &b_jet_AntiKt4TruthJets_Offset);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_EMJES", &jet_AntiKt4TruthJets_EMJES, &b_jet_AntiKt4TruthJets_EMJES);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_EMJES_EtaCorr", &jet_AntiKt4TruthJets_EMJES_EtaCorr, &b_jet_AntiKt4TruthJets_EMJES_EtaCorr);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_EMJESnooffset", &jet_AntiKt4TruthJets_EMJESnooffset, &b_jet_AntiKt4TruthJets_EMJESnooffset);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_GCWJES", &jet_AntiKt4TruthJets_GCWJES, &b_jet_AntiKt4TruthJets_GCWJES);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_GCWJES_EtaCorr", &jet_AntiKt4TruthJets_GCWJES_EtaCorr, &b_jet_AntiKt4TruthJets_GCWJES_EtaCorr);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_CB", &jet_AntiKt4TruthJets_CB, &b_jet_AntiKt4TruthJets_CB);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_LCJES", &jet_AntiKt4TruthJets_LCJES, &b_jet_AntiKt4TruthJets_LCJES);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_emscale_E", &jet_AntiKt4TruthJets_emscale_E, &b_jet_AntiKt4TruthJets_emscale_E);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_emscale_pt", &jet_AntiKt4TruthJets_emscale_pt, &b_jet_AntiKt4TruthJets_emscale_pt);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_emscale_m", &jet_AntiKt4TruthJets_emscale_m, &b_jet_AntiKt4TruthJets_emscale_m);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_emscale_eta", &jet_AntiKt4TruthJets_emscale_eta, &b_jet_AntiKt4TruthJets_emscale_eta);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_emscale_phi", &jet_AntiKt4TruthJets_emscale_phi, &b_jet_AntiKt4TruthJets_emscale_phi);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_n", &jet_AntiKt4TopoNewEM_n, &b_jet_AntiKt4TopoNewEM_n);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_E", &jet_AntiKt4TopoNewEM_E, &b_jet_AntiKt4TopoNewEM_E);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_pt", &jet_AntiKt4TopoNewEM_pt, &b_jet_AntiKt4TopoNewEM_pt);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_m", &jet_AntiKt4TopoNewEM_m, &b_jet_AntiKt4TopoNewEM_m);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_eta", &jet_AntiKt4TopoNewEM_eta, &b_jet_AntiKt4TopoNewEM_eta);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_phi", &jet_AntiKt4TopoNewEM_phi, &b_jet_AntiKt4TopoNewEM_phi);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_EtaOrigin", &jet_AntiKt4TopoNewEM_EtaOrigin, &b_jet_AntiKt4TopoNewEM_EtaOrigin);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_PhiOrigin", &jet_AntiKt4TopoNewEM_PhiOrigin, &b_jet_AntiKt4TopoNewEM_PhiOrigin);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_MOrigin", &jet_AntiKt4TopoNewEM_MOrigin, &b_jet_AntiKt4TopoNewEM_MOrigin);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_EtaOriginEM", &jet_AntiKt4TopoNewEM_EtaOriginEM, &b_jet_AntiKt4TopoNewEM_EtaOriginEM);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_PhiOriginEM", &jet_AntiKt4TopoNewEM_PhiOriginEM, &b_jet_AntiKt4TopoNewEM_PhiOriginEM);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_MOriginEM", &jet_AntiKt4TopoNewEM_MOriginEM, &b_jet_AntiKt4TopoNewEM_MOriginEM);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_WIDTH", &jet_AntiKt4TopoNewEM_WIDTH, &b_jet_AntiKt4TopoNewEM_WIDTH);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_n90", &jet_AntiKt4TopoNewEM_n90, &b_jet_AntiKt4TopoNewEM_n90);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_Timing", &jet_AntiKt4TopoNewEM_Timing, &b_jet_AntiKt4TopoNewEM_Timing);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_LArQuality", &jet_AntiKt4TopoNewEM_LArQuality, &b_jet_AntiKt4TopoNewEM_LArQuality);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_nTrk", &jet_AntiKt4TopoNewEM_nTrk, &b_jet_AntiKt4TopoNewEM_nTrk);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_sumPtTrk", &jet_AntiKt4TopoNewEM_sumPtTrk, &b_jet_AntiKt4TopoNewEM_sumPtTrk);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_OriginIndex", &jet_AntiKt4TopoNewEM_OriginIndex, &b_jet_AntiKt4TopoNewEM_OriginIndex);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_HECQuality", &jet_AntiKt4TopoNewEM_HECQuality, &b_jet_AntiKt4TopoNewEM_HECQuality);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_NegativeE", &jet_AntiKt4TopoNewEM_NegativeE, &b_jet_AntiKt4TopoNewEM_NegativeE);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_AverageLArQF", &jet_AntiKt4TopoNewEM_AverageLArQF, &b_jet_AntiKt4TopoNewEM_AverageLArQF);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_YFlip12", &jet_AntiKt4TopoNewEM_YFlip12, &b_jet_AntiKt4TopoNewEM_YFlip12);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_YFlip23", &jet_AntiKt4TopoNewEM_YFlip23, &b_jet_AntiKt4TopoNewEM_YFlip23);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_BCH_CORR_CELL", &jet_AntiKt4TopoNewEM_BCH_CORR_CELL, &b_jet_AntiKt4TopoNewEM_BCH_CORR_CELL);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_BCH_CORR_DOTX", &jet_AntiKt4TopoNewEM_BCH_CORR_DOTX, &b_jet_AntiKt4TopoNewEM_BCH_CORR_DOTX);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_BCH_CORR_JET", &jet_AntiKt4TopoNewEM_BCH_CORR_JET, &b_jet_AntiKt4TopoNewEM_BCH_CORR_JET);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_BCH_CORR_JET_FORCELL", &jet_AntiKt4TopoNewEM_BCH_CORR_JET_FORCELL, &b_jet_AntiKt4TopoNewEM_BCH_CORR_JET_FORCELL);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_ENG_BAD_CELLS", &jet_AntiKt4TopoNewEM_ENG_BAD_CELLS, &b_jet_AntiKt4TopoNewEM_ENG_BAD_CELLS);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_N_BAD_CELLS", &jet_AntiKt4TopoNewEM_N_BAD_CELLS, &b_jet_AntiKt4TopoNewEM_N_BAD_CELLS);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_N_BAD_CELLS_CORR", &jet_AntiKt4TopoNewEM_N_BAD_CELLS_CORR, &b_jet_AntiKt4TopoNewEM_N_BAD_CELLS_CORR);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_BAD_CELLS_CORR_E", &jet_AntiKt4TopoNewEM_BAD_CELLS_CORR_E, &b_jet_AntiKt4TopoNewEM_BAD_CELLS_CORR_E);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_NumTowers", &jet_AntiKt4TopoNewEM_NumTowers, &b_jet_AntiKt4TopoNewEM_NumTowers);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_ootFracCells5", &jet_AntiKt4TopoNewEM_ootFracCells5, &b_jet_AntiKt4TopoNewEM_ootFracCells5);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_ootFracCells10", &jet_AntiKt4TopoNewEM_ootFracCells10, &b_jet_AntiKt4TopoNewEM_ootFracCells10);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_ootFracClusters5", &jet_AntiKt4TopoNewEM_ootFracClusters5, &b_jet_AntiKt4TopoNewEM_ootFracClusters5);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_ootFracClusters10", &jet_AntiKt4TopoNewEM_ootFracClusters10, &b_jet_AntiKt4TopoNewEM_ootFracClusters10);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_SamplingMax", &jet_AntiKt4TopoNewEM_SamplingMax, &b_jet_AntiKt4TopoNewEM_SamplingMax);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_fracSamplingMax", &jet_AntiKt4TopoNewEM_fracSamplingMax, &b_jet_AntiKt4TopoNewEM_fracSamplingMax);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_hecf", &jet_AntiKt4TopoNewEM_hecf, &b_jet_AntiKt4TopoNewEM_hecf);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_tgap3f", &jet_AntiKt4TopoNewEM_tgap3f, &b_jet_AntiKt4TopoNewEM_tgap3f);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_isUgly", &jet_AntiKt4TopoNewEM_isUgly, &b_jet_AntiKt4TopoNewEM_isUgly);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_isBadLooseMinus", &jet_AntiKt4TopoNewEM_isBadLooseMinus, &b_jet_AntiKt4TopoNewEM_isBadLooseMinus);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_isBadLoose", &jet_AntiKt4TopoNewEM_isBadLoose, &b_jet_AntiKt4TopoNewEM_isBadLoose);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_isBadMedium", &jet_AntiKt4TopoNewEM_isBadMedium, &b_jet_AntiKt4TopoNewEM_isBadMedium);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_isBadTight", &jet_AntiKt4TopoNewEM_isBadTight, &b_jet_AntiKt4TopoNewEM_isBadTight);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_emfrac", &jet_AntiKt4TopoNewEM_emfrac, &b_jet_AntiKt4TopoNewEM_emfrac);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_Offset", &jet_AntiKt4TopoNewEM_Offset, &b_jet_AntiKt4TopoNewEM_Offset);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_EMJES", &jet_AntiKt4TopoNewEM_EMJES, &b_jet_AntiKt4TopoNewEM_EMJES);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_EMJES_EtaCorr", &jet_AntiKt4TopoNewEM_EMJES_EtaCorr, &b_jet_AntiKt4TopoNewEM_EMJES_EtaCorr);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_EMJESnooffset", &jet_AntiKt4TopoNewEM_EMJESnooffset, &b_jet_AntiKt4TopoNewEM_EMJESnooffset);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_GCWJES", &jet_AntiKt4TopoNewEM_GCWJES, &b_jet_AntiKt4TopoNewEM_GCWJES);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_GCWJES_EtaCorr", &jet_AntiKt4TopoNewEM_GCWJES_EtaCorr, &b_jet_AntiKt4TopoNewEM_GCWJES_EtaCorr);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_CB", &jet_AntiKt4TopoNewEM_CB, &b_jet_AntiKt4TopoNewEM_CB);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_LCJES", &jet_AntiKt4TopoNewEM_LCJES, &b_jet_AntiKt4TopoNewEM_LCJES);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_emscale_E", &jet_AntiKt4TopoNewEM_emscale_E, &b_jet_AntiKt4TopoNewEM_emscale_E);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_emscale_pt", &jet_AntiKt4TopoNewEM_emscale_pt, &b_jet_AntiKt4TopoNewEM_emscale_pt);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_emscale_m", &jet_AntiKt4TopoNewEM_emscale_m, &b_jet_AntiKt4TopoNewEM_emscale_m);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_emscale_eta", &jet_AntiKt4TopoNewEM_emscale_eta, &b_jet_AntiKt4TopoNewEM_emscale_eta);
  fChain->SetBranchAddress("jet_AntiKt4TopoNewEM_emscale_phi", &jet_AntiKt4TopoNewEM_emscale_phi, &b_jet_AntiKt4TopoNewEM_emscale_phi);
  fChain->SetBranchAddress("mc_n", &mc_n, &b_mc_n);
  fChain->SetBranchAddress("mc_E", &mc_E, &b_mc_E);
  fChain->SetBranchAddress("mc_pt", &mc_pt, &b_mc_pt);
  fChain->SetBranchAddress("mc_m", &mc_m, &b_mc_m);
  fChain->SetBranchAddress("mc_eta", &mc_eta, &b_mc_eta);
  fChain->SetBranchAddress("mc_phi", &mc_phi, &b_mc_phi);
  fChain->SetBranchAddress("mc_px", &mc_px, &b_mc_px);
  fChain->SetBranchAddress("mc_py", &mc_py, &b_mc_py);
  fChain->SetBranchAddress("mc_pz", &mc_pz, &b_mc_pz);
  fChain->SetBranchAddress("mc_status", &mc_status, &b_mc_status);
  fChain->SetBranchAddress("mc_barcode", &mc_barcode, &b_mc_barcode);
  fChain->SetBranchAddress("mc_pdgId", &mc_pdgId, &b_mc_pdgId);
  fChain->SetBranchAddress("mc_charge", &mc_charge, &b_mc_charge);
  fChain->SetBranchAddress("mc_vx_x", &mc_vx_x, &b_mc_vx_x);
  fChain->SetBranchAddress("mc_vx_y", &mc_vx_y, &b_mc_vx_y);
  fChain->SetBranchAddress("mc_vx_z", &mc_vx_z, &b_mc_vx_z);
  fChain->SetBranchAddress("mc_vx_barcode", &mc_vx_barcode, &b_mc_vx_barcode);
  fChain->SetBranchAddress("mc_child_index", &mc_child_index, &b_mc_child_index);
  fChain->SetBranchAddress("mc_parent_index", &mc_parent_index, &b_mc_parent_index);
  fChain->SetBranchAddress("MET_Truth_NonInt_etx", &MET_Truth_NonInt_etx, &b_MET_Truth_NonInt_etx);
  fChain->SetBranchAddress("MET_Truth_NonInt_ety", &MET_Truth_NonInt_ety, &b_MET_Truth_NonInt_ety);
  fChain->SetBranchAddress("MET_Truth_Int_etx", &MET_Truth_Int_etx, &b_MET_Truth_Int_etx);
  fChain->SetBranchAddress("MET_Truth_Int_ety", &MET_Truth_Int_ety, &b_MET_Truth_Int_ety);
  fChain->SetBranchAddress("MET_Truth_IntCentral_etx", &MET_Truth_IntCentral_etx, &b_MET_Truth_IntCentral_etx);
  fChain->SetBranchAddress("MET_Truth_IntCentral_ety", &MET_Truth_IntCentral_ety, &b_MET_Truth_IntCentral_ety);
  fChain->SetBranchAddress("MET_Truth_IntFwd_etx", &MET_Truth_IntFwd_etx, &b_MET_Truth_IntFwd_etx);
  fChain->SetBranchAddress("MET_Truth_IntFwd_ety", &MET_Truth_IntFwd_ety, &b_MET_Truth_IntFwd_ety);
  fChain->SetBranchAddress("MET_Truth_IntOutCover_etx", &MET_Truth_IntOutCover_etx, &b_MET_Truth_IntOutCover_etx);
  fChain->SetBranchAddress("MET_Truth_IntOutCover_ety", &MET_Truth_IntOutCover_ety, &b_MET_Truth_IntOutCover_ety);
  fChain->SetBranchAddress("MET_Truth_IntMuons_etx", &MET_Truth_IntMuons_etx, &b_MET_Truth_IntMuons_etx);
  fChain->SetBranchAddress("MET_Truth_IntMuons_ety", &MET_Truth_IntMuons_ety, &b_MET_Truth_IntMuons_ety);
  fChain->SetBranchAddress("el_n", &el_n, &b_el_n);
  fChain->SetBranchAddress("el_E", &el_E, &b_el_E);
  fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
  fChain->SetBranchAddress("el_m", &el_m, &b_el_m);
  fChain->SetBranchAddress("el_px", &el_px, &b_el_px);
  fChain->SetBranchAddress("el_py", &el_py, &b_el_py);
  fChain->SetBranchAddress("el_pz", &el_pz, &b_el_pz);
  fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
  fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
  fChain->SetBranchAddress("el_status", &el_status, &b_el_status);
  fChain->SetBranchAddress("el_barcode", &el_barcode, &b_el_barcode);
  fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
  fChain->SetBranchAddress("mu_muid_n", &mu_muid_n, &b_mu_muid_n);
  fChain->SetBranchAddress("mu_muid_E", &mu_muid_E, &b_mu_muid_E);
  fChain->SetBranchAddress("mu_muid_pt", &mu_muid_pt, &b_mu_muid_pt);
  fChain->SetBranchAddress("mu_muid_m", &mu_muid_m, &b_mu_muid_m);
  fChain->SetBranchAddress("mu_muid_px", &mu_muid_px, &b_mu_muid_px);
  fChain->SetBranchAddress("mu_muid_py", &mu_muid_py, &b_mu_muid_py);
  fChain->SetBranchAddress("mu_muid_pz", &mu_muid_pz, &b_mu_muid_pz);
  fChain->SetBranchAddress("mu_muid_eta", &mu_muid_eta, &b_mu_muid_eta);
  fChain->SetBranchAddress("mu_muid_phi", &mu_muid_phi, &b_mu_muid_phi);
  fChain->SetBranchAddress("mu_muid_status", &mu_muid_status, &b_mu_muid_status);
  fChain->SetBranchAddress("mu_muid_barcode", &mu_muid_barcode, &b_mu_muid_barcode);
  fChain->SetBranchAddress("mu_muid_charge", &mu_muid_charge, &b_mu_muid_charge);
  fChain->SetBranchAddress("mu_staco_n", &mu_staco_n, &b_mu_staco_n);
  fChain->SetBranchAddress("mu_staco_E", &mu_staco_E, &b_mu_staco_E);
  fChain->SetBranchAddress("mu_staco_pt", &mu_staco_pt, &b_mu_staco_pt);
  fChain->SetBranchAddress("mu_staco_m", &mu_staco_m, &b_mu_staco_m);
  fChain->SetBranchAddress("mu_staco_px", &mu_staco_px, &b_mu_staco_px);
  fChain->SetBranchAddress("mu_staco_py", &mu_staco_py, &b_mu_staco_py);
  fChain->SetBranchAddress("mu_staco_pz", &mu_staco_pz, &b_mu_staco_pz);
  fChain->SetBranchAddress("mu_staco_eta", &mu_staco_eta, &b_mu_staco_eta);
  fChain->SetBranchAddress("mu_staco_phi", &mu_staco_phi, &b_mu_staco_phi);
  fChain->SetBranchAddress("mu_staco_status", &mu_staco_status, &b_mu_staco_status);
  fChain->SetBranchAddress("mu_staco_barcode", &mu_staco_barcode, &b_mu_staco_barcode);
  fChain->SetBranchAddress("mu_staco_charge", &mu_staco_charge, &b_mu_staco_charge);
  fChain->SetBranchAddress("tau_n", &tau_n, &b_tau_n);
  fChain->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
  fChain->SetBranchAddress("tau_m", &tau_m, &b_tau_m);
  fChain->SetBranchAddress("tau_eta", &tau_eta, &b_tau_eta);
  fChain->SetBranchAddress("tau_phi", &tau_phi, &b_tau_phi);
  fChain->SetBranchAddress("tau_status", &tau_status, &b_tau_status);
  fChain->SetBranchAddress("tau_barcode", &tau_barcode, &b_tau_barcode);
  fChain->SetBranchAddress("tau_charge", &tau_charge, &b_tau_charge);
  fChain->SetBranchAddress("ph_n", &ph_n, &b_ph_n);
  fChain->SetBranchAddress("ph_E", &ph_E, &b_ph_E);
  fChain->SetBranchAddress("ph_pt", &ph_pt, &b_ph_pt);
  fChain->SetBranchAddress("ph_m", &ph_m, &b_ph_m);
  fChain->SetBranchAddress("ph_px", &ph_px, &b_ph_px);
  fChain->SetBranchAddress("ph_py", &ph_py, &b_ph_py);
  fChain->SetBranchAddress("ph_pz", &ph_pz, &b_ph_pz);
  fChain->SetBranchAddress("ph_eta", &ph_eta, &b_ph_eta);
  fChain->SetBranchAddress("ph_phi", &ph_phi, &b_ph_phi);
  fChain->SetBranchAddress("ph_status", &ph_status, &b_ph_status);
  fChain->SetBranchAddress("ph_barcode", &ph_barcode, &b_ph_barcode);
  fChain->SetBranchAddress("el_truth_n", &el_truth_n, &b_el_truth_n);
  fChain->SetBranchAddress("el_truth_E", &el_truth_E, &b_el_truth_E);
  fChain->SetBranchAddress("el_truth_pt", &el_truth_pt, &b_el_truth_pt);
  fChain->SetBranchAddress("el_truth_m", &el_truth_m, &b_el_truth_m);
  fChain->SetBranchAddress("el_truth_px", &el_truth_px, &b_el_truth_px);
  fChain->SetBranchAddress("el_truth_py", &el_truth_py, &b_el_truth_py);
  fChain->SetBranchAddress("el_truth_pz", &el_truth_pz, &b_el_truth_pz);
  fChain->SetBranchAddress("el_truth_eta", &el_truth_eta, &b_el_truth_eta);
  fChain->SetBranchAddress("el_truth_phi", &el_truth_phi, &b_el_truth_phi);
  fChain->SetBranchAddress("el_truth_status", &el_truth_status, &b_el_truth_status);
  fChain->SetBranchAddress("el_truth_barcode", &el_truth_barcode, &b_el_truth_barcode);
  fChain->SetBranchAddress("el_truth_charge", &el_truth_charge, &b_el_truth_charge);
  fChain->SetBranchAddress("mu_muid_truth_n", &mu_muid_truth_n, &b_mu_muid_truth_n);
  fChain->SetBranchAddress("mu_muid_truth_E", &mu_muid_truth_E, &b_mu_muid_truth_E);
  fChain->SetBranchAddress("mu_muid_truth_pt", &mu_muid_truth_pt, &b_mu_muid_truth_pt);
  fChain->SetBranchAddress("mu_muid_truth_m", &mu_muid_truth_m, &b_mu_muid_truth_m);
  fChain->SetBranchAddress("mu_muid_truth_px", &mu_muid_truth_px, &b_mu_muid_truth_px);
  fChain->SetBranchAddress("mu_muid_truth_py", &mu_muid_truth_py, &b_mu_muid_truth_py);
  fChain->SetBranchAddress("mu_muid_truth_pz", &mu_muid_truth_pz, &b_mu_muid_truth_pz);
  fChain->SetBranchAddress("mu_muid_truth_eta", &mu_muid_truth_eta, &b_mu_muid_truth_eta);
  fChain->SetBranchAddress("mu_muid_truth_phi", &mu_muid_truth_phi, &b_mu_muid_truth_phi);
  fChain->SetBranchAddress("mu_muid_truth_status", &mu_muid_truth_status, &b_mu_muid_truth_status);
  fChain->SetBranchAddress("mu_muid_truth_barcode", &mu_muid_truth_barcode, &b_mu_muid_truth_barcode);
  fChain->SetBranchAddress("mu_muid_truth_charge", &mu_muid_truth_charge, &b_mu_muid_truth_charge);
  fChain->SetBranchAddress("mu_staco_truth_n", &mu_staco_truth_n, &b_mu_staco_truth_n);
  fChain->SetBranchAddress("mu_staco_truth_E", &mu_staco_truth_E, &b_mu_staco_truth_E);
  fChain->SetBranchAddress("mu_staco_truth_pt", &mu_staco_truth_pt, &b_mu_staco_truth_pt);
  fChain->SetBranchAddress("mu_staco_truth_m", &mu_staco_truth_m, &b_mu_staco_truth_m);
  fChain->SetBranchAddress("mu_staco_truth_px", &mu_staco_truth_px, &b_mu_staco_truth_px);
  fChain->SetBranchAddress("mu_staco_truth_py", &mu_staco_truth_py, &b_mu_staco_truth_py);
  fChain->SetBranchAddress("mu_staco_truth_pz", &mu_staco_truth_pz, &b_mu_staco_truth_pz);
  fChain->SetBranchAddress("mu_staco_truth_eta", &mu_staco_truth_eta, &b_mu_staco_truth_eta);
  fChain->SetBranchAddress("mu_staco_truth_phi", &mu_staco_truth_phi, &b_mu_staco_truth_phi);
  fChain->SetBranchAddress("mu_staco_truth_status", &mu_staco_truth_status, &b_mu_staco_truth_status);
  fChain->SetBranchAddress("mu_staco_truth_barcode", &mu_staco_truth_barcode, &b_mu_staco_truth_barcode);
  fChain->SetBranchAddress("mu_staco_truth_charge", &mu_staco_truth_charge, &b_mu_staco_truth_charge);
  fChain->SetBranchAddress("trueTau_n", &trueTau_n, &b_trueTau_n);
  fChain->SetBranchAddress("trueTau_pt", &trueTau_pt, &b_trueTau_pt);
  fChain->SetBranchAddress("trueTau_m", &trueTau_m, &b_trueTau_m);
  fChain->SetBranchAddress("trueTau_eta", &trueTau_eta, &b_trueTau_eta);
  fChain->SetBranchAddress("trueTau_phi", &trueTau_phi, &b_trueTau_phi);
  fChain->SetBranchAddress("trueTau_status", &trueTau_status, &b_trueTau_status);
  fChain->SetBranchAddress("trueTau_barcode", &trueTau_barcode, &b_trueTau_barcode);
  fChain->SetBranchAddress("trueTau_charge", &trueTau_charge, &b_trueTau_charge);
  fChain->SetBranchAddress("ph_truth_n", &ph_truth_n, &b_ph_truth_n);
  fChain->SetBranchAddress("ph_truth_E", &ph_truth_E, &b_ph_truth_E);
  fChain->SetBranchAddress("ph_truth_pt", &ph_truth_pt, &b_ph_truth_pt);
  fChain->SetBranchAddress("ph_truth_m", &ph_truth_m, &b_ph_truth_m);
  fChain->SetBranchAddress("ph_truth_px", &ph_truth_px, &b_ph_truth_px);
  fChain->SetBranchAddress("ph_truth_py", &ph_truth_py, &b_ph_truth_py);
  fChain->SetBranchAddress("ph_truth_pz", &ph_truth_pz, &b_ph_truth_pz);
  fChain->SetBranchAddress("ph_truth_eta", &ph_truth_eta, &b_ph_truth_eta);
  fChain->SetBranchAddress("ph_truth_phi", &ph_truth_phi, &b_ph_truth_phi);
  fChain->SetBranchAddress("ph_truth_status", &ph_truth_status, &b_ph_truth_status);
  fChain->SetBranchAddress("ph_truth_barcode", &ph_truth_barcode, &b_ph_truth_barcode);
  Notify();
}

// -----------------------------------------------------------------------------
Bool_t TruthNtuple::TruthNtupleLooper::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

// -----------------------------------------------------------------------------
void TruthNtuple::TruthNtupleLooper::Loop()
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  ProgressBar progress_bar(nentries, 100);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    progress_bar.checkProgress(jentry);

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    processEvent();
  }
}

// -----------------------------------------------------------------------------
void TruthNtuple::TruthNtupleLooper::processEvent()
{
  // do nothing
  for (int el_index = 0; el_index != el_n; ++el_index) {
    Electron test_el = Electron(this, el_index);
    std::cout << "electron parent: " << test_el.getParentPdgid() << "\n";
    if (test_el.getParentPdgid() == 0) std::cout << "\tel pt: " << test_el.getPt() << "\n";
  }
  for (int mu_index = 0; mu_index != mu_staco_n; ++mu_index) {
    Muon test_mu = Muon(this, mu_index);
    std::cout << "muon parent: " << test_mu.getParentPdgid() << "\n";
    if (test_mu.getParentPdgid() == 0) std::cout << "\tmu pt: " << test_mu.getPt() << "\n";
  }
  for (int jet_index = 0; jet_index != jet_AntiKt4TruthJets_n; ++jet_index) {
    Jet test_jet = Jet(this, jet_index);
    std::cout << "is b-jet: " << test_jet.getIsBJet() << "\n";
  }
}
