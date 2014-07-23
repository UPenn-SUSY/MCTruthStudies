#include "TruthNtupleLooper/include/TruthNtupleLooper.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <dirent.h>
#include <math.h>
#include <algorithm>

// #include "TH2.h"

#include "TruthNtupleLooper/include/TruthNtupleEnums.h"
#include "TruthNtupleLooper/include/ObjectDefs.h"
#include "TruthNtupleLooper/include/OverlapRemoval.h"
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
  jet_AntiKt4TruthJets_E = 0;
  jet_AntiKt4TruthJets_pt = 0;
  jet_AntiKt4TruthJets_m = 0;
  jet_AntiKt4TruthJets_eta = 0;
  jet_AntiKt4TruthJets_phi = 0;
  jet_AntiKt4TruthJets_flavor_partonDR = 0;
  jet_AntiKt4TruthJets_flavor_partonFlavor = 0;
  jet_AntiKt4TruthJets_flavor_hadronFlavor = 0;
  jet_AntiKt4TruthJets_flavor_hadronPDGID = 0;
  jet_AntiKt4TruthJets_WZ_E = 0;
  jet_AntiKt4TruthJets_WZ_pt = 0;
  jet_AntiKt4TruthJets_WZ_m = 0;
  jet_AntiKt4TruthJets_WZ_eta = 0;
  jet_AntiKt4TruthJets_WZ_phi = 0;
  jet_AntiKt4TruthJets_WZ_flavor_partonDR = 0;
  jet_AntiKt4TruthJets_WZ_flavor_partonFlavor = 0;
  jet_AntiKt4TruthJets_WZ_flavor_hadronFlavor = 0;
  jet_AntiKt4TruthJets_WZ_flavor_hadronPDGID = 0;
  jet_AntiKt6TruthJets_WZ_E = 0;
  jet_AntiKt6TruthJets_WZ_pt = 0;
  jet_AntiKt6TruthJets_WZ_m = 0;
  jet_AntiKt6TruthJets_WZ_eta = 0;
  jet_AntiKt6TruthJets_WZ_phi = 0;
  jet_AntiKt6TruthJets_WZ_flavor_partonDR = 0;
  jet_AntiKt6TruthJets_WZ_flavor_partonFlavor = 0;
  jet_AntiKt6TruthJets_WZ_flavor_hadronFlavor = 0;
  jet_AntiKt6TruthJets_WZ_flavor_hadronPDGID = 0;
  mc_pt = 0;
  mc_m = 0;
  mc_eta = 0;
  mc_phi = 0;
  mc_status = 0;
  mc_barcode = 0;
  mc_pdgId = 0;
  mc_charge = 0;
  mc_parents = 0;
  mc_children = 0;
  mc_vx_x = 0;
  mc_vx_y = 0;
  mc_vx_z = 0;
  mc_vx_barcode = 0;
  mc_child_index = 0;
  mc_parent_index = 0;
  el_pt = 0;
  el_m = 0;
  el_eta = 0;
  el_phi = 0;
  el_status = 0;
  el_barcode = 0;
  el_charge = 0;
  el_parent_index = 0;
  mu_pt = 0;
  mu_m = 0;
  mu_eta = 0;
  mu_phi = 0;
  mu_status = 0;
  mu_barcode = 0;
  mu_charge = 0;
  mu_parent_index = 0;
  tau_pt = 0;
  tau_m = 0;
  tau_eta = 0;
  tau_phi = 0;
  tau_status = 0;
  tau_barcode = 0;
  tau_charge = 0;
  tau_parent_index = 0;
  tau_decay_index = 0;
  ph_pt = 0;
  ph_m = 0;
  ph_eta = 0;
  ph_phi = 0;
  ph_status = 0;
  ph_barcode = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
  fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
  fChain->SetBranchAddress("mc_channel_number", &mc_channel_number, &b_mc_channel_number);
  fChain->SetBranchAddress("mc_event_number", &mc_event_number, &b_mc_event_number);
  fChain->SetBranchAddress("mc_event_weight", &mc_event_weight, &b_mc_event_weight);
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
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_n", &jet_AntiKt4TruthJets_n, &b_jet_AntiKt4TruthJets_n);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_E", &jet_AntiKt4TruthJets_E, &b_jet_AntiKt4TruthJets_E);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_pt", &jet_AntiKt4TruthJets_pt, &b_jet_AntiKt4TruthJets_pt);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_m", &jet_AntiKt4TruthJets_m, &b_jet_AntiKt4TruthJets_m);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_eta", &jet_AntiKt4TruthJets_eta, &b_jet_AntiKt4TruthJets_eta);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_phi", &jet_AntiKt4TruthJets_phi, &b_jet_AntiKt4TruthJets_phi);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_flavor_partonDR", &jet_AntiKt4TruthJets_flavor_partonDR, &b_jet_AntiKt4TruthJets_flavor_partonDR);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_flavor_partonFlavor", &jet_AntiKt4TruthJets_flavor_partonFlavor, &b_jet_AntiKt4TruthJets_flavor_partonFlavor);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_flavor_hadronFlavor", &jet_AntiKt4TruthJets_flavor_hadronFlavor, &b_jet_AntiKt4TruthJets_flavor_hadronFlavor);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_flavor_hadronPDGID", &jet_AntiKt4TruthJets_flavor_hadronPDGID, &b_jet_AntiKt4TruthJets_flavor_hadronPDGID);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_WZ_n", &jet_AntiKt4TruthJets_WZ_n, &b_jet_AntiKt4TruthJets_WZ_n);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_WZ_E", &jet_AntiKt4TruthJets_WZ_E, &b_jet_AntiKt4TruthJets_WZ_E);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_WZ_pt", &jet_AntiKt4TruthJets_WZ_pt, &b_jet_AntiKt4TruthJets_WZ_pt);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_WZ_m", &jet_AntiKt4TruthJets_WZ_m, &b_jet_AntiKt4TruthJets_WZ_m);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_WZ_eta", &jet_AntiKt4TruthJets_WZ_eta, &b_jet_AntiKt4TruthJets_WZ_eta);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_WZ_phi", &jet_AntiKt4TruthJets_WZ_phi, &b_jet_AntiKt4TruthJets_WZ_phi);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_WZ_flavor_partonDR", &jet_AntiKt4TruthJets_WZ_flavor_partonDR, &b_jet_AntiKt4TruthJets_WZ_flavor_partonDR);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_WZ_flavor_partonFlavor", &jet_AntiKt4TruthJets_WZ_flavor_partonFlavor, &b_jet_AntiKt4TruthJets_WZ_flavor_partonFlavor);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_WZ_flavor_hadronFlavor", &jet_AntiKt4TruthJets_WZ_flavor_hadronFlavor, &b_jet_AntiKt4TruthJets_WZ_flavor_hadronFlavor);
  fChain->SetBranchAddress("jet_AntiKt4TruthJets_WZ_flavor_hadronPDGID", &jet_AntiKt4TruthJets_WZ_flavor_hadronPDGID, &b_jet_AntiKt4TruthJets_WZ_flavor_hadronPDGID);
  fChain->SetBranchAddress("jet_AntiKt6TruthJets_WZ_n", &jet_AntiKt6TruthJets_WZ_n, &b_jet_AntiKt6TruthJets_WZ_n);
  fChain->SetBranchAddress("jet_AntiKt6TruthJets_WZ_E", &jet_AntiKt6TruthJets_WZ_E, &b_jet_AntiKt6TruthJets_WZ_E);
  fChain->SetBranchAddress("jet_AntiKt6TruthJets_WZ_pt", &jet_AntiKt6TruthJets_WZ_pt, &b_jet_AntiKt6TruthJets_WZ_pt);
  fChain->SetBranchAddress("jet_AntiKt6TruthJets_WZ_m", &jet_AntiKt6TruthJets_WZ_m, &b_jet_AntiKt6TruthJets_WZ_m);
  fChain->SetBranchAddress("jet_AntiKt6TruthJets_WZ_eta", &jet_AntiKt6TruthJets_WZ_eta, &b_jet_AntiKt6TruthJets_WZ_eta);
  fChain->SetBranchAddress("jet_AntiKt6TruthJets_WZ_phi", &jet_AntiKt6TruthJets_WZ_phi, &b_jet_AntiKt6TruthJets_WZ_phi);
  fChain->SetBranchAddress("jet_AntiKt6TruthJets_WZ_flavor_partonDR", &jet_AntiKt6TruthJets_WZ_flavor_partonDR, &b_jet_AntiKt6TruthJets_WZ_flavor_partonDR);
  fChain->SetBranchAddress("jet_AntiKt6TruthJets_WZ_flavor_partonFlavor", &jet_AntiKt6TruthJets_WZ_flavor_partonFlavor, &b_jet_AntiKt6TruthJets_WZ_flavor_partonFlavor);
  fChain->SetBranchAddress("jet_AntiKt6TruthJets_WZ_flavor_hadronFlavor", &jet_AntiKt6TruthJets_WZ_flavor_hadronFlavor, &b_jet_AntiKt6TruthJets_WZ_flavor_hadronFlavor);
  fChain->SetBranchAddress("jet_AntiKt6TruthJets_WZ_flavor_hadronPDGID", &jet_AntiKt6TruthJets_WZ_flavor_hadronPDGID, &b_jet_AntiKt6TruthJets_WZ_flavor_hadronPDGID);
  fChain->SetBranchAddress("mc_n", &mc_n, &b_mc_n);
  fChain->SetBranchAddress("mc_pt", &mc_pt, &b_mc_pt);
  fChain->SetBranchAddress("mc_m", &mc_m, &b_mc_m);
  fChain->SetBranchAddress("mc_eta", &mc_eta, &b_mc_eta);
  fChain->SetBranchAddress("mc_phi", &mc_phi, &b_mc_phi);
  fChain->SetBranchAddress("mc_status", &mc_status, &b_mc_status);
  fChain->SetBranchAddress("mc_barcode", &mc_barcode, &b_mc_barcode);
  fChain->SetBranchAddress("mc_pdgId", &mc_pdgId, &b_mc_pdgId);
  fChain->SetBranchAddress("mc_charge", &mc_charge, &b_mc_charge);
  fChain->SetBranchAddress("mc_parents", &mc_parents, &b_mc_parents);
  fChain->SetBranchAddress("mc_children", &mc_children, &b_mc_children);
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
  fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
  fChain->SetBranchAddress("el_m", &el_m, &b_el_m);
  fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
  fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
  fChain->SetBranchAddress("el_status", &el_status, &b_el_status);
  fChain->SetBranchAddress("el_barcode", &el_barcode, &b_el_barcode);
  fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
  fChain->SetBranchAddress("el_parent_index", &el_parent_index, &b_el_parent_index);
  fChain->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
  fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
  fChain->SetBranchAddress("mu_m", &mu_m, &b_mu_m);
  fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
  fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
  fChain->SetBranchAddress("mu_status", &mu_status, &b_mu_status);
  fChain->SetBranchAddress("mu_barcode", &mu_barcode, &b_mu_barcode);
  fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
  fChain->SetBranchAddress("mu_parent_index", &mu_parent_index, &b_mu_parent_index);
  fChain->SetBranchAddress("tau_n", &tau_n, &b_tau_n);
  fChain->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
  fChain->SetBranchAddress("tau_m", &tau_m, &b_tau_m);
  fChain->SetBranchAddress("tau_eta", &tau_eta, &b_tau_eta);
  fChain->SetBranchAddress("tau_phi", &tau_phi, &b_tau_phi);
  fChain->SetBranchAddress("tau_status", &tau_status, &b_tau_status);
  fChain->SetBranchAddress("tau_barcode", &tau_barcode, &b_tau_barcode);
  fChain->SetBranchAddress("tau_charge", &tau_charge, &b_tau_charge);
  fChain->SetBranchAddress("tau_parent_index", &tau_parent_index, &b_tau_parent_index);
  fChain->SetBranchAddress("tau_decay_index", &tau_decay_index, &b_tau_decay_index);
  fChain->SetBranchAddress("ph_n", &ph_n, &b_ph_n);
  fChain->SetBranchAddress("ph_pt", &ph_pt, &b_ph_pt);
  fChain->SetBranchAddress("ph_m", &ph_m, &b_ph_m);
  fChain->SetBranchAddress("ph_eta", &ph_eta, &b_ph_eta);
  fChain->SetBranchAddress("ph_phi", &ph_phi, &b_ph_phi);
  fChain->SetBranchAddress("ph_status", &ph_status, &b_ph_status);
  fChain->SetBranchAddress("ph_barcode", &ph_barcode, &b_ph_barcode);
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

  nentries = fChain->GetEntries();

  ProgressBar progress_bar(nentries, 100);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    progress_bar.checkProgress(jentry);
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    clearObjects();
    constructObjects();
    processEvent();
  }
}

// -----------------------------------------------------------------------------
void TruthNtuple::TruthNtupleLooper::clearObjects()
{
  m_particle_list.clear();
  m_el_list.clear();
  m_mu_list.clear();
  m_jet_list.clear();
  m_met.clear();
}

// -----------------------------------------------------------------------------
void TruthNtuple::TruthNtupleLooper::constructObjects()
{
  for (int mc_index = 0; mc_index != mc_n; ++mc_index) {
    Particle this_particle = Particle(this, mc_index);
    m_particle_list.push_back(this_particle);
  }
  for (int el_index = 0; el_index != el_n; ++el_index) {
    Electron this_el = Electron(this, el_index, true);
    m_el_list.push_back(this_el);
  }
  for (int mu_index = 0; mu_index != mu_n; ++mu_index) {
    Muon this_mu = Muon(this, mu_index, true);
    m_mu_list.push_back(this_mu);
  }
  for (int jet_index = 0; jet_index != jet_AntiKt4TruthJets_WZ_n; ++jet_index) {
    Jet this_jet = Jet(this, jet_index);
    m_jet_list.push_back(this_jet);
  }
  m_met.setMetNoint(MET_Truth_NonInt_etx, MET_Truth_NonInt_ety);
}

// -----------------------------------------------------------------------------
void TruthNtuple::TruthNtupleLooper::processEvent()
{
  std::vector<TruthNtuple::Electron*> el_baseline_list;
  std::vector<TruthNtuple::Electron*> el_signal_list;
  for (size_t el_it = 0; el_it != m_el_list.size(); ++el_it)
    el_baseline_list.push_back(&m_el_list.at(el_it));

  std::vector<TruthNtuple::Muon*> mu_baseline_list;
  std::vector<TruthNtuple::Muon*> mu_signal_list;
  for (size_t mu_it = 0; mu_it != m_mu_list.size(); ++mu_it)
    mu_baseline_list.push_back(&m_mu_list.at(mu_it));

  std::vector<TruthNtuple::Jet*> jet_baseline_list;
  std::vector<TruthNtuple::Jet*> jet_signal_list;
  for (size_t jet_it = 0; jet_it != m_jet_list.size(); ++jet_it)
    jet_baseline_list.push_back(&m_jet_list.at(jet_it));

  TruthNtuple::OverlapRemoval overlap_removal;
  overlap_removal.setDrEE(0.05);
  overlap_removal.setDrEJ(0.20);
  overlap_removal.setDrJE(0.40);
  overlap_removal.setDrJM(0.40);
  overlap_removal.setDrEM(0.01);
  overlap_removal.setDrMM(0.05);
  overlap_removal.doOverlapRemoval( el_baseline_list, mu_baseline_list, jet_baseline_list
                                  , el_signal_list  , mu_signal_list  , jet_signal_list
                                  );

  std::cout << "num el  -- (nominal): " << m_el_list.size() 
            << " -- (before): " << el_baseline_list.size()
            << " -- (after): "  << el_signal_list.size()
            << "\n";
  std::cout << "num mu  -- (nominal): " << m_mu_list.size()
            << " -- (before): " << mu_baseline_list.size()
            << " -- (after): "  << mu_signal_list.size()
            << "\n";
  std::cout << "num jet -- (nominal): " << m_jet_list.size()
            << " -- (before): " << jet_baseline_list.size()
            << " -- (after): " << jet_signal_list.size()
            << "\n";
}

// ------------------------------------------------------------------------------
void TruthNtuple::TruthNtupleLooper::cleanParticleList(
    std::vector<TruthNtuple::Particle*>& part_list, bool get_first_in_list)
{
  // total nubmer of particles before cleaning
  size_t num_part = part_list.size();

  // vector to tell us if this particle is to be removed
  std::vector<bool> to_remove(num_part, false);

  // loop over all particles to look for a particle who is its own parent
  for (size_t part_it = 0; part_it != num_part; ++part_it) {
    if (part_list.at(part_it)->getM() < 0) {
      to_remove.at(part_it) = true;
      continue;
    }

    // store the pdgid and parent index of this particle
    int this_pdgid = part_list.at(part_it)->getPdgid();
    int this_barcode = part_list.at(part_it)->getBarcode();
    int this_parent_index = part_list.at(part_it)->getImmediateParentMCIndex();
    int this_parent_barcode = part_list.at(part_it)->getImmediateParentBarcode();

    // loop over particles again to look for parent. if this particle and its
    // parent have the same pdgid, we will say they are the same and flag this
    // particle for removal
    for (size_t inner_part_it = 0; inner_part_it != num_part; ++inner_part_it) {
      if (part_list.at(inner_part_it)->getM() < 0 ) {
        continue;
      }

      if ( part_list.at(inner_part_it)->getMCIndex() == this_parent_index ) {
        // we found the parent, do the pdgid's match?
        int parent_pdgid = part_list.at(inner_part_it)->getPdgid();
        if ( parent_pdgid == this_pdgid ) {
          if (get_first_in_list) to_remove.at(part_it)       = true;
          else                   to_remove.at(inner_part_it) = true;
        }
        if ( parent_pdgid == -this_pdgid ) {
          to_remove.at(inner_part_it) = true;
        }
        // we already found the parent, we don't need to continue the inner loop
        break;
      }
    }
  }

  // loop backward through particles. if a particle is flagged for removal, remove it
  for (size_t part_it = 0; part_it != num_part; ++part_it) {
    size_t this_index = num_part - part_it - 1;
    if (to_remove.at(this_index)) {
      part_list.erase(part_list.begin() + this_index);
    }
  }
}
