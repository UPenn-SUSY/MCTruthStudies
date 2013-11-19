#include "BMinusLCutflow/include/BMinusLCutflow.h"

#include <iostream>
#include <math.h>

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include "TruthNtupleLooper/include/ObjectDefs.h"

// -----------------------------------------------------------------------------
BMinusL::Cutflow::Cutflow(TTree* tree) : TruthNtuple::TruthNtupleLooper(tree)
{
  // TODO set up hist handle classes
  m_h_num_lep    = new TH1D("h_num_lep"   , "num leptons        ; Lepton multiplicity; Entries        ", 6, -0.5, 5.5);
  m_h_num_jet    = new TH1D("h_num_jet"   , "num b-jets         ; B-jet multiplicity;  Entries        ", 10, -0.5, 9.5);

  m_h_lep_pt_0   = new TH1D("h_lep_pt_0"  , "lep pt - leading   ; p_{T}^{0} [GeV];     Entries        ", 50, 0, 500);
  m_h_lep_pt_1   = new TH1D("h_lep_pt_1"  , "lep pt - subleading; p_{T}^{1} [GeV];     Entries        ", 50, 0, 500);
  m_h_lep_pt_2d  = new TH2D("h_lep_pt_2d" , "lep pt             ; p_{T}^{0} [GeV];     p_{T}^{1} [GeV]", 50, 0, 500, 50, 0, 500);

  m_h_lep_eta_0   = new TH1D("h_lep_eta_0"  , "lep eta - leading   ; #eta^{0};     Entries ", 50, -5, 5);
  m_h_lep_eta_1   = new TH1D("h_lep_eta_1"  , "lep eta - subleading; #eta^{1};     Entries ", 50, -5, 5);
  m_h_lep_eta_2d  = new TH2D("h_lep_eta_2d" , "lep eta             ; #eta^{0};     #eta^{1}", 50, -5, 5, 50, -5, 5);

  m_h_lep_phi_0   = new TH1D("h_lep_phi_0"  , "lep phi - leading   ; #phi^{0};     Entries ", 32, -3.2, 3.2);
  m_h_lep_phi_1   = new TH1D("h_lep_phi_1"  , "lep phi - subleading; #phi^{1};     Entries ", 32, -3.2, 3.2);
  m_h_lep_phi_2d  = new TH2D("h_lep_phi_2d" , "lep phi             ; #phi^{0};     #phi^{1}", 32, -3.2, 3.2, 32, -3.2, 3.2);

  m_h_jet_pt_0   = new TH1D("h_jet_pt_0"  , "jet pt - leading   ; p_{T}^{0} [GeV];     Entries", 50, 0, 500);
  m_h_jet_pt_1   = new TH1D("h_jet_pt_1"  , "jet pt - subleading; p_{T}^{1} [GeV];     Entries", 50, 0, 500);
  m_h_jet_pt_2d  = new TH2D("h_jet_pt_2d" , "jet pt             ; p_{T}^{0} [GeV];     p_{T}^{1} [GeV]", 50, 0, 500, 50, 0, 500);

  m_h_jet_eta_0   = new TH1D("h_jet_eta_0"  , "jet eta - leading   ; #eta^{0};     Entries ", 50, -5, 5);
  m_h_jet_eta_1   = new TH1D("h_jet_eta_1"  , "jet eta - subleading; #eta^{1};     Entries ", 50, -5, 5);
  m_h_jet_eta_2d  = new TH2D("h_jet_eta_2d" , "jet eta             ; #eta^{0};     #eta^{1}", 50, -5, 5, 50, -5, 5);

  m_h_jet_phi_0   = new TH1D("h_jet_phi_0"  , "jet phi - leading   ; #phi^{0};     Entries ", 32, -3.2, 3.2);
  m_h_jet_phi_1   = new TH1D("h_jet_phi_1"  , "jet phi - subleading; #phi^{1};     Entries ", 32, -3.2, 3.2);
  m_h_jet_phi_2d  = new TH2D("h_jet_phi_2d" , "jet phi             ; #phi^{0};     #phi^{1}", 32, -3.2, 3.2, 32, -3.2, 3.2);

  m_h_mbl_truth  = new TH1D("h_mbl_truth" , "mbl truth          ; m_{bl} [GeV];        Entries", 50, 0, 500);
  m_h_mbl_paired = new TH1D("h_mbl_paired", "mbl paired         ; m_{bl} [GeV];        Entries", 50, 0, 500);
}

// -----------------------------------------------------------------------------
BMinusL::Cutflow::~Cutflow()
{
  // do nothing
}

// -----------------------------------------------------------------------------
void BMinusL::Cutflow::clearObjects()
{
  TruthNtuple::TruthNtupleLooper::clearObjects();

  m_flavor_channel = TruthNtuple::FLAVOR_NONE;

  m_selected_el.clear();
  m_selected_mu.clear();
  m_selected_jet.clear();

  m_daughter_el.clear();
  m_daughter_mu.clear();
  m_daughter_jet.clear();
}

// -----------------------------------------------------------------------------
void BMinusL::Cutflow::processEvent()
{
  doObjectSelection();

  size_t num_el = m_daughter_el.size();
  size_t num_mu = m_daughter_mu.size();
  size_t num_jet = m_daughter_jet.size();

  if (num_el == 2 && num_mu == 0) {
    m_flavor_channel = TruthNtuple::FLAVOR_EE;
  }
  else if (num_el == 0 && num_mu == 2) {
    m_flavor_channel = TruthNtuple::FLAVOR_MM;
  }
  else if (num_el == 1 && num_mu == 1) {
    if (m_daughter_el.at(0)->getPt() >= m_daughter_mu.at(0)->getPt())
      m_flavor_channel = TruthNtuple::FLAVOR_EM;
    else
      m_flavor_channel = TruthNtuple::FLAVOR_ME;
  }
  else {
    m_flavor_channel = TruthNtuple::FLAVOR_NONE;
  }

  m_h_num_lep->Fill(num_el + num_mu);
  m_h_num_jet->Fill(num_jet);

  double lep_0_pt = 0;
  double lep_1_pt = 0;

  double lep_0_eta = 0;
  double lep_1_eta = 0;

  double lep_0_phi = 0;
  double lep_1_phi = 0;

  if (m_flavor_channel == TruthNtuple::FLAVOR_EE) {
    lep_0_pt = m_daughter_el.at(0)->getPt()/1.e3;
    lep_1_pt = m_daughter_el.at(1)->getPt()/1.e3;

    lep_0_eta = m_daughter_el.at(0)->getEta();
    lep_1_eta = m_daughter_el.at(1)->getEta();

    lep_0_phi = m_daughter_el.at(0)->getPhi();
    lep_1_phi = m_daughter_el.at(1)->getPhi();
  }
  else if (m_flavor_channel == TruthNtuple::FLAVOR_MM) {
    lep_0_pt = m_daughter_mu.at(0)->getPt()/1.e3;
    lep_1_pt = m_daughter_mu.at(1)->getPt()/1.e3;

    lep_0_eta = m_daughter_mu.at(0)->getEta();
    lep_1_eta = m_daughter_mu.at(1)->getEta();

    lep_0_phi = m_daughter_mu.at(0)->getPhi();
    lep_1_phi = m_daughter_mu.at(1)->getPhi();
  }
  else if (  m_flavor_channel == TruthNtuple::FLAVOR_EM
          || m_flavor_channel == TruthNtuple::FLAVOR_ME
          ) {
    lep_0_pt = m_daughter_el.at(0)->getPt()/1.e3;
    lep_1_pt = m_daughter_mu.at(0)->getPt()/1.e3;

    lep_0_eta = m_daughter_el.at(0)->getEta();
    lep_1_eta = m_daughter_mu.at(0)->getEta();

    lep_0_phi = m_daughter_el.at(0)->getPhi();
    lep_1_phi = m_daughter_mu.at(0)->getPhi();
  }

  if (lep_0_pt < lep_1_pt) {
    double tmp_pt = lep_0_pt;
    lep_0_pt = lep_1_pt;
    lep_1_pt = tmp_pt;

    double tmp_eta = lep_0_eta;
    lep_0_eta = lep_1_eta;
    lep_1_eta = tmp_eta;

    double tmp_phi = lep_0_phi;
    lep_0_phi = lep_1_phi;
    lep_1_phi = tmp_phi;
  }

  m_h_lep_pt_0->Fill(lep_0_pt);
  m_h_lep_pt_1->Fill(lep_1_pt);
  m_h_lep_pt_2d->Fill(lep_0_pt, lep_1_pt);

  m_h_lep_eta_0->Fill(lep_0_eta);
  m_h_lep_eta_1->Fill(lep_1_eta);
  m_h_lep_eta_2d->Fill(lep_0_eta, lep_1_eta);

  m_h_lep_phi_0->Fill(lep_0_phi);
  m_h_lep_phi_1->Fill(lep_1_phi);
  m_h_lep_phi_2d->Fill(lep_0_phi, lep_1_phi);

  double jet_0_pt = 0;
  double jet_1_pt = 0;

  double jet_0_eta = 0;
  double jet_1_eta = 0;

  double jet_0_phi = 0;
  double jet_1_phi = 0;

  if (num_jet > 0) {
    jet_0_pt  = m_daughter_jet.at(0)->getPt()/1.e3;
    jet_0_eta = m_daughter_jet.at(0)->getEta();
    jet_0_phi = m_daughter_jet.at(0)->getPhi();

    m_h_jet_pt_0->Fill(jet_0_pt);
    m_h_jet_eta_0->Fill(jet_0_eta);
    m_h_jet_phi_0->Fill(jet_0_phi);
  }
  if (num_jet > 1) {
    jet_1_pt = m_daughter_jet.at(1)->getPt()/1.e3;
    jet_1_eta = m_daughter_jet.at(1)->getEta();
    jet_1_phi = m_daughter_jet.at(1)->getPhi();

    m_h_jet_pt_1->Fill(jet_1_pt);
    m_h_jet_eta_1->Fill(jet_1_eta);
    m_h_jet_phi_1->Fill(jet_1_phi);

    m_h_jet_pt_2d->Fill( jet_0_pt , jet_1_pt );
    m_h_jet_eta_2d->Fill(jet_0_eta, jet_1_eta);
    m_h_jet_phi_2d->Fill(jet_0_phi, jet_1_phi);
  }
}

// -----------------------------------------------------------------------------
void BMinusL::Cutflow::writeToFile()
{
  TFile* f = new TFile("output_hists.root", "RECREATE");
  f->cd();

  m_h_num_lep->Write();
  m_h_num_jet->Write();

  m_h_lep_pt_0->Write();
  m_h_lep_pt_1->Write();
  m_h_lep_pt_2d->Write();

  m_h_lep_eta_0->Write();
  m_h_lep_eta_1->Write();
  m_h_lep_eta_2d->Write();

  m_h_lep_phi_0->Write();
  m_h_lep_phi_1->Write();
  m_h_lep_phi_2d->Write();

  m_h_jet_pt_0->Write();
  m_h_jet_pt_1->Write();
  m_h_jet_pt_2d->Write();

  m_h_jet_eta_0->Write();
  m_h_jet_eta_1->Write();
  m_h_jet_eta_2d->Write();

  m_h_jet_phi_0->Write();
  m_h_jet_phi_1->Write();
  m_h_jet_phi_2d->Write();

  m_h_mbl_truth->Write();
  m_h_mbl_paired->Write();

  f->Close();
}

// -----------------------------------------------------------------------------
void BMinusL::Cutflow::doObjectSelection()
{
  std::vector<TruthNtuple::Electron*> tmp_el;
  std::vector<TruthNtuple::Muon*>     tmp_mu;
  std::vector<TruthNtuple::Jet*>      tmp_jet;

  m_selected_el.resize(m_el_list.size());
  for (size_t el_it = 0; el_it != m_el_list.size(); ++el_it) {
    m_selected_el.push_back(&m_el_list.at(el_it));

    if (fabs(m_el_list.at(el_it).getParentPdgid()) > 1e6)
      m_daughter_el.push_back(&m_el_list.at(el_it));
  }

  m_selected_mu.resize(m_mu_list.size());
  for (size_t mu_it = 0; mu_it != m_mu_list.size(); ++mu_it) {
    m_selected_mu.push_back(&m_mu_list.at(mu_it));

    if (fabs(m_mu_list.at(mu_it).getParentPdgid()) > 1e6)
      m_daughter_mu.push_back(&m_mu_list.at(mu_it));
  }

  m_selected_jet.resize(m_jet_list.size());
  for (size_t jet_it = 0; jet_it != m_jet_list.size(); ++jet_it) {
    if (!m_jet_list.at(jet_it).getIsBJet())
      continue;

    m_selected_jet.push_back(&m_jet_list.at(jet_it));

    if (fabs(m_jet_list.at(jet_it).getParentPdgid()) > 1e6)
      m_daughter_jet.push_back(&m_jet_list.at(jet_it));
  }
}
