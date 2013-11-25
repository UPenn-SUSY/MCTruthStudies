#include "BMinusLCutflow/include/BMinusLCutflow.h"

#include <iostream>
#include <math.h>

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include "TruthNtupleLooper/include/ObjectDefs.h"

#include "HistogramHandlers/include/HistogramHandlers.h"

// -----------------------------------------------------------------------------
BMinusL::Cutflow::Cutflow(TTree* tree) : TruthNtuple::TruthNtupleLooper(tree)
{
  // construct histogram list
  m_histograms.push_back(new HistogramHandlers::FlavorChannel());
  m_histograms.push_back(new HistogramHandlers::ObjectMultiplicity());
  m_histograms.push_back(new HistogramHandlers::LeptonPt());
  m_histograms.push_back(new HistogramHandlers::LeptonEta());
  m_histograms.push_back(new HistogramHandlers::LeptonPhi());
  m_histograms.push_back(new HistogramHandlers::JetPt());
  m_histograms.push_back(new HistogramHandlers::JetEta());
  m_histograms.push_back(new HistogramHandlers::JetPhi());
  m_histograms.push_back(new HistogramHandlers::Met());
  m_histograms.push_back(new HistogramHandlers::Mll());
  m_histograms.push_back(new HistogramHandlers::Mbl());
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

  m_met.clear();
}

// -----------------------------------------------------------------------------
void BMinusL::Cutflow::processEvent()
{
  doObjectSelection();

  size_t num_el  = m_daughter_el.size();
  size_t num_mu  = m_daughter_mu.size();
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

  std::cout << "========================================"
            << "\nevent number: " << EventNumber
            << "\n\tnum el: " << num_el
            << "\n\tnum mu: " << num_mu
            << "\n\tnum jet: " << num_jet
            << "\n----------------------------------------"
            << "\n";

  for (size_t el_it = 0; el_it != num_el; ++el_it) {
    m_daughter_el.at(el_it)->print(this);
  }
  for (size_t mu_it = 0; mu_it != num_mu; ++mu_it) {
    m_daughter_mu.at(mu_it)->print(this);
  }
  for (size_t jet_it = 0; jet_it != num_jet; ++jet_it) {
    m_daughter_jet.at(jet_it)->print(this);
  }

  size_t num_hists = m_histograms.size();
  for (size_t hist_it = 0; hist_it != num_hists; ++hist_it) {
    m_histograms.at(hist_it);
    m_histograms.at(hist_it)->Fill( m_flavor_channel
                                  , m_daughter_el
                                  , m_daughter_mu
                                  , m_daughter_jet
                                  , m_met
                                  );
  }
}

// -----------------------------------------------------------------------------
void BMinusL::Cutflow::writeToFile()
{
  TFile* f = new TFile("output_hists.root", "RECREATE");
  f->cd();

  size_t num_hists = m_histograms.size();
  for (size_t hist_it = 0; hist_it != num_hists; ++hist_it) {
    m_histograms.at(hist_it)->write(f);
  }
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

  m_met.setMetNoint(MET_Truth_NonInt_etx, MET_Truth_NonInt_ety);
  m_met.calculateMetRelNoint( m_daughter_el, m_daughter_mu, m_daughter_jet);
}
