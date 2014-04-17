#include "BMinusLCutflow/include/BMinusLCutflow.h"

#include <iostream>
#include <math.h>

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include "TruthNtupleLooper/include/ObjectDefs.h"

#include "HistogramHandlers/include/HistogramHandlers.h"
#include "BMinusLCutflow/include/BMinusLHistogramHandlers.h"

// -----------------------------------------------------------------------------
BMinusL::Cutflow::Cutflow(TTree* tree) : TruthNtuple::TruthNtupleLooper(tree)
{
  // construct histogram list
  m_histograms.push_back(new HistogramHandlers::FlavorChannel());
  m_histograms.push_back(new HistogramHandlers::ObjectMultiplicity());
  m_histograms.push_back(new HistogramHandlers::LeptonKinematics());
  m_histograms.push_back(new HistogramHandlers::JetKinematics());
  m_histograms.push_back(new HistogramHandlers::Met());
  m_histograms.push_back(new HistogramHandlers::Mll());
  // m_histograms.push_back(new HistogramHandlers::Mjl());

  m_h_mbl                = new HistogramHandlers::Mbl();
  m_h_bl_pair_kinematics = new HistogramHandlers::BLPairKinematics();
  m_h_quark_kinematics   = new HistogramHandlers::QuarkKinematics();
  m_h_stop_kinematics    = new HistogramHandlers::StopKinematics();

  m_sa_hists = new BMinusL::BMinusLStandAloneHistograms();
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

  m_truth_stops.clear();
  m_truth_electrons.clear();
  m_truth_muons.clear();
  m_truth_taus.clear();
  m_truth_b_quarks.clear();
  m_b_jets.clear();

  m_daughter_el.clear();
  m_daughter_mu.clear();
  m_daughter_tau.clear();
  m_daughter_b_quarks.clear();
  // m_daughter_jet.clear();

  m_leading_b_jets.clear();

  m_met.clear();
}

// -----------------------------------------------------------------------------
void BMinusL::Cutflow::processEvent()
{
  doObjectSelection();

  size_t num_el  = m_daughter_el.size();
  size_t num_mu  = m_daughter_mu.size();
  // size_t num_tau = m_daughter_tau.size();
  // size_t num_jet = m_daughter_jet.size();
  // size_t num_truth_b_quarks = m_daughter_b_quarks.size();
  // size_t num_b_jet = m_b_jets.size();

  // count the number of light leptons which come from tau decays
  size_t num_tau = 0;
  for (size_t el_it = 0 ; el_it != num_el; ++el_it) {
    if (fabs(m_daughter_el.at(el_it)->getParentPdgid()) == 15) {
      ++num_tau;
    }
  }
  for (size_t mu_it = 0 ; mu_it != num_mu; ++mu_it) {
    if (fabs(m_daughter_mu.at(mu_it)->getParentPdgid()) == 15) {
      ++num_tau;
    }
  }

  // if (num_tau == 0) return;

  // define the flavor channel based on the number of each lepton flavor
  // std::cout << "\n--------------------------------------------------------------------------------"
  //           << "\nEvent number: " << EventNumber
  //           << "\nchecking the flavor channel -- num_el: " << num_el << " -- num_mu: " << num_mu
  //           << "\n";
  if (num_el == 2 && num_mu == 0) {
    // std::cout << "\tthis event is EE\n";
    m_flavor_channel = TruthNtuple::FLAVOR_EE;
  }
  else if (num_el == 0 && num_mu == 2) {
    // std::cout << "\tthis event is MM\n";
    m_flavor_channel = TruthNtuple::FLAVOR_MM;
  }
  else if (num_el == 1 && num_mu == 1) {
    // std::cout << "\tthis event is EM or ME\n";
    if (m_daughter_el.at(0)->getPt() >= m_daughter_mu.at(0)->getPt()) {
      // std::cout << "\t\tthis event is really EM\n";
      m_flavor_channel = TruthNtuple::FLAVOR_EM;
    }
    else {
      // std::cout << "\t\tthis event is really ME\n";
      m_flavor_channel = TruthNtuple::FLAVOR_ME;
    }
  }
  else {
    // std::cout << "\tthis event is NONE\n";
    m_flavor_channel = TruthNtuple::FLAVOR_NONE;
  }

  // // TODO decide if we want this
  // if (  m_flavor_channel == TruthNtuple::FLAVOR_NONE
  //    || num_truth_b_quarks != 2
  //    )
  // {
  //   std::cout << "\nskipping event -- flavor: " << m_flavor_channel
  //             << " - num b quarks: " << num_truth_b_quarks
  //             << "\n";
  //   return;
  // }

  // print();

  // Fill "normal histograms"
  size_t num_hists = m_histograms.size();
  for (size_t hist_it = 0; hist_it != num_hists; ++hist_it) {
    m_histograms.at(hist_it)->Fill( m_flavor_channel
                                  , m_daughter_el
                                  , m_daughter_mu
                                  , m_leading_b_jets
                                  , m_met
                                  );
  }

  // fill special histograms
  m_h_mbl->FillSpecial( m_flavor_channel
                      , m_daughter_el
                      , m_daughter_mu
                      , m_leading_b_jets
                      );
  m_h_bl_pair_kinematics->FillSpecial( m_flavor_channel
                                     , m_daughter_el
                                     , m_daughter_mu
                                     , m_daughter_b_quarks
                                     );
  m_h_quark_kinematics->FillSpecial( m_flavor_channel
                                   , m_daughter_b_quarks
                                   );
  m_h_stop_kinematics->FillSpecial( m_flavor_channel
                                  , m_truth_stops
                                  );

  m_sa_hists->Fill(this);
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

  m_h_mbl->write(f);
  m_h_bl_pair_kinematics->write(f);
  m_h_quark_kinematics->write(f);
  m_h_stop_kinematics->write(f);
  m_sa_hists->write(f);
}

// -----------------------------------------------------------------------------
void BMinusL::Cutflow::doObjectSelection()
{
  // look for stops, leptons, and b quarks in full truth record
  for ( size_t particle_it = 0
      ; particle_it != m_particle_list.size()
      ; ++particle_it
      ) {
    if (fabs(m_particle_list.at(particle_it).getPdgid()) == 1e6+6) {
      m_truth_stops.push_back(&m_particle_list.at(particle_it));
    }
    if (  fabs(m_particle_list.at(particle_it).getPdgid()) == 11
       && (  m_particle_list.at(particle_it).getStatus() == 3
          || m_particle_list.at(particle_it).getStatus() == 1
          || m_particle_list.at(particle_it).getStatus() == 11 // herwig++
          )
       ) {
      m_truth_electrons.push_back(&m_particle_list.at(particle_it));
    }
    if (  fabs(m_particle_list.at(particle_it).getPdgid()) == 13
       && (  m_particle_list.at(particle_it).getStatus() == 3
          || m_particle_list.at(particle_it).getStatus() == 1
          || m_particle_list.at(particle_it).getStatus() == 11 // herwig++
          )
       ) {
      m_truth_muons.push_back(&m_particle_list.at(particle_it));
    }
    if (  fabs(m_particle_list.at(particle_it).getPdgid()) == 5
       && (  m_particle_list.at(particle_it).getStatus() == 3
          || m_particle_list.at(particle_it).getStatus() == 11 // herwig++
          )
       ) {
      m_truth_b_quarks.push_back(&m_particle_list.at(particle_it));
    }
  }

  // pick electrons coming from a stop
  m_daughter_el.reserve(m_truth_electrons.size());
  for (size_t el_it = 0; el_it != m_truth_electrons.size(); ++el_it) {
    if (  (  (  m_truth_electrons.at(el_it)->getStatus() == 3
             || m_truth_electrons.at(el_it)->getStatus() == 11 // herwig++
             )
          && fabs(m_truth_electrons.at(el_it)->getParentPdgid()) == 1e6+6
          )
        // TODO reimpliment the stop decays thgrouh a tau
       // || (  (  m_truth_electrons.at(el_it)->getStatus() == 1
       //       || m_truth_electrons.at(el_it)->getStatus() == 11 // herwig++
       //       )
       //    && fabs(m_truth_electrons.at(el_it)->getParentPdgid()) == 15
       //    )
       ) {
      m_daughter_el.push_back(m_truth_electrons.at(el_it));
    }
  }

  // pick muons coming from a stop
  m_daughter_mu.reserve(m_truth_muons.size());
  for (size_t mu_it = 0; mu_it != m_truth_muons.size(); ++mu_it) {
    if (  (  (  m_truth_muons.at(mu_it)->getStatus() == 3
             || m_truth_muons.at(mu_it)->getStatus() == 11 // herwig++
             )
          && fabs(m_truth_muons.at(mu_it)->getParentPdgid()) == 1e6+6
          )
        // TODO reimpliment the stop decays thgrouh a tau
       // || (  (  m_truth_muons.at(mu_it)->getStatus() == 1
       //       || m_truth_muons.at(mu_it)->getStatus() == 11 // herwig++
       //       )
       //    && fabs(m_truth_muons.at(mu_it)->getParentPdgid()) == 15
       //    )
       ) {
      m_daughter_mu.push_back(m_truth_muons.at(mu_it));
    }
  }

  // pick b quarks coming from a stop
  m_daughter_b_quarks.reserve(m_truth_b_quarks.size());
  for (size_t b_quarks_it = 0; b_quarks_it != m_truth_b_quarks.size(); ++b_quarks_it) {
    if (fabs(m_truth_b_quarks.at(b_quarks_it)->getParentPdgid()) == 1e6+6)
      m_daughter_b_quarks.push_back(m_truth_b_quarks.at(b_quarks_it));
  }

  size_t num_stops_before_cleaning = m_truth_stops.size();
  cleanParticleList(m_truth_stops);
  size_t num_stops_after_cleaning = m_truth_stops.size();
  // std::cout << "num stops: " << m_truth_stops.size()
  // std::cout << "num stops before cleaning: " << num_stops_before_cleaning
  //           << " -- num_stops after cleaning: " << num_stops_after_cleaning
  //           << " -- num daughter e: " << m_daughter_el.size()
  //           << " -- num daughter mu: " << m_daughter_mu.size()
  //           << " -- num_daughter b: " << m_daughter_b_quarks.size()
  //           << "\n";

  // pick b jets
  m_b_jets.reserve(m_jet_list.size());
  for (size_t jet_it = 0; jet_it != m_jet_list.size(); ++jet_it) {
    if (  jet_it > 0
       && m_jet_list.at(jet_it).getPt() > m_jet_list.at(jet_it-1).getPt()
       ) {
      std::cout << "\nWARNING!!! Jets are not pt ordered!\n";
    }
    if (m_jet_list.at(jet_it).getIsBJet()) {
      m_b_jets.push_back(&m_jet_list.at(jet_it));
    }
  }

  // fill leading b jets list
  m_leading_b_jets.reserve(m_b_jets.size());
  for (size_t jet_it = 0; jet_it != m_b_jets.size() && jet_it != 2; ++jet_it) {
    m_leading_b_jets.push_back(m_b_jets.at(jet_it));
  }

  // calculate met and metrel
  m_met.setMetNoint(MET_Truth_NonInt_etx, MET_Truth_NonInt_ety);
  // TODO calculate met-rel or remove
  // m_met.calculateMetRelNoint( m_daughter_el, m_daughter_mu, m_daughter_jet);
}

// -----------------------------------------------------------------------------
void  BMinusL::Cutflow::print()
{
  size_t num_el  = m_daughter_el.size();
  size_t num_mu  = m_daughter_mu.size();
  size_t num_truth_b_quarks = m_daughter_b_quarks.size();
  size_t num_b_jet = m_b_jets.size();

  std::cout << "========================================"
            << "\nevent number: " << EventNumber
            << "\n\ttotal num el: " << m_el_list.size()
            << "\n\ttotal num mu: " << m_mu_list.size()
            << "\n\ttotal num jet: " << m_jet_list.size()
            << "\n\tnum daughter el: " << num_el
            << "\n\tnum daughter mu: " << num_mu
            << "\n\tnum b quarks: " << num_truth_b_quarks
            << "\n\tnum b jets: " << num_b_jet
            << "\n----------------------------------------"
            << "\n";

  /*
  std::cout << "----------------------------------------"
            << "\nall truth particles"
            << "\n----------------------------------------"
            << "\n";
  for (size_t mc_it = 0; mc_it != m_particle_list.size(); ++mc_it) {
    m_particle_list.at(mc_it).printGeneralInfo();
  }
  std::cout << "----------------------------------------"
            << "\nall electrons"
            << "\n----------------------------------------"
            << "\n";
  for (size_t el_it = 0; el_it != m_el_list.size(); ++el_it) {
    m_el_list.at(el_it).print(this);
  }
  std::cout << "----------------------------------------"
            << "\nall muons"
            << "\n----------------------------------------"
            << "\n";
  for (size_t mu_it = 0; mu_it != m_mu_list.size(); ++mu_it) {
    m_mu_list.at(mu_it).print(this);
  }
  std::cout << "----------------------------------------"
            << "\ndaughter electrons"
            << "\n----------------------------------------"
            << "\n";
  for (size_t el_it = 0; el_it != num_el; ++el_it) {
    m_daughter_el.at(el_it)->print(this);
  }
  std::cout << "----------------------------------------"
            << "\ndaughter muons"
            << "\n----------------------------------------"
            << "\n";
  for (size_t mu_it = 0; mu_it != num_mu; ++mu_it) {
    m_daughter_mu.at(mu_it)->print(this);
  }
  std::cout << "----------------------------------------"
            << "\ndaughter b quarks"
            << "\n----------------------------------------"
            << "\n";
  for (size_t q_it = 0; q_it != num_truth_b_quarks; ++q_it) {
    m_daughter_b_quarks.at(q_it)->printGeneralInfo();
  }
  */
}

// =============================================================================
// = BMinusLStandAlone
// =============================================================================
BMinusL::BMinusLStandAloneHistograms::BMinusLStandAloneHistograms()
{
  // initialize histograms here
  m_h_meff = new TH1D( "meff"
                     , "m_{eff} ; m_{eff} [GeV] ; Entries"
                     , 500 , 0, 5000
                     );
}

// -----------------------------------------------------------------------------
void BMinusL::BMinusLStandAloneHistograms::Fill( const BMinusL::Cutflow* cutflow)
{
  // calculate event variables used to fill histograms
  float m_eff = 0;

  size_t num_b  = cutflow->m_daughter_b_quarks.size();
  for (size_t b_it = 0; b_it != num_b; ++b_it) {
    m_eff += cutflow->m_daughter_b_quarks.at(b_it)->getPt();
  }

  size_t num_el = cutflow->m_daughter_el.size();
  for (size_t el_it = 0; el_it != num_el; ++el_it) {
    m_eff += cutflow->m_daughter_el.at(el_it)->getPt();
  }


  size_t num_mu = cutflow->m_daughter_mu.size();
  for (size_t mu_it = 0; mu_it != num_mu; ++mu_it) {
    m_eff += cutflow->m_daughter_mu.at(mu_it)->getPt();
  }

  // Fill histograms based on this event
  m_h_meff->Fill(m_eff/1.e3);
}

// -----------------------------------------------------------------------------
void BMinusL::BMinusLStandAloneHistograms::write(TFile* f)
{
  // change directory to output file
  f->cd();

  // write histograms to output histogram file
  m_h_meff->Write();
}

