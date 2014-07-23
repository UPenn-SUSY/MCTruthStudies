#include "BMinusLCutflow/include/BMinusLCutflow.h" 
#include <iostream>
#include <math.h>
#include <algorithm>

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom.h"

#include "TruthNtupleLooper/include/ObjectDefs.h"
#include "TruthNtupleLooper/include/Calculators.h"

#include "TruthRecordHelpers/include/ParentFinder.h"

#include "HistogramHandlers/include/HistogramHandlers.h"
#include "BMinusLCutflow/include/BMinusLHistogramHandlers.h"
//#include "PDFTool_standalone/PDFTool.h"

// -----------------------------------------------------------------------------
BMinusL::Cutflow::Cutflow(TTree* tree, bool isSignal) : TruthNtuple::TruthNtupleLooper(tree)
{
  m_is_signal = isSignal;

  // construct histogram list
  m_histograms.push_back(new HistogramHandlers::FlavorChannel());
  m_histograms.push_back(new HistogramHandlers::ObjectMultiplicity());
  m_histograms.push_back(new HistogramHandlers::LeptonKinematics());
  m_histograms.push_back(new HistogramHandlers::JetKinematics());
  m_histograms.push_back(new HistogramHandlers::Met());
  m_histograms.push_back(new HistogramHandlers::Mll());
  m_histograms.push_back(new HistogramHandlers::Ht());
  m_histograms.push_back(new HistogramHandlers::Dr());
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

  double m_event_weight = 1./nentries;
  if (m_is_signal == false) {
    // only need to scale bkgd samples
    m_event_weight = scalePDF();
  }
  std::cout<<"About to fill normal histograms.\n";
  // Fill "normal histograms"
  size_t num_hists = m_histograms.size();
  std::cout<<"num_hists = "<<num_hists<<"\n";
  for (size_t hist_it = 0; hist_it != num_hists; ++hist_it) {
    m_histograms.at(hist_it)->Fill( m_flavor_channel
                                  , m_daughter_el
                                  , m_daughter_mu
                                  , m_leading_b_jets
                                  , m_daughter_b_quarks
                                  , m_met
                                  , m_event_weight
                                  );
    std::cout<<"Filling the "<<hist_it<<"th histogram\n";
  }
  std::cout<<"About to fill special histograms.\n";
  // fill special histograms
  m_h_mbl->FillSpecial( m_flavor_channel
                      , m_daughter_el
                      , m_daughter_mu
                      , m_leading_b_jets
		      , m_event_weight
                      );
  std::cout<<"Filled mbl histos.\n";
  m_h_bl_pair_kinematics->FillSpecial( m_flavor_channel
                                     , m_daughter_el
                                     , m_daughter_mu
                                     , m_daughter_b_quarks
				       , m_event_weight
                                     );
  std::cout<<"Filled blpair histos.\n";
  m_h_quark_kinematics->FillSpecial( m_flavor_channel
                                   , m_daughter_b_quarks
				     , m_event_weight
                                   );
  std::cout<<"Filled quark histos.\n";
  m_h_stop_kinematics->FillSpecial( m_flavor_channel
                                  , m_truth_stops
				    , m_event_weight
                                  );
  std::cout<<"Filled stop histos.\n";
  m_sa_hists->Fill(this);
  std::cout<<"Fill the special histograms successfully.\n";
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
     // std::cout << "\nfound a truth electron!\n";
     // std::cout << "\tstatus code: "  << m_truth_electrons.at(el_it)->getStatus() << "\n";
     // std::cout << "\tparent pdgid: " << m_truth_electrons.at(el_it)->getParentPdgid() << "\n";

    if (  (  m_truth_electrons.at(el_it)->getStatus() == 3
          || m_truth_electrons.at(el_it)->getStatus() == 11 // herwig++
          )
	  && ( fabs(m_truth_electrons.at(el_it)->getParentPdgid()) == 1e6+6
	       || m_is_signal == false 
	       )
       ) {
      m_daughter_el.push_back(m_truth_electrons.at(el_it));
    }
    else if (  (  m_truth_electrons.at(el_it)->getStatus() == 1
               || m_truth_electrons.at(el_it)->getStatus() == 3
               || true
               )
	       && fabs(m_truth_electrons.at(el_it)->getParentPdgid()) == 15
            ) {
      // std::cout << "\nthis leptons has a tau parent -- checking the tau's parents\n";

      if (isLeptonFromTauFromStop(m_truth_electrons.at(el_it))) {
        m_daughter_el.push_back(m_truth_electrons.at(el_it));
      }
    }
  }

  // pick muons coming from a stop
  m_daughter_mu.reserve(m_truth_muons.size());
  for (size_t mu_it = 0; mu_it != m_truth_muons.size(); ++mu_it) {

     // std::cout << "\nfound a truth muon!\n";
     // std::cout << "\tstatus code: " << m_truth_muons.at(mu_it)->getStatus() << "\n";
     // std::cout << "\tparent pdgid: " << m_truth_muons.at(mu_it)->getParentPdgid() << "\n";

    if (  (  m_truth_muons.at(mu_it)->getStatus() == 3
          || m_truth_muons.at(mu_it)->getStatus() == 11 // herwig++
          )
	  && ( fabs(m_truth_muons.at(mu_it)->getParentPdgid()) == 1e6+6 
	       || m_is_signal == false
	       )
       ) {
      m_daughter_mu.push_back(m_truth_muons.at(mu_it));
    }
    else if (  (  m_truth_muons.at(mu_it)->getStatus() == 1
               || m_truth_muons.at(mu_it)->getStatus() == 3
               || true
               )
	       && fabs(m_truth_muons.at(mu_it)->getParentPdgid()) == 15
            ) {
      //       std::cout << "\nthis leptons has a tau parent -- checking the tau's parents\n";

      if (isLeptonFromTauFromStop(m_truth_muons.at(mu_it))) {
        m_daughter_mu.push_back(m_truth_muons.at(mu_it));
      }
    }
  }

  // std::cout << "total number el: " << m_daughter_el.size() << "\n";
  // std::cout << "total number mu: " << m_daughter_mu.size() << "\n";
  // std::cout << "total number leptons: " << m_daughter_el.size() + m_daughter_mu.size() << "\n";

  // pick b quarks coming from a stop
  m_daughter_b_quarks.reserve(m_truth_b_quarks.size());
  for (size_t b_quarks_it = 0; b_quarks_it != m_truth_b_quarks.size(); ++b_quarks_it) {
    if ( fabs(m_truth_b_quarks.at(b_quarks_it)->getParentPdgid()) == 1e6+6
	 || m_is_signal == false
	 ) {
      m_daughter_b_quarks.push_back(m_truth_b_quarks.at(b_quarks_it));
    }
  }

  // artificially broaden bquark energy distribution
  broadenResolution();

  size_t num_stops_before_cleaning = m_truth_stops.size();
  // cleanParticleList(m_truth_stops);
  cleanParticleList(m_truth_stops, false);
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

// -----------------------------------------------------------------------------
bool BMinusL::Cutflow::isLeptonFromTauFromStop(const TruthNtuple::Particle* part)
{
  // get the tau index using the barcode
  int tau_index = TruthRecordHelpers::getParentIndexFromBarcode( part->getBarcode()
                                                                , mc_barcode
                                                                , mc_pdgId
                                                                , mc_parent_index
                                                                );
  // find the tau's parent
  int tau_parent_pdgid = TruthRecordHelpers::getParentPdgId( tau_index
                                                           , mc_pdgId
                                                           , mc_parent_index
                                                           );

  // if the parent is 0, we probably lost part of the truth record -- do dR matching to find the correct parent 
  if ( tau_parent_pdgid == 0 ) {
    // do dR matching to find the first tau in this chain
    int first_tau_mc_index = TruthRecordHelpers::doDrMatchForParent( tau_index
                                                                   , mc_pdgId
                                                                   , mc_status
                                                                   , mc_eta
                                                                   , mc_phi
                                                                   , 3
                                                                   , 0.1
                                                                   );
    tau_parent_pdgid = TruthRecordHelpers::getParentPdgId( first_tau_mc_index
                                                         , mc_pdgId
                                                         , mc_parent_index
                                                         );
  }

  // check if taur parent is a tau
  if ( fabs(tau_parent_pdgid) == 1e6+6 ) {
    return true;
  }

  // the parent is not a stop
  return false;
}

// =============================================================================
// = scalePDF
// =============================================================================
double BMinusL::Cutflow::scalePDF()
{
  double m_event_weight = 1.;
  PDFTool* pdfTool = new PDFTool(7000000, 1, -1, 10042);
  pdfTool->setEventInfo( pow(mcevt_pdf_scale->at(0),2)
			,mcevt_pdf_x1->at(0)
			,mcevt_pdf_x2->at(0)
			,mcevt_pdf_id1->at(0)
			,mcevt_pdf_id2->at(0)
			);
  double tmp_pdf = pdfTool->weight(-1,14./13.);     if(tmp_pdf > 500.) tmp_pdf = 0.;
  m_event_weight *= tmp_pdf;
  return m_event_weight;

  // initialize histograms here
  m_h_meff = new TH1D( "meff"
                     , "m_{eff} ; m_{eff} [GeV] ; Entries"
                     , 500 , 0, 5000
                     );

  m_h_pt_b1vsl1 = new TH2D( "pt_b1vsl1"
                          , "p_{t} of subleading bl ; p_{t}^{l} [GeV] ; p_{t}^{b} [GeV]"
                          , 150, 0, 1500
                          , 150, 0, 1500
                          );

  m_h_pt_b1vse1 = new TH2D( "pt_b1vse1"
                          , "p_{t} of subleading be ; p_{t}^{e} [GeV] ; p_{t}^{b} [GeV]"
                          , 150, 0, 1500
                          , 150, 0, 1500
                          ); // for ee or me events

  m_h_pt_b1vsm1 = new TH2D( "pt_b1vsm1"
                          , "p_{t} of subleading bm ; p_{t}^{m} [GeV] ; p_{t}^{b} [GeV]"
                          , 150, 0, 1500
                          , 150, 0, 1500
                          ); // for mm or em events
}

// =============================================================================
// = broadenResolution
// =============================================================================
void BMinusL::Cutflow::broadenResolution()
{
  double b_quark_resolution = 0.3; // THIS IS A GUESS--???

  for (size_t b_quark_it = 0; b_quark_it != m_daughter_b_quarks.size(); b_quark_it++) {
    double rand_scale = gRandom->Gaus(); //default: mean=0., sigma=1.

    double mass  = m_daughter_b_quarks.at(b_quark_it)->getM();
    double theta = m_daughter_b_quarks.at(b_quark_it)->getTheta();
    double phi   = m_daughter_b_quarks.at(b_quark_it)->getPhi();

    double old_e  = m_daughter_b_quarks.at(b_quark_it)->getE();
    double new_e  = old_e * (1. + rand_scale*b_quark_resolution);
    double new_p  = TruthNtuple::pFromEM(new_e, mass);
    double new_pt = TruthNtuple::ptFromPTheta(new_p, theta);
    double new_px = TruthNtuple::pxFromPtPhi(new_pt, phi);
    double new_py = TruthNtuple::pyFromPtPhi(new_pt, phi);
    double new_pz = TruthNtuple::pzFromPtTheta(new_pt, theta);

    m_daughter_b_quarks.at(b_quark_it)->setE(new_e);
    m_daughter_b_quarks.at(b_quark_it)->setPt(new_pt);
    m_daughter_b_quarks.at(b_quark_it)->setPx(new_px);
    m_daughter_b_quarks.at(b_quark_it)->setPy(new_py);
    m_daughter_b_quarks.at(b_quark_it)->setPz(new_pz);
  }
}

// =============================================================================
// = BMinusLStandAlone
// =============================================================================
 BMinusL::BMinusLStandAloneHistograms::BMinusLStandAloneHistograms()
 {
//   // initialize histograms here
//   m_h_meff = new TH1D( "meff"
//                      , "m_{eff} ; m_{eff} [GeV] ; Entries"
//                      , 500 , 0, 5000
//                      );

//   m_h_fc_all__pt_event_b1vl1 = new TH2D( "fc_all__pt_event_b1vl1"
// 		          , "p_{t} of subleading bl ; p_{t}^{l} [GeV] ; p_{t}^{b} [GeV]"
// 			    , 150, 0, 1500
// 			    , 150, 0, 1500
// 			    );

//   m_h_fc_all__pt_event_b1ve1 = new TH2D( "fc_all__pt_event_b1ve1"
// 			    , "p_{t} of subleading be ; p_{t}^{e} [GeV] ; p_{t}^{b} [GeV]"
// 			    , 150, 0, 1500
// 			    , 150, 0, 1500
// 			    ); // for ee or me events

//   m_h_fc_all__pt_event_b1vm1 = new TH2D( "fc_all__pt_event_b1vm1"
// 			    , "p_{t} of subleading bm ; p_{t}^{m} [GeV] ; p_{t}^{b} [GeV]"
// 			    , 150, 0, 1500
//      			    , 150, 0, 1500
// 			    ); // for mm or em events

//   m_h_fc_all__eta_event_eff = new TH1D( "fc_all__eta_event_eff"
// 				, "pass (1) or fail (0) of eta event cuts ; Events ; Passability "
// 				, 2, 0, 2
// 				);
//   m_h_fc_ee__eta_event_eff = new TH1D( "fc_ee__eta_event_eff"
// 				, "pass (1) or fail (0) of eta event cuts ; Events ; Passability "
// 				, 2, 0, 2
// 				);
//   m_h_fc_em__eta_event_eff = new TH1D( "fc_em__eta_event_eff"
// 				, "pass (1) or fail (0) of eta event cuts ; Events ; Passability "
// 				, 2, 0, 2
// 				);
//   m_h_fc_me__eta_event_eff = new TH1D( "fc_me__eta_event_eff"
// 				, "pass (1) or fail (0) of eta event cuts ; Events ; Passability "
// 				, 2, 0, 2
// 				);
//   m_h_fc_mm__eta_event_eff = new TH1D( "fc_mm__eta_event_eff"
// 				, "pass (1) or fail (0) of eta event cuts ; Events ; Passability "
// 				, 2, 0, 2
// 				);

//   m_h_fc_all__lep_eta_event_eff = new TH1D( "fc_all__lep_eta_event_eff"
// 					    , "pass (1) or fail (0) of eta event cuts on leps ; Events ; Passability "
// 					    , 2, 0, 2
// 					    );
//   m_h_fc_all__quark_eta_event_eff = new TH1D( "fc_all__quark_eta_event_eff"
// 					    , "pass (1) or fail (0) of eta event cuts on quarks ; Events ; Passability "
// 					    , 2, 0, 2
// 					    );

//   m_h_fc_all__lep_fiducial_eventall_pass = new TH1D( "fc_all__lep_fiducial_eventall_pass"
// 				, "p_{t} of subleading leps for events in which every object passes #eta<2.4 cuts; p_{t}^{l} [GeV]; Events"
// 				, 150, 0, 1500
// 				);
//   m_h_fc_all__lep_fiducial_eventall_fail = new TH1D( "fc_all__lep_fiducial_eventall_fail"
// 				, "p_{t} of subleading leps for events in which at least one object fails #eta<2.4 cuts; p_{t}^{l} [GeV]; Events"
// 				, 150, 0, 1500
// 				);
//   m_h_fc_all__lep_fiducial_eventlep_pass = new TH1D( "fc_all__lep_fiducial_eventlep_pass"
// 				, "p_{t} of subleading leps for events in which all leps pass #eta<2.4 cuts; p_{t}^{l} [GeV]; Events"
// 				, 150, 0, 1500
// 				);
//   m_h_fc_all__lep_fiducial_eventlep_fail = new TH1D( "fc_all__lep_fiducial_eventlep_fail"
// 				, "p_{t} of subleading leps for events in which at least one lep fails #eta<2.4 cuts; p_{t}^{l} [GeV]; Events"
// 				, 150, 0, 1500
// 				);

//   m_h_fc_all__lep_fiducial_single_pass = new TH1D( "fc_all__lep_fiducial_single_pass"
// 				  , "p_{t} of all leps which pass #eta<2.4 cuts; p_{t}^{l} [GeV]; Events"
// 				  , 150, 0, 1500
// 				  );
//   m_h_fc_all__lep_fiducial_single_fail = new TH1D( "fc_all__lep_fiducial_single_fail"
// 				  , "p_{t} of all leps which fail #eta<2.4 cuts; p_{t}^{l} [GeV]; Events"
// 				  , 150, 0, 1500
// 				  );

//   m_h_fc_all__quark_fiducial_eventall_pass = new TH1D( "fc_all__quark_fiducial_eventall_pass"
// 				, "p_{t} of subleading bs for events in which every object passes #eta<2.4 cuts; p_{t}^{b} [GeV]; Events"
// 				, 150, 0, 1500
// 				);
//   m_h_fc_all__quark_fiducial_eventall_fail = new TH1D( "fc_all__quark_fiducial_eventall_fail"
// 				, "p_{t} of subleading bs for events in which at least one object fails #eta<2.4 cuts; p_{t}^{b} [GeV]; Events"
// 				, 150, 0, 1500
// 				);
//   m_h_fc_all__quark_fiducial_eventquark_pass = new TH1D( "fc_all__quark_fiducial_eventquark_pass"
// 				, "p_{t} of subleading bs for events in which all quarks pass #eta<2.4 cuts; p_{t}^{b} [GeV]; Events"
// 				, 150, 0, 1500
// 				);
//   m_h_fc_all__quark_fiducial_eventquark_fail = new TH1D( "fc_all__quark_fiducial_eventquark_fail"
// 				, "p_{t} of subleading bs for events in which at least one quark fails #eta<2.4 cuts; p_{t}^{b} [GeV]; Events"
// 				, 150, 0, 1500
// 				);

//   m_h_fc_all__quark_fiducial_single_pass = new TH1D( "fc_all__quark_fiducial_single_pass"
// 				  , "p_{t} of all bs which pass #eta<2.4 cuts; p_{t}^{b} [GeV]; Events"
// 				  , 150, 0, 1500
// 				  );
//   m_h_fc_all__quark_fiducial_single_fail = new TH1D( "fc_all__quark_fiducial_single_fail"
// 				  , "p_{t} of all bs which fail #eta<2.4 cuts; p_{t}^{b} [GeV]; Events"
// 				  , 150, 0, 1500
// 				  );

//   /* m_h_fc_all__fiducial_event_e_pass = new TH1D("fc_all__fiducial_event_e_pass"
// 					       ,"p_{t} of subleading es which pass #eta<2.4 cuts; p_{t}^{e} [GeV]; Events"
// 					       , 150, 0, 1500
// 					       );
//   m_h_fc_all__fiducial_event_e_fail = new TH1D("fc_all__fiducial_event_e_fail"
// 					       ,"p_{t} of subleading es which fail #eta<2.4 cuts; p_{t}^{e} [GeV]; Events"
// 					       , 150, 0, 1500
// 					       );
//   m_h_fc_all__fiducial_event_m_pass = new TH1D("fc_all__fiducial_event_m_pass"
// 					       ,"p_{t} of subleading ms which pass #eta<2.4 cuts; p_{t}^{m} [GeV]; Events"
// 					       , 150, 0, 1500
// 					       );
//   m_h_fc_all__fiducial_event_m_fail = new TH1D("fc_all__fiducial_event_m_fail"
// 					       ,"p_{t} of subleading ms which fail #eta<2.4 cuts; p_{t}^{m} [GeV]; Events"
// 					       , 150, 0, 1500
// 					       );
//   */
//  m_h_fc_all__fiducial_event_b1vl1_pass = new TH2D( "fc_all__fiducial_event_b1vl1_pass"
// 				  , "p_{t} of events which pass #{eta}<2.4 cuts; p_{t}^{l} [GeV]; p_{t}^{b} [GeV]"
// 				, 150, 0, 1500
// 				, 150, 0, 1500
// 				);
//   m_h_fc_all__fiducial_event_b1vl1_fail = new TH2D( "fc_all__fiducial_event_b1vl1_fail"
// 				, "p_{t} of events which fail #{eta}<2.4 cuts; p_{t}^{l} [GeV]; p_{t}^{b} [GeV]"
// 				, 150, 0, 1500
// 				, 150, 0, 1500
// 				);

//   m_h_fc_all__fiducial_event_b1vl1_eff = new TH2D("fc_all__fiducial_event_b1vl1_eff",
// 						  "Efficiency of #eta<2.4 and p_{t} cuts; p_{t}^{l} [GeV]; p_{t}^{b} [GeV]",
// 						  150, 0, 1500,
// 						  150, 0, 1500
// 						  );

//   m_h_fc_all__lep_deltaRq = new TH1D( "fc_all__lep_deltaRq"
// 			    , "deltaR of lep and closer b ; deltaR ; Events"
// 			    , 100, 0, 10
// 			    );
//   m_h_fc_all__lep_deltaRq_l0vl1 = new TH2D( "fc_all__lep_deltaRq_l0vl1"
// 			       ,"deltaR of leading lepton and closer b vs subleading lepton and closer b; deltaR^{1} ; deltaR^{0}"
// 			       , 100, 0, 10
// 			       , 100, 0, 10
// 			       );
//   m_h_fc_all__lep_deltaRsameq0 = new TH1D("fc_all__lep_deltaRsameq0"
// 					 , "How often are leps closer to the leading b? ; Diff (0), Same (1) ; Events"
// 					 ,2,0,2
// 					 );
//   m_h_fc_all__lep_deltaRsameq1 = new TH1D("fc_all__lep_deltaRsameq1"
// 					 , "How often are leps closer to the subleading b? ; Diff (0), Same (1) ; Events"
// 					 ,2,0,2
// 					 );
//   m_h_fc_all__lep_deltaRl = new TH1D("fc_all__lep_deltaRl"
// 				       , "deltaR of leps ; deltaR ; Events"
// 				       ,100, 0, 10
// 				       );
//   m_h_fc_all__quark_deltaRq = new TH1D("fc_all__quark_deltaRq"
// 				       , "deltaR of bs ; deltaR ; Events"
// 				       , 100, 0, 10
// 				       );

//   m_h_fc_all__pt_event_b1vl1_eff = new TH2D( "fc_all__pt_event_b1vl1_eff"
// 				,"Efficiency of p_{t} cuts; pt_{l} [GeV]; pt_{b} [GeV]"
// 				,150, 0, 1500
// 				,150, 0, 1500
// 				);
//   m_h_fc_all__pt_event_b1ve1_eff = new TH2D( "fc_all__pt_event_b1ve1_eff"
// 				,"Efficiency of p_{t} cuts; pt_{e} [GeV]; pt_{b} [GeV]"
// 				,150, 0, 1500
// 				,150, 0, 1500
// 				);
//   m_h_fc_all__pt_event_b1vm1_eff = new TH2D( "fc_all__pt_event_b1vm1_eff"
// 				,"Efficiency of p_{t} cuts; pt_{m} [GeV]; pt_{b} [GeV]"
// 				,150, 0, 1500
// 				,150, 0, 1500
// 				);
 }

// -----------------------------------------------------------------------------
void BMinusL::BMinusLStandAloneHistograms::Fill( const BMinusL::Cutflow* cutflow)
{
//   // calculate event variables used to fill histograms
//   float m_eff = 0;

//   size_t num_b  = cutflow->m_daughter_b_quarks.size();
//   for (size_t b_it = 0; b_it != num_b; ++b_it) {
//     m_eff += cutflow->m_daughter_b_quarks.at(b_it)->getPt();
//   }

//   size_t num_el = cutflow->m_daughter_el.size();
//   for (size_t el_it = 0; el_it != num_el; ++el_it) {
//     m_eff += cutflow->m_daughter_el.at(el_it)->getPt();
//   }


//   size_t num_mu = cutflow->m_daughter_mu.size();
//   for (size_t mu_it = 0; mu_it != num_mu; ++mu_it) {
//     m_eff += cutflow->m_daughter_mu.at(mu_it)->getPt();
//   }

//   // Fill histograms based on this event
//   m_h_meff->Fill(m_eff/1.e3);


//   if (num_b == 2 && (num_el+num_mu) ==2) {

//     // subleading pt_b vs pt_l histos
//     float pt_b0, pt_b1, pt_l0, pt_l1;
//     pt_b0 = cutflow->m_daughter_b_quarks.at(0)->getPt();
//     pt_b1 = cutflow->m_daughter_b_quarks.at(1)->getPt();
//     // eta distribution
//     float etab0, etab1, etal0, etal1;
//     etab0 = abs(cutflow->m_daughter_b_quarks.at(0)->getEta());
//     etab1 = abs(cutflow->m_daughter_b_quarks.at(1)->getEta());

//     // deltaR distribution
//     double deltaRq00,deltaRq01,deltaRq10,deltaRq11,deltaRq0,deltaRq1,deltaRq,deltaRl;
//     deltaRq =  TruthNtuple::deltaR(cutflow->m_daughter_b_quarks.at(0),cutflow->m_daughter_b_quarks.at(1));
//     m_h_fc_all__quark_deltaRq->Fill(deltaRq);

//     if (num_el ==2) {
//       pt_l0 = cutflow->m_daughter_el.at(0)->getPt();
//       pt_l1 = cutflow->m_daughter_el.at(1)->getPt();
//       m_h_fc_all__pt_event_b1ve1->Fill(pt_l1/1.e3,pt_b1/1.e3);

//       etal0 = abs(cutflow->m_daughter_el.at(0)->getEta());
//       etal1 = abs(cutflow->m_daughter_el.at(1)->getEta());
//       if (etab0 < 2.4 && etab1 < 2.4 && etal0 < 2.4 && etal1 < 2.4) {
// 	m_h_fc_ee__eta_event_eff->Fill(1);
//       }
//       else {
// 	m_h_fc_ee__eta_event_eff->Fill(0);
//       }

//       deltaRCalc(cutflow->m_daughter_b_quarks.at(0)
// 		 ,cutflow->m_daughter_el.at(0)
// 		 ,cutflow->m_daughter_b_quarks.at(1)
// 		 ,cutflow->m_daughter_el.at(1)
// 		 ,deltaRq00
// 		 ,deltaRq01
// 		 ,deltaRq10
// 		 ,deltaRq11
// 		 );
//       deltaRl = TruthNtuple::deltaR(cutflow->m_daughter_el.at(0),cutflow->m_daughter_el.at(1));
//     }
//     else if (num_mu ==2) {
//       pt_l0 = cutflow->m_daughter_mu.at(0)->getPt();
//       pt_l1 = cutflow->m_daughter_mu.at(1)->getPt();
//       m_h_fc_all__pt_event_b1vm1->Fill(pt_l1/1.e3,pt_b1/1.e3);

//       etal0 = abs(cutflow->m_daughter_mu.at(0)->getEta());
//       etal1 = abs(cutflow->m_daughter_mu.at(1)->getEta());
//       if (etab0 < 2.4 && etab1 < 2.4 && etal0 < 2.4 && etal1 < 2.4) {
// 	m_h_fc_mm__eta_event_eff->Fill(1);
//       }
//       else {
// 	m_h_fc_mm__eta_event_eff->Fill(0);
//       }

//       deltaRCalc(cutflow->m_daughter_b_quarks.at(0)
// 		 ,cutflow->m_daughter_mu.at(0)
// 		 ,cutflow->m_daughter_b_quarks.at(1)
// 		 ,cutflow->m_daughter_mu.at(1)
// 		 ,deltaRq00
// 		 ,deltaRq01
// 		 ,deltaRq10
// 		 ,deltaRq11
// 		 );
//       deltaRl = TruthNtuple::deltaR(cutflow->m_daughter_mu.at(0),cutflow->m_daughter_mu.at(1));
//     }
//     else {
//       float pt_e = cutflow->m_daughter_el.at(0)->getPt();
//       float pt_m = cutflow->m_daughter_mu.at(0)->getPt();
//       if (pt_e < pt_m) {
// 	pt_l0 = pt_m;
// 	pt_l1 = pt_e;
// 	m_h_fc_all__pt_event_b1ve1->Fill(pt_l1/1.e3,pt_b1/1.e3);
// 	etal1 = abs(cutflow->m_daughter_el.at(0)->getEta());
// 	etal0 = abs(cutflow->m_daughter_mu.at(0)->getEta());
// 	if (etab0 < 2.4 && etab1 < 2.4 && etal0 < 2.4 && etal1 < 2.4) {
// 	  m_h_fc_me__eta_event_eff->Fill(1);
// 	}
// 	else {
// 	  m_h_fc_me__eta_event_eff->Fill(0);
// 	}
//       }
//       else {
// 	pt_l0 = pt_e;
// 	pt_l1 = pt_m;
// 	m_h_fc_all__pt_event_b1vm1->Fill(pt_l1/1.e3,pt_b1/1.e3);
// 	etal0 = abs(cutflow->m_daughter_el.at(0)->getEta());
// 	etal1 = abs(cutflow->m_daughter_mu.at(0)->getEta());
// 	if (etab0 < 2.4 && etab1 < 2.4 && etal0 < 2.4 && etal1 < 2.4) {
// 	  m_h_fc_em__eta_event_eff->Fill(1);
// 	}
// 	else {
// 	  m_h_fc_em__eta_event_eff->Fill(0);
// 	}
//       }
//       deltaRCalc(cutflow->m_daughter_b_quarks.at(0)
// 		 ,cutflow->m_daughter_el.at(0)
// 		 ,cutflow->m_daughter_b_quarks.at(1)
// 		 ,cutflow->m_daughter_mu.at(0)
// 		 ,deltaRq00
// 		 ,deltaRq01
// 		 ,deltaRq10
// 		 ,deltaRq11
// 		 );
//       deltaRl = TruthNtuple::deltaR(cutflow->m_daughter_el.at(0),cutflow->m_daughter_mu.at(0));
//     }

//     m_h_fc_all__pt_event_b1vl1->Fill(pt_l1/1.e3,pt_b1/1.e3);

//     // eta event histos (binary: pass(1)/fail(0) ... we get efficiency in python)
//     // and fiducial event histos
//     if (etab0 < 2.4 && etab1 < 2.4 && etal0 < 2.4 && etal1 < 2.4) {
//       m_h_fc_all__eta_event_eff->Fill(1);
//       m_h_fc_all__lep_fiducial_eventall_pass->Fill(pt_l1/1.e3);
//       m_h_fc_all__quark_fiducial_eventall_pass->Fill(pt_b1/1.e3);
//       m_h_fc_all__fiducial_event_b1vl1_pass->Fill(pt_l1/1.e3,pt_b1/1.e3);
//     }
//     else {
//       m_h_fc_all__eta_event_eff->Fill(0);
//       m_h_fc_all__lep_fiducial_eventall_fail->Fill(pt_l1/1.e3);
//       m_h_fc_all__quark_fiducial_eventall_fail->Fill(pt_b1/1.e3);
//       m_h_fc_all__fiducial_event_b1vl1_fail->Fill(pt_l1/1.e3,pt_b1/1.e3);
//     }
//     if (etab0 < 2.4 && etab1 < 2.4) {
//       m_h_fc_all__quark_eta_event_eff->Fill(1);
//       m_h_fc_all__quark_fiducial_eventquark_pass->Fill(pt_b1/1.e3);
//     }
//     else {
//       m_h_fc_all__quark_eta_event_eff->Fill(0);
//       m_h_fc_all__quark_fiducial_eventquark_fail->Fill(pt_b1/1.e3);
//     }
//     if (etal0 < 2.4 && etal1 < 2.4) {
//       m_h_fc_all__lep_eta_event_eff->Fill(1);
//       m_h_fc_all__lep_fiducial_eventlep_pass->Fill(pt_l1/1.e3);
//     }
//     else {
//       m_h_fc_all__lep_eta_event_eff->Fill(0);
//       m_h_fc_all__lep_fiducial_eventlep_fail->Fill(pt_l1/1.e3);
//     }

//     // fiducial single object cut histos
//     if (etab0 < 2.4) {
//       m_h_fc_all__quark_fiducial_single_pass->Fill(pt_b0/1.e3);
//     }
//     else {
//       m_h_fc_all__quark_fiducial_single_fail->Fill(pt_b0/1.e3);
//     }
//     if (etab1 < 2.4) {
//       m_h_fc_all__quark_fiducial_single_pass->Fill(pt_b1/1.e3);
//     }
//     else {
//       m_h_fc_all__quark_fiducial_single_fail->Fill(pt_b1/1.e3);
//     }
//     if (etal0 < 2.4) {
//       m_h_fc_all__lep_fiducial_single_pass->Fill(pt_l0/1.e3);
//     }
//     else {
//       m_h_fc_all__lep_fiducial_single_fail->Fill(pt_l0/1.e3);
//     }
//     if (etal1 < 2.4) {
//       m_h_fc_all__lep_fiducial_single_pass->Fill(pt_l1/1.e3);
//     }
//     else {
//       m_h_fc_all__lep_fiducial_single_fail->Fill(pt_l1/1.e3);
//     }
    
//     // deltaR distribution
//     // take the smaller 2 deltaR's 
//     // for a given lep, which b is closer?
//     // what if both leps are closer to a single b? is this interesting?
//     // probably not, but if it happens a lot then maybe
//     deltaRq0 = std::min(deltaRq00,deltaRq10);
//     deltaRq1 = std::min(deltaRq01,deltaRq11);
//     m_h_fc_all__lep_deltaRq->Fill(deltaRq0);
//     m_h_fc_all__lep_deltaRq->Fill(deltaRq1);
//     if (deltaRq0 == deltaRq00 && deltaRq1 == deltaRq01) {
//       m_h_fc_all__lep_deltaRsameq0->Fill(1);
//     }
//     else {
//       m_h_fc_all__lep_deltaRsameq0->Fill(0);
//     }
//     if (deltaRq0 == deltaRq10 && deltaRq1 == deltaRq11) {
//       m_h_fc_all__lep_deltaRsameq1->Fill(1);
//     }
//     else {
//       m_h_fc_all__lep_deltaRsameq1->Fill(0);
//     }
//     m_h_fc_all__lep_deltaRq_l0vl1->Fill(deltaRq1,deltaRq0);
//     m_h_fc_all__lep_deltaRl->Fill(deltaRl);
//   }

//   size_t num_mu = cutflow->m_daughter_mu.size();
//   for (size_t mu_it = 0; mu_it != num_mu; ++mu_it) {
//     m_eff += cutflow->m_daughter_mu.at(mu_it)->getPt();
//   }

//   // Fill histograms based on this event
//   m_h_meff->Fill(m_eff/1.e3);
}

// -----------------------------------------------------------------------------
// void BMinusL::BMinusLStandAloneHistograms::deltaRCalc(const TruthNtuple::Particle* b0
// 						      , const TruthNtuple::Particle* l0
//  						      , const TruthNtuple::Particle* b1
// 						      , const TruthNtuple::Particle* l1
//  						      , double& deltaR00 
//  						      , double& deltaR01
//  						      , double& deltaR10
// 						      , double& deltaR11
//  						      )
//  {
// //   deltaR00 = TruthNtuple::deltaR(b0,l0);
// //   deltaR01 = TruthNtuple::deltaR(b0,l1);
// //   deltaR10 = TruthNtuple::deltaR(b1,l0);
// //   deltaR11 = TruthNtuple::deltaR(b1,l1);


//  }
// -----------------------------------------------------------------------------
void BMinusL::BMinusLStandAloneHistograms::write(TFile* f)
 {
//   // change directory to output file
//   f->cd();

//   // write histograms to output histogram file
//   m_h_meff->Write();
//   m_h_fc_all__pt_event_b1vl1->Write();
//   m_h_fc_all__pt_event_b1ve1->Write();
//   m_h_fc_all__pt_event_b1vm1->Write();
//   m_h_fc_all__eta_event_eff->Write();
//   m_h_fc_ee__eta_event_eff->Write();
//   m_h_fc_em__eta_event_eff->Write();
//   m_h_fc_me__eta_event_eff->Write();
//   m_h_fc_mm__eta_event_eff->Write();
//   m_h_fc_all__lep_eta_event_eff->Write();
//   m_h_fc_all__quark_eta_event_eff->Write();
//   m_h_fc_all__lep_fiducial_eventall_pass->Write();
//   m_h_fc_all__lep_fiducial_eventall_fail->Write();
//   m_h_fc_all__quark_fiducial_eventall_pass->Write();
//   m_h_fc_all__quark_fiducial_eventall_fail->Write();
//   m_h_fc_all__lep_fiducial_eventlep_pass->Write();
//   m_h_fc_all__lep_fiducial_eventlep_fail->Write();
//   m_h_fc_all__quark_fiducial_eventquark_pass->Write();
//   m_h_fc_all__quark_fiducial_eventquark_fail->Write();
//   m_h_fc_all__lep_fiducial_single_pass->Write();
//   m_h_fc_all__lep_fiducial_single_fail->Write();
//   m_h_fc_all__quark_fiducial_single_pass->Write();
//   m_h_fc_all__quark_fiducial_single_fail->Write();
//   m_h_fc_all__fiducial_event_b1vl1_pass->Write();
//   m_h_fc_all__fiducial_event_b1vl1_fail->Write();
//   m_h_fc_all__lep_deltaRq->Write();
//   m_h_fc_all__lep_deltaRq_l0vl1->Write();
//   m_h_fc_all__lep_deltaRsameq0->Write();
//   m_h_fc_all__lep_deltaRsameq1->Write();
//   m_h_fc_all__lep_deltaRl->Write();
//   m_h_fc_all__quark_deltaRq->Write();

//   // This really shouldn't be done here, but...
//   // calculate efficiency of pt cuts from 2d pt histos

//   // efficiency of b1vsl1:
//   m_h_fc_all__pt_event_b1vl1_eff = calcEffPt(m_h_fc_all__pt_event_b1vl1, "l");
//   m_h_fc_all__pt_event_b1vl1_eff->Write();
//   m_h_fc_all__pt_event_b1ve1_eff = calcEffPt(m_h_fc_all__pt_event_b1ve1, "e");
//   m_h_fc_all__pt_event_b1ve1_eff->Write();
//   m_h_fc_all__pt_event_b1vm1_eff = calcEffPt(m_h_fc_all__pt_event_b1vm1, "m");
//   m_h_fc_all__pt_event_b1vm1_eff->Write();
//   m_h_fc_all__fiducial_event_b1vl1_eff = calcEffFid(m_h_fc_all__fiducial_event_b1vl1_pass
// 						    ,m_h_fc_all__fiducial_event_b1vl1_fail
// 						    );
//   m_h_fc_all__fiducial_event_b1vl1_eff->Write();
 }

// TH2D* BMinusL::BMinusLStandAloneHistograms::calcEffPt(TH2D* h, std::string tag)
// {
//   std::string histo_name = "fc_all__pt_event_b1v";
//   histo_name += tag;
//   histo_name += "1_eff";
//   TH2D* h_eff = new TH2D((histo_name).c_str()
// 			 , " ; p_{t}^{b} [GeV}; p_{t}^{l} [GeV]"
// 			 , 150, 0, 1500
// 			 , 150, 0, 1500
// 			 );
  
//   float denom = h->Integral();
//   for (int ix=0; ix!=h->GetXaxis()->GetNbins(); ix++) {
//     for (int iy=0; iy!=h->GetXaxis()->GetNbins(); iy++) {
//       float cutvaluex = h->GetXaxis()->GetBinLowEdge(ix+1);
//       float cutvaluey = h->GetYaxis()->GetBinLowEdge(iy+1);
//       float numer = h->Integral(ix+1, h->GetXaxis()->GetNbins(), iy+1, h->GetYaxis()->GetNbins());
//       float efficiency = numer/denom;
//       h_eff->Fill(cutvaluex,cutvaluey,efficiency);
//     }
//   }
//   return h_eff;
// }

// TH2D* BMinusL::BMinusLStandAloneHistograms::calcEffFid(TH2D* h_pass, TH2D* h_fail)
// {
//   std::string histo_name = "fc_all__fiducial_event_b1vl1_eff";
//   TH2D* h_eff = new TH2D((histo_name).c_str()
// 			 , " ; p_{t}^{b} [GeV}; p_{t}^{l} [GeV]"
// 			 , 150, 0, 1500
// 			 , 150, 0, 1500
// 			 );
  
//   float denom = h_pass->Integral() + h_fail->Integral();
//   for (int ix=0; ix!=h_pass->GetXaxis()->GetNbins(); ix++) {
//     for (int iy=0; iy!=h_pass->GetXaxis()->GetNbins(); iy++) {
//       float cutvaluex = h_pass->GetXaxis()->GetBinLowEdge(ix+1);
//       float cutvaluey = h_pass->GetYaxis()->GetBinLowEdge(iy+1);
//       float numer = h_pass->Integral(ix+1, h_pass->GetXaxis()->GetNbins(), iy+1, h_pass->GetYaxis()->GetNbins());
//       float efficiency = numer/denom;
//       h_eff->Fill(cutvaluex,cutvaluey,efficiency);
//     }
//   }
//   return h_eff;
// }
