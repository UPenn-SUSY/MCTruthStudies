#include "BMinusLCutflow/include/BMinusLHistogramHandlers.h"
#include "HistogramHandlers/include/HistogramHandlers.h"

#include "TruthNtupleLooper/include/Calculators.h"
#include "TruthNtupleLooper/include/ObjectDefs.h"
#include "TFile.h"

#include <iostream>

// =============================================================================
// = StopKinematics
// =============================================================================
HistogramHandlers::StopKinematics::StopKinematics() : HistogramHandlers::Handle()
{
  const int pt_bins   = 50;
  const double pt_min = 0.;
  const double pt_max = 500.;

  const int eta_bins   = 50;
  const double eta_min = -5.;
  const double eta_max = +5.;

  const int phi_bins   = 64;
  const double phi_min = -3.2;
  const double phi_max = +3.2;

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    m_h_pt_all.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__stop_pt_all"
                                    ).c_str()
                                  , ( "p_{T} - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; p_{T} [GeV] ; Entries"
                                    ).c_str()
                                  , pt_bins, pt_min, pt_max
                                  )
                        );
    m_h_pt_stop.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_pt_stop"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{#tilde{t}} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_astp.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_pt_astp"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{#tilde{t}*} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_pt_diff"
                                     ).c_str()
                                   , ( "p_{T} diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; p_{T}^{#tilde{t}} - p_{T}^{#tilde{t}*} [GeV] ; Entries"
                                     ).c_str()
                                   , 2*pt_bins, -pt_max, pt_max
                                   )
                         );
    m_h_pt_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__stop_pt_2d"
                                   ).c_str()
                                 , ( "p_{T} map - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; p_{T}^{#tilde{t}} [GeV] ; p_{T}^{#tilde{t}*} [GeV]"
                                   ).c_str()
                                 , pt_bins, pt_min, pt_max
                                 , pt_bins, pt_min, pt_max
                                 )
                       );

    m_h_eta_all.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_eta_all"
                                     ).c_str()
                                   , ( "#eta - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #eta ; Entries"
                                     ).c_str()
                                   , eta_bins, eta_min, eta_max
                                   )
                         );
    m_h_eta_stop.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_eta_stop"
                                  ).c_str()
                                , ( "#eta - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; #eta^{#tilde{t}} ; Entries"
                                  ).c_str()
                                , eta_bins, eta_min, eta_max
                                )
                      );
    m_h_eta_astp.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_eta_astp"
                                  ).c_str()
                                , ( "#eta - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; #eta^{#tilde{t}*} ; Entries"
                                  ).c_str()
                                , eta_bins, eta_min, eta_max
                                )
                      );
    m_h_eta_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_eta_diff"
                                     ).c_str()
                                   , ( "#eta diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #eta^{#tilde{t}} - #eta^{#tilde{t}*} ; Entries"
                                     ).c_str()
                                   , 2*eta_bins, -eta_max, eta_max
                                   )
                         );
    m_h_eta_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__stop_eta_2d"
                                   ).c_str()
                                 , ( "#eta map - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #eta^{#tilde{t}} ; #eta^{#tilde{t}*}"
                                   ).c_str()
                                 , eta_bins, eta_min, eta_max
                                 , eta_bins, eta_min, eta_max
                                 )
                       );

    m_h_phi_all.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__stop_phi_all"
                                    ).c_str()
                                  , ( "#phi - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; #phi ; Entries"
                                    ).c_str()
                                  , phi_bins, phi_min, phi_max
                                  )
                        );
    m_h_phi_stop.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_phi_stop"
                                  ).c_str()
                                , ( "#phi - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; #phi^{#tilde{t}} ; Entries"
                                  ).c_str()
                                , phi_bins, phi_min, phi_max
                                )
                      );
    m_h_phi_astp.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_phi_astp"
                                  ).c_str()
                                , ( "#phi - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; #phi^{#tilde{t}*} ; Entries"
                                  ).c_str()
                                , phi_bins, phi_min, phi_max
                                )
                      );
    m_h_phi_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_phi_diff"
                                     ).c_str()
                                   , ( "#phi diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #phi^{#tilde{t}} - #phi^{#tilde{t}*} ; Entries"
                                     ).c_str()
                                   , 2*phi_bins, -phi_max, phi_max
                                   )
                         );
    m_h_phi_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__stop_phi_2d"
                                   ).c_str()
                                 , ( "p_{T} map - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #phi^{#tilde{t}} ; #phi^{#tilde{t}*}"
                                   ).c_str()
                                 , phi_bins, phi_min, phi_max
                                 , phi_bins, phi_min, phi_max
                                 )
                       );
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::StopKinematics::FillSpecial( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Particle*>& stop_list
         )
{
  if (stop_list.size() != 2) {
    std::cout << "ERROR! Number stops (" << stop_list.size() << " != 2\n";
    return;
  }

  double pt_stop = 0;
  double pt_astp = 0;

  double eta_stop = 0;
  double eta_astp = 0;

  double phi_stop = 0;
  double phi_astp = 0;

  if (  stop_list.at(0)->getPdgid() == +(1e6+6)
     && stop_list.at(1)->getPdgid() == -(1e6+6)
     ) {
    pt_stop = stop_list.at(0)->getPt()/1.e3;
    pt_astp = stop_list.at(1)->getPt()/1.e3;

    eta_stop = stop_list.at(0)->getEta();
    eta_astp = stop_list.at(1)->getEta();

    phi_stop = stop_list.at(0)->getPhi();
    phi_astp = stop_list.at(1)->getPhi();
  }
  else if (  stop_list.at(0)->getPdgid() == -(1e6+6)
          && stop_list.at(1)->getPdgid() == +(1e6+6)
          ) {
    pt_stop = stop_list.at(1)->getPt()/1.e3;
    pt_astp = stop_list.at(0)->getPt()/1.e3;

    eta_stop = stop_list.at(1)->getEta();
    eta_astp = stop_list.at(0)->getEta();

    phi_stop = stop_list.at(1)->getPhi();
    phi_astp = stop_list.at(0)->getPhi();
  }
  else {
    std::cout << "ERROR:\n"
              << "\tfirst stop pdgid: " << stop_list.at(0)->getPdgid()
              << "\tsecond stop pdgid: " << stop_list.at(1)->getPdgid()
              << "\n";
    return;
  }

  // fill histograms
  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {
      // fill pt
      m_h_pt_all.at(fc)->Fill(pt_stop);
      m_h_pt_all.at(fc)->Fill(pt_astp);

      m_h_pt_stop.at(fc)->Fill(pt_stop);
      m_h_pt_astp.at(fc)->Fill(pt_astp);

      m_h_pt_diff.at(fc)->Fill(pt_stop - pt_astp);
      m_h_pt_2d.at(fc)->Fill(pt_stop, pt_astp);

      // Fill eta
      m_h_eta_all.at(fc)->Fill(eta_stop);
      m_h_eta_all.at(fc)->Fill(eta_astp);

      m_h_eta_stop.at(fc)->Fill(eta_stop);
      m_h_eta_astp.at(fc)->Fill(eta_astp);

      m_h_eta_diff.at(fc)->Fill(eta_stop - eta_astp);
      m_h_eta_2d.at(fc)->Fill(eta_stop, eta_astp);

      // Fill phi
      m_h_phi_all.at(fc)->Fill(phi_stop);
      m_h_phi_all.at(fc)->Fill(phi_astp);

      m_h_phi_stop.at(fc)->Fill(phi_stop);
      m_h_phi_astp.at(fc)->Fill(phi_astp);

      m_h_phi_diff.at(fc)->Fill(phi_stop - phi_astp);
      m_h_phi_2d.at(fc)->Fill(phi_stop, phi_astp);
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::StopKinematics::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
      m_h_pt_all.at(fc_it)->Write();
      m_h_pt_stop.at(fc_it)->Write();
      m_h_pt_astp.at(fc_it)->Write();
      m_h_pt_diff.at(fc_it)->Write();
      m_h_pt_2d.at(fc_it)->Write();

      m_h_eta_all.at(fc_it)->Write();
      m_h_eta_stop.at(fc_it)->Write();
      m_h_eta_astp.at(fc_it)->Write();
      m_h_eta_diff.at(fc_it)->Write();
      m_h_eta_2d.at(fc_it)->Write();

      m_h_phi_all.at(fc_it)->Write();
      m_h_phi_stop.at(fc_it)->Write();
      m_h_phi_astp.at(fc_it)->Write();
      m_h_phi_diff.at(fc_it)->Write();
      m_h_phi_2d.at(fc_it)->Write();
  }
}


// =============================================================================
// = Mbl
// =============================================================================
HistogramHandlers::Mbl::Mbl() : HistogramHandlers::Handle()
{
  const int bins   = 50;
  const double min = 0;
  const double max = 500;

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    m_h_mbl_truth.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__mbl_truth"
                                       ).c_str()
                                     , ( "m_{bl}^{truth} - "
                                       + TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "; m_{bl}^{truth} [GeV] ; Entries"
                                       ).c_str()
                                     , bins, min, max
                                     )
                           );
    m_h_mbl_dphi_matching.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__mbl_dphi_matching"
                                               ).c_str()
                                             , ( "m_{bl} - "
                                               + TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "; m_{bl} [GeV] ; Entries"
                                               ).c_str()
                                             , bins, min, max
                                             )
                                   );
    m_h_mbl_dr_matching.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                             + "__mbl_dr_matching"
                                             ).c_str()
                                           , ( "m_{bl} - "
                                             + TruthNtuple::FlavorChannelStrings[fc_it]
                                             + "; m_{bl} [GeV] ; Entries"
                                             ).c_str()
                                           , bins, min, max
                                           )
                                 );
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::Mbl::FillSpecial( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Electron*>& el_list
         , const std::vector<TruthNtuple::Muon*>& mu_list
         , const std::vector<TruthNtuple::Particle*>& quark_list
         )
{
  if (flavor_channel == TruthNtuple::FLAVOR_NONE) return;

  // merge el and mu lists to lepton list
  std::vector<TruthNtuple::Lepton*> lep_list;
  lep_list.reserve(el_list.size() + mu_list.size());
  lep_list.insert(lep_list.end(), el_list.begin(), el_list.end());
  lep_list.insert(lep_list.end(), mu_list.begin(), mu_list.end());

  // get lists of mbl values for different matching methods
  std::vector<double> mbl_list_truth         = TruthNtuple::getInvariantMassList(lep_list, quark_list, 0);
  std::vector<double> mbl_list_dphi_matching = TruthNtuple::getInvariantMassList(lep_list, quark_list, 1);
  std::vector<double> mbl_list_dr_matching   = TruthNtuple::getInvariantMassList(lep_list, quark_list, 2);

  // loop over flavor channels
  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {
      // fill mbl truth histogram
      for (size_t mbl_it = 0; mbl_it != mbl_list_truth.size(); ++mbl_it) {
        m_h_mbl_truth.at(fc)->Fill(mbl_list_truth.at(mbl_it)/1.e3);
      }
      // fill mbl dphi matching histogram
      for (size_t mbl_it = 0; mbl_it != mbl_list_dphi_matching.size(); ++mbl_it) {
        m_h_mbl_dphi_matching.at(fc)->Fill(mbl_list_dphi_matching.at(mbl_it)/1.e3);
      }
      // fill mbl dr matching histogram
      for (size_t mbl_it = 0; mbl_it != mbl_list_dr_matching.size(); ++mbl_it) {
        m_h_mbl_dr_matching.at(fc)->Fill(mbl_list_dr_matching.at(mbl_it)/1.e3);
      }
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::Mbl::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
      m_h_mbl_truth.at(fc_it)->Write();
      m_h_mbl_dphi_matching.at(fc_it)->Write();
      m_h_mbl_dr_matching.at(fc_it)->Write();
  }
}

