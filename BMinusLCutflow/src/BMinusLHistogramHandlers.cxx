#include "BMinusLCutflow/include/BMinusLHistogramHandlers.h"
#include "HistogramHandlers/include/HistogramHandlers.h"

#include "TruthNtupleLooper/include/Calculators.h"
#include "TruthNtupleLooper/include/ObjectDefs.h"
#include "TFile.h"

#include <iostream>
#include <math.h>
#include <algorithm>

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
// = QuarkKinematics
// =============================================================================
HistogramHandlers::QuarkKinematics::QuarkKinematics() : HistogramHandlers::Handle()
{
  const int    pt_bins = 50;
  const double pt_min  = 0.;
  const double pt_max  = 500.;

  const int    eta_bins = 50;
  const double eta_min  = -5.;
  const double eta_max  = +5.;

  const int    phi_bins = 64;
  const double phi_min  = -3.2;
  const double phi_max  = +3.2;

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // initialize pt histograms
    m_h_pt_all.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__quark_pt_all"
                                    ).c_str()
                                  , ( "p_{T} - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; p_{T} [GeV] ; Entries"
                                    ).c_str()
                                  , pt_bins, pt_min, pt_max
                                  )
                        );
    m_h_pt_0.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__quark_pt_0"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{0} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_1.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__quark_pt_1"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{1} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__quark_pt_diff"
                                     ).c_str()
                                   , ( "p_{T} diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; p_{T}^{0} - p_{T}^{1} [GeV] ; Entries"
                                     ).c_str()
                                   , pt_bins, pt_min, pt_max
                                   )
                         );
    m_h_pt_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__quark_pt_2d"
                                   ).c_str()
                                 , ( "p_{T} map - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; p_{T}^{0} [GeV] ; p_{T}^{1} [GeV]"
                                   ).c_str()
                                 , pt_bins, pt_min, pt_max
                                 , pt_bins, pt_min, pt_max
                                 )
                       );


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // initialize eta histograms
    m_h_eta_all.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__quark_eta_all"
                                     ).c_str()
                                   , ( "#eta - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #eta ; Entries"
                                     ).c_str()
                                   , eta_bins, eta_min, eta_max
                                   )
                         );
    m_h_eta_0.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__quark_eta_0"
                                   ).c_str()
                                 , ( "#eta - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #eta^{0} ; Entries"
                                   ).c_str()
                                 , eta_bins, eta_min, eta_max
                                 )
                       );
    m_h_eta_1.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__quark_eta_1"
                                   ).c_str()
                                 , ( "#eta - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #eta^{1} ; Entries"
                                   ).c_str()
                                 , eta_bins, eta_min, eta_max
                                 )
                       );
    m_h_eta_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__quark_eta_diff"
                                      ).c_str()
                                    , ( "#eta diff - "
                                      + TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "; #eta^{0} - #eta^{1} ; Entries"
                                      ).c_str()
                                    , eta_bins/2, 0, eta_max
                                    )
                          );
    m_h_eta_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__quark_eta_2d"
                                    ).c_str()
                                  , ( "#eta map - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; #eta^{0} ; #eta^{1}"
                                    ).c_str()
                                  , eta_bins, eta_min, eta_max
                                  , eta_bins, eta_min, eta_max
                                  )
                        );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // initialize phi histograms
    m_h_phi_all.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__quark_phi_all"
                                     ).c_str()
                                   , ( "#phi - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #phi ; Entries"
                                     ).c_str()
                                   , phi_bins, phi_min, phi_max
                                   )
                         );
    m_h_phi_0.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__quark_phi_0"
                                   ).c_str()
                                 , ( "#phi - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #phi^{0} ; Entries"
                                   ).c_str()
                                 , phi_bins, phi_min, phi_max
                                 )
                       );
    m_h_phi_1.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__quark_phi_1"
                                   ).c_str()
                                 , ( "#phi - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #phi^{1} ; Entries"
                                   ).c_str()
                                 , phi_bins, phi_min, phi_max
                                 )
                       );
    m_h_phi_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__quark_phi_diff"
                                      ).c_str()
                                    , ( "#phi diff - "
                                      + TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "; #phi^{0} - #phi^{1} ; Entries"
                                      ).c_str()
                                    , phi_bins/2, 0, phi_max
                                    )
                          );
    m_h_phi_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__quark_phi_2d"
                                    ).c_str()
                                  , ( "#phi map - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; #phi^{0} ; #phi^{1}"
                                    ).c_str()
                                  , phi_bins, phi_min, phi_max
                                  , phi_bins, phi_min, phi_max
                                  )
                        );
    
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::QuarkKinematics::FillSpecial( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Particle*>& quark_list
         )
{
  if (quark_list.size() != 2) {
    std::cout << "ERROR! Number quarks (" << quark_list.size() << " != 2\n";
    return;
  }

  double pt_0 = quark_list.at(0)->getPt()/1.e3;
  double pt_1 = quark_list.at(1)->getPt()/1.e3;

  double eta_0 = quark_list.at(0)->getEta();
  double eta_1 = quark_list.at(1)->getEta();

  double phi_0 = quark_list.at(0)->getPhi();
  double phi_1 = quark_list.at(1)->getPhi();

  // check pt ordering
  if (pt_0 < pt_1) {
    std::swap(pt_0, pt_1);
    std::swap(eta_0, eta_1);
    std::swap(phi_0, phi_1);

    // double tmp_pt = pt_1;
    // pt_1 = pt_0;
    // pt_0 = tmp_pt;

    // double tmp_eta = eta_1;
    // eta_1 = eta_0;
    // eta_0 = tmp_eta;

    // double tmp_phi = phi_1;
    // phi_1 = phi_0;
    // phi_0 = tmp_phi;
  }

  // fill histograms
  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {
      // fill pt histograms
      m_h_pt_all.at(fc)->Fill(pt_0);
      m_h_pt_all.at(fc)->Fill(pt_1);

      m_h_pt_0.at(fc)->Fill(pt_0);
      m_h_pt_1.at(fc)->Fill(pt_1);

      m_h_pt_diff.at(fc)->Fill(pt_0 - pt_1);
      m_h_pt_2d.at(fc)->Fill(pt_0, pt_1);

      // fill eta histograms
      m_h_eta_all.at(fc)->Fill(eta_0);
      m_h_eta_all.at(fc)->Fill(eta_1);

      m_h_eta_0.at(fc)->Fill(eta_0);
      m_h_eta_1.at(fc)->Fill(eta_1);

      m_h_eta_diff.at(fc)->Fill(TruthNtuple::deltaEta(eta_0, eta_1));
      m_h_eta_2d.at(fc)->Fill(eta_0, eta_1);

      // fill phi histograms
      m_h_phi_all.at(fc)->Fill(phi_0);
      m_h_phi_all.at(fc)->Fill(phi_1);

      m_h_phi_0.at(fc)->Fill(phi_0);
      m_h_phi_1.at(fc)->Fill(phi_1);

      m_h_phi_diff.at(fc)->Fill(TruthNtuple::deltaPhi(phi_0, phi_1));
      m_h_phi_2d.at(fc)->Fill(phi_0, phi_1);
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::QuarkKinematics::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
      m_h_pt_all.at(fc_it)->Write();
      m_h_pt_0.at(fc_it)->Write();
      m_h_pt_1.at(fc_it)->Write();
      m_h_pt_diff.at(fc_it)->Write();
      m_h_pt_2d.at(fc_it)->Write();

      m_h_eta_all.at(fc_it)->Write();
      m_h_eta_0.at(fc_it)->Write();
      m_h_eta_1.at(fc_it)->Write();
      m_h_eta_diff.at(fc_it)->Write();
      m_h_eta_2d.at(fc_it)->Write();

      m_h_phi_all.at(fc_it)->Write();
      m_h_phi_0.at(fc_it)->Write();
      m_h_phi_1.at(fc_it)->Write();
      m_h_phi_diff.at(fc_it)->Write();
      m_h_phi_2d.at(fc_it)->Write();
  }
}



// =============================================================================
// = BLPairKinematics
// =============================================================================
HistogramHandlers::BLPairKinematics::BLPairKinematics() : m_l_from_stop(0)
                                                        , m_l_from_astp(0)
                                                        , m_b_from_stop(0)
                                                        , m_b_from_astp(0)
{
  const int    pt_bins = 50;
  const double pt_min  = 0.;
  const double pt_max  = 500.;

  const int    eta_bins = 50;
  const double eta_min  = -5.;
  const double eta_max  = +5.;

  const int    phi_bins = 64;
  const double phi_min  = -3.2;
  const double phi_max  = +3.2;

  const int    mbl_bins = 50;
  const double mbl_min  = 0.;
  const double mbl_max  = 500.;

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_l_pt_stop.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__l_pt_stop" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - p_{T}" // title suffix
                                       + " ; p_{T} [GeV]" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , pt_bins, pt_min, pt_max
                                     )
                           );

    m_h_l_pt_astp.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__l_pt_astp" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - p_{T}" // title suffix
                                       + " ; p_{T} [GeV]" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , pt_bins, pt_min, pt_max
                                     )
                           );

    m_h_l_eta_stop.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__l_eta_stop" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #eta" // title suffix
                                       + " ; #eta" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , eta_bins, eta_min, eta_max
                                     )
                           );

    m_h_l_eta_astp.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__l_eta_astp" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #eta" // title suffix
                                       + " ; #eta" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , eta_bins, eta_min, eta_max
                                     )
                           );

    m_h_l_phi_stop.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__l_phi_stop" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #phi" // title suffix
                                       + " ; #phi" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , phi_bins, phi_min, phi_max
                                     )
                           );

    m_h_l_phi_astp.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__l_phi_astp" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #phi" // title suffix
                                       + " ; #phi" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , phi_bins, phi_min, phi_max
                                     )
                           );


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_b_pt_stop.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__b_pt_stop" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - p_{T}" // title suffix
                                       + " ; p_{T} [GeV]" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , pt_bins, pt_min, pt_max
                                     )
                           );

    m_h_b_pt_astp.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__b_pt_astp" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - p_{T}" // title suffix
                                       + " ; p_{T} [GeV]" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , pt_bins, pt_min, pt_max
                                     )
                           );

    m_h_b_eta_stop.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__b_eta_stop" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #eta" // title suffix
                                       + " ; #eta" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , eta_bins, eta_min, eta_max
                                     )
                           );

    m_h_b_eta_astp.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__b_eta_astp" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #eta" // title suffix
                                       + " ; #eta" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , eta_bins, eta_min, eta_max
                                     )
                           );

    m_h_b_phi_stop.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__b_phi_stop" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #phi" // title suffix
                                       + " ; #phi" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , phi_bins, phi_min, phi_max
                                     )
                           );

    m_h_b_phi_astp.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__b_phi_astp" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #phi" // title suffix
                                       + " ; #phi" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , phi_bins, phi_min, phi_max
                                     )
                           );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_right_pair_bl_pt_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__right_pair_bl_pt_2d" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - p_{T}" // title suffix
                                                 + " ; p_{T}^{l} [GeV]" // x-axis label
                                                 + " ; p_{T}^{b} [GeV]" // y-axis label
                                                 ).c_str()
                                               , pt_bins, pt_min, pt_max
                                               , pt_bins, pt_min, pt_max
                                               )
                                     );

    m_h_right_pair_bl_eta_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__right_pair_bl_eta_2d" // name suffix
                                                  ).c_str()
                                                , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - #eta" // title suffix
                                                  + " ; #eta^{l}" // x-axis label
                                                  + " ; #eta^{b}" // y-axis label
                                                  ).c_str()
                                                , eta_bins, eta_min, eta_max
                                                , eta_bins, eta_min, eta_max
                                                )
                                      );

    m_h_right_pair_bl_phi_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__right_pair_bl_phi_2d" // name suffix
                                                  ).c_str()
                                                , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - #phi" // title suffix
                                                  + "; #phi^{l}" // x-axis label
                                                  + "; #phi^{b}" // y-axis label
                                                  ).c_str()
                                                , phi_bins, phi_min, phi_max
                                                , phi_bins, phi_min, phi_max
                                                )
                                      );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_pair_bl_pt_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__wrong_pair_bl_pt_2d" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - p_{T}" // title suffix
                                                 + " ; p_{T}^{l} [GeV]" // x-axis label
                                                 + " ; p_{T}^{b} [GeV]" // y-axis label
                                                 ).c_str()
                                               , pt_bins, pt_min, pt_max
                                               , pt_bins, pt_min, pt_max
                                               )
                                     );

    m_h_wrong_pair_bl_eta_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__wrong_pair_bl_eta_2d" // name suffix
                                                  ).c_str()
                                                , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - #eta" // title suffix
                                                  + " ; #eta^{l}" // x-axis label
                                                  + " ; #eta^{b}" // y-axis label
                                                  ).c_str()
                                                , eta_bins, eta_min, eta_max
                                                , eta_bins, eta_min, eta_max
                                                )
                                      );

    m_h_wrong_pair_bl_phi_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__wrong_pair_bl_phi_2d" // name suffix
                                                  ).c_str()
                                                , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - #phi" // title suffix
                                                  + " ; #phi^{l}" // x-axis label
                                                  + " ; #phi^{b}" // y-axis label
                                                  ).c_str()
                                                , phi_bins, phi_min, phi_max
                                                , phi_bins, phi_min, phi_max
                                                )
                                      );


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_right_pair_mbl_all.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__right_pair_mbl_all" // name suffix
                                                ).c_str()
                                              , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + " - mbl" // title suffix
                                                + "; m_{bl} [GeV]" // x-axis label
                                                + "; Entries" // y-axis label
                                                ).c_str()
                                              , mbl_bins, mbl_min, mbl_max
                                              )
                                    );

    m_h_right_pair_mbl_stop.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__right_pair_mbl_stop" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - mbl" // title suffix
                                                 + "; m_{bl}^{#tilde{t}} [GeV]" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , mbl_bins, mbl_min, mbl_max
                                               )
                                     );

    m_h_right_pair_mbl_astp.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__right_pair_mbl_astp" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - mbl" // title suffix
                                                 + "; m_{bl}^{#tilde{t}*} [GeV]" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , mbl_bins, mbl_min, mbl_max
                                               )
                                     );

    m_h_right_pair_mbl_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__right_pair_mbl_diff" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - mbl" // title suffix
                                                 + "; |m_{bl}^{#tilde{t}} - m_{bl}^{#tilde{t}*}| [GeV]" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , mbl_bins, mbl_min, mbl_max
                                               )
                                     );

    m_h_right_pair_mbl_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__right_pair_mbl_2d" // name suffix
                                               ).c_str()
                                             , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + " - mbl" // title suffix
                                               + "; m_{bl}^{#tilde{t}} [GeV]" // x-axis label
                                               + "; #m_{bl}^{#tilde{t}*} [GeV]" // y-axis label
                                               ).c_str()
                                             , mbl_bins, mbl_min, mbl_max
                                             , mbl_bins, mbl_min, mbl_max
                                             )
                                   );


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_pair_mbl_all.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__wrong_pair_mbl_all" // name suffix
                                                ).c_str()
                                              , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + " - mbl" // title suffix
                                                + "; m_{bl} [GeV]" // x-axis label
                                                + "; Entries" // y-axis label
                                                ).c_str()
                                              , mbl_bins, mbl_min, mbl_max
                                              )
                                    );

    m_h_wrong_pair_mbl_0.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + "__wrong_pair_mbl_0" // name suffix
                                              ).c_str()
                                            , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + " - mbl" // title suffix
                                              + "; m_{bl}^{0} [GeV]" // x-axis label
                                              + "; Entries" // y-axis label
                                              ).c_str()
                                            , mbl_bins, mbl_min, mbl_max
                                            )
                                  );

    m_h_wrong_pair_mbl_1.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + "__wrong_pair_mbl_1" // name suffix
                                              ).c_str()
                                            , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + " - mbl" // title suffix
                                              + "; m_{bl}^{1} [GeV]" // x-axis label
                                              + "; Entries" // y-axis label
                                              ).c_str()
                                            , mbl_bins, mbl_min, mbl_max
                                            )
                                  );

    m_h_wrong_pair_mbl_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__wrong_pair_mbl_diff" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - mbl" // title suffix
                                                 + "; |m_{bl}^{0} - m_{bl}^{1}| [GeV]" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , mbl_bins, mbl_min, mbl_max
                                               )
                                     );

    m_h_wrong_pair_mbl_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__wrong_pair_mbl_2d" // name suffix
                                               ).c_str()
                                             , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + " - mbl" // title suffix
                                               + "; m_{bl}^{0} [GeV]" // x-axis label
                                               + "; #m_{bl}^{1} [GeV]" // y-axis label
                                               ).c_str()
                                             , mbl_bins, mbl_min, mbl_max
                                             , mbl_bins, mbl_min, mbl_max
                                             )
                                   );


  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::BLPairKinematics::FillSpecial( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Electron*>& el_list
         , const std::vector<TruthNtuple::Muon*>& mu_list
         , const std::vector<TruthNtuple::Particle*>& quark_list
         )
{
  if (flavor_channel == TruthNtuple::FLAVOR_NONE) return;

  // sort objects based on parent particle - if sorting fails, exis the function
  if (sortObjects(el_list, mu_list, quark_list) == false) return;

  // get and store the mbl for all combinations
  double mbl_stop = TruthNtuple::invariantMass(m_l_from_stop, m_b_from_stop)/1.e3;
  double mbl_astp = TruthNtuple::invariantMass(m_l_from_astp, m_b_from_astp)/1.e3;

  double mbl_wrong_0 = TruthNtuple::invariantMass(m_l_from_stop, m_b_from_astp)/1.e3;
  double mbl_wrong_1 = TruthNtuple::invariantMass(m_l_from_astp, m_b_from_stop)/1.e3;
  if (mbl_wrong_0 < mbl_wrong_1) std::swap(mbl_wrong_0, mbl_wrong_1);


  // loop through flavor channels
  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // fill lepton kinematic plots
      m_h_l_pt_stop.at(fc)->Fill(m_l_from_stop->getPt()/1.e3);
      m_h_l_pt_astp.at(fc)->Fill(m_l_from_astp->getPt()/1.e3);

      m_h_l_eta_stop.at(fc)->Fill(m_l_from_stop->getEta());
      m_h_l_eta_astp.at(fc)->Fill(m_l_from_astp->getEta());

      m_h_l_phi_stop.at(fc)->Fill(m_l_from_stop->getPhi());
      m_h_l_phi_astp.at(fc)->Fill(m_l_from_astp->getPhi());

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // fill b quark kinematic plots
      m_h_b_pt_stop.at(fc)->Fill(m_b_from_stop->getPt()/1.e3);
      m_h_b_pt_astp.at(fc)->Fill(m_b_from_astp->getPt()/1.e3);

      m_h_b_eta_stop.at(fc)->Fill(m_b_from_stop->getEta());
      m_h_b_eta_astp.at(fc)->Fill(m_b_from_astp->getEta());

      m_h_b_phi_stop.at(fc)->Fill(m_b_from_stop->getPhi());
      m_h_b_phi_astp.at(fc)->Fill(m_b_from_astp->getPhi());

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // fill b-vs-l histograms for right pairing
      m_h_right_pair_bl_pt_2d.at(fc)->Fill( m_l_from_stop->getPt()/1.e3
                                          , m_b_from_stop->getPt()/1.e3
                                          );
      m_h_right_pair_bl_eta_2d.at(fc)->Fill( m_l_from_stop->getEta()
                                           , m_b_from_stop->getEta()
                                           );
      m_h_right_pair_bl_phi_2d.at(fc)->Fill( m_l_from_stop->getPhi()
                                           , m_b_from_stop->getPhi()
                                           );

      m_h_right_pair_bl_pt_2d.at(fc)->Fill( m_l_from_astp->getPt()/1.e3
                                          , m_b_from_astp->getPt()/1.e3
                                          );
      m_h_right_pair_bl_eta_2d.at(fc)->Fill( m_l_from_astp->getEta()
                                           , m_b_from_astp->getEta()
                                           );
      m_h_right_pair_bl_phi_2d.at(fc)->Fill( m_l_from_astp->getPhi()
                                           , m_b_from_astp->getPhi()
                                           );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // fill b-vs-l histograms for wrong pairing
      m_h_wrong_pair_bl_pt_2d.at(fc)->Fill( m_l_from_stop->getPt()/1.e3
                                          , m_b_from_astp->getPt()/1.e3
                                          );
      m_h_wrong_pair_bl_eta_2d.at(fc)->Fill( m_l_from_stop->getEta()
                                           , m_b_from_astp->getEta()
                                           );
      m_h_wrong_pair_bl_phi_2d.at(fc)->Fill( m_l_from_stop->getPhi()
                                           , m_b_from_astp->getPhi()
                                           );

      m_h_wrong_pair_bl_pt_2d.at(fc)->Fill( m_l_from_astp->getPt()/1.e3
                                          , m_b_from_stop->getPt()/1.e3
                                          );
      m_h_wrong_pair_bl_eta_2d.at(fc)->Fill( m_l_from_astp->getEta()
                                           , m_b_from_stop->getEta()
                                           );
      m_h_wrong_pair_bl_phi_2d.at(fc)->Fill( m_l_from_astp->getPhi()
                                           , m_b_from_stop->getPhi()
                                           );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // fill mbl histograms for right pairing
      m_h_right_pair_mbl_all.at(fc)->Fill(mbl_stop);
      m_h_right_pair_mbl_all.at(fc)->Fill(mbl_astp);
      m_h_right_pair_mbl_stop.at(fc)->Fill(mbl_stop);
      m_h_right_pair_mbl_astp.at(fc)->Fill(mbl_astp);
      m_h_right_pair_mbl_diff.at(fc)->Fill(fabs(mbl_stop - mbl_astp));
      m_h_right_pair_mbl_2d.at(fc)->Fill(mbl_stop, mbl_astp);

      // fill mbl histograms for wrong pairing
      m_h_wrong_pair_mbl_all.at(fc)->Fill(mbl_wrong_0);
      m_h_wrong_pair_mbl_all.at(fc)->Fill(mbl_wrong_1);
      m_h_wrong_pair_mbl_0.at(fc)->Fill(mbl_wrong_0);
      m_h_wrong_pair_mbl_1.at(fc)->Fill(mbl_wrong_1);
      m_h_wrong_pair_mbl_diff.at(fc)->Fill(fabs(mbl_wrong_0 - mbl_wrong_1));
      m_h_wrong_pair_mbl_2d.at(fc)->Fill(mbl_wrong_0, mbl_wrong_1);
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::BLPairKinematics::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    m_h_l_pt_stop.at(fc_it)->Write();
    m_h_l_pt_astp.at(fc_it)->Write();
    m_h_l_eta_stop.at(fc_it)->Write();
    m_h_l_eta_astp.at(fc_it)->Write();
    m_h_l_phi_stop.at(fc_it)->Write();
    m_h_l_phi_astp.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_b_pt_stop.at(fc_it)->Write();
    m_h_b_pt_astp.at(fc_it)->Write();
    m_h_b_eta_stop.at(fc_it)->Write();
    m_h_b_eta_astp.at(fc_it)->Write();
    m_h_b_phi_stop.at(fc_it)->Write();
    m_h_b_phi_astp.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_right_pair_bl_pt_2d.at(fc_it)->Write();
    m_h_right_pair_bl_eta_2d.at(fc_it)->Write();
    m_h_right_pair_bl_phi_2d.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_pair_bl_pt_2d.at(fc_it)->Write();
    m_h_wrong_pair_bl_eta_2d.at(fc_it)->Write();
    m_h_wrong_pair_bl_phi_2d.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_right_pair_mbl_all.at(fc_it)->Write();
    m_h_right_pair_mbl_stop.at(fc_it)->Write();
    m_h_right_pair_mbl_astp.at(fc_it)->Write();
    m_h_right_pair_mbl_diff.at(fc_it)->Write();
    m_h_right_pair_mbl_2d.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_pair_mbl_all.at(fc_it)->Write();
    m_h_wrong_pair_mbl_0.at(fc_it)->Write();
    m_h_wrong_pair_mbl_1.at(fc_it)->Write();
    m_h_wrong_pair_mbl_diff.at(fc_it)->Write();
    m_h_wrong_pair_mbl_2d.at(fc_it)->Write();
  }
}

// -----------------------------------------------------------------------------
bool HistogramHandlers::BLPairKinematics::sortObjects( const std::vector<TruthNtuple::Electron*>& el_list
                                                     , const std::vector<TruthNtuple::Muon*>& mu_list
                                                     , const std::vector<TruthNtuple::Particle*>& quark_list
                                                     )
{
  // merge el and mu lists to lepton list
  std::vector<TruthNtuple::Lepton*> lep_list;
  lep_list.reserve(el_list.size() + mu_list.size());
  lep_list.insert(lep_list.end(), el_list.begin(), el_list.end());
  lep_list.insert(lep_list.end(), mu_list.begin(), mu_list.end());

  // Find lepton from stop and anti-stop
  m_l_from_stop = 0;
  m_l_from_astp = 0;

  // I'm probably being overly careful here
  for (size_t lep_it = 0; lep_it != lep_list.size(); ++lep_it) {
    // if parent is stop
    if (lep_list.at(lep_it)->getParentPdgid() == +(1e6+6))  {
      if (m_l_from_stop == 0) m_l_from_stop = lep_list.at(lep_it);
      else
        std::cout << "WARNING! Found multiple leptons paired to stop!\n";
    }
    // if parent is anti-stop
    if (lep_list.at(lep_it)->getParentPdgid() == -(1e6+6)) {
      if (m_l_from_astp == 0) m_l_from_astp = lep_list.at(lep_it);
      else
        std::cout << "WARNING! Found multiple leptons paired to anti-stop!\n";
    }
  }
  if (m_l_from_stop == 0 || m_l_from_astp == 0) {
    std::cout << "ERROR!"
              << "\tLepton from stop: "      << (m_l_from_stop == 0 ? "not found" : "found")
              << "\tLepton from anti-stop: " << (m_l_from_astp == 0 ? "not found" : "found")
              << "\n";
    return false;
  }

  // Find b quark from stop and anti-stop
  m_b_from_stop = 0;
  m_b_from_astp = 0;

  // I'm probably being overly careful here
  for (size_t quark_it = 0; quark_it != quark_list.size(); ++quark_it) {
    // if parent is stop
    if (quark_list.at(quark_it)->getParentPdgid() == +(1e6+6))  {
      if (m_b_from_stop == 0) m_b_from_stop = quark_list.at(quark_it);
      else
        std::cout << "WARNING! Found multiple quarks paired to stop!\n";
    }
    // if parent is anti-stop
    if (quark_list.at(quark_it)->getParentPdgid() == -(1e6+6)) {
      if (m_b_from_astp == 0) m_b_from_astp = quark_list.at(quark_it);
      else
        std::cout << "WARNING! Found multiple quarks paired to anti-stop!\n";
    }
  }
  if (m_b_from_stop == 0 || m_b_from_astp == 0) {
    std::cout << "ERROR!"
              << "\tB quark from stop: "      << (m_b_from_stop == 0 ? "not found" : "found")
              << "\tB quark from anti-stop: " << (m_b_from_astp == 0 ? "not found" : "found")
              << "\n";
    return false;
  }

  // check that the objects from the stop have the same parent barcode
  if (m_l_from_stop->getParentBarcode() != m_b_from_stop->getParentBarcode()) {
    std::cout << "ERROR! Lepton and b from stop have different parent barcodes\n";
    return false;
  }
  // check that the objects from the anti-stop have the same parent barcode
  if (m_l_from_astp->getParentBarcode() != m_b_from_astp->getParentBarcode()) {
    std::cout << "ERROR! Lepton and b from anti-stop have different parent barcodes\n";
    return false;
  }

  return true;
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

