#include "HistogramHandlers/include/HistogramHandlers.h"

#include "TruthNtupleLooper/include/Calculators.h"
#include "TruthNtupleLooper/include/ObjectDefs.h"
#include "TFile.h"

#include <iostream>
#include <algorithm>

// =============================================================================
// = Handle
// =============================================================================
HistogramHandlers::Handle::Handle()
{
  // do nothing
}

// -----------------------------------------------------------------------------
void HistogramHandlers::Handle::Fill( const TruthNtuple::FLAVOR_CHANNEL
                                    , const std::vector<TruthNtuple::Electron*>&
                                    , const std::vector<TruthNtuple::Muon*>&
                                    , const std::vector<TruthNtuple::Jet*>&
                                    , const TruthNtuple::Met&
                                    )
{
  // do nothing
}

// -----------------------------------------------------------------------------
void HistogramHandlers::Handle::write(TFile*)
{
  // do nothing
}

// =============================================================================
// = Flavor channel
// =============================================================================
HistogramHandlers::FlavorChannel::FlavorChannel() : HistogramHandlers::Handle()
{
  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    m_h_flavor_channel.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                            + "__flavor_channel"
                                            ).c_str()
                                          , ( "Flavor Channel - "
                                            + TruthNtuple::FlavorChannelStrings[fc_it]
                                            + "; Flavor Channel ; Entries"
                                            ).c_str()
                                          , TruthNtuple::FLAVOR_N, -0.5, TruthNtuple::FLAVOR_N - 0.5
                                          )
                                );
    for (int flavor_it = 0; flavor_it != TruthNtuple::FLAVOR_N; ++flavor_it) {
      m_h_flavor_channel.at(fc_it)->GetXaxis()->SetBinLabel( flavor_it+1
                                                           , TruthNtuple::FlavorChannelStrings[flavor_it].c_str()
                                                           );
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::FlavorChannel::Fill( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Electron*>&
         , const std::vector<TruthNtuple::Muon*>&
         , const std::vector<TruthNtuple::Jet*>&
         , const TruthNtuple::Met&
         )
{
  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {
      m_h_flavor_channel.at(fc)->Fill(flavor_channel);
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::FlavorChannel::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
      m_h_flavor_channel.at(fc_it)->Write();
  }
}

// =============================================================================
// = Object Multiplicity
// =============================================================================
HistogramHandlers::ObjectMultiplicity::ObjectMultiplicity() : HistogramHandlers::Handle()
{
  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    m_h_num_lep.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__num_lep"
                                     ).c_str()
                                   , ( "Lepton Multiplicity - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; Lepton multiplicity; Entries"
                                     ).c_str()
                                   , 5, -0.5, 4.5
                                   )
                         );

    m_h_num_jet.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__num_jet"
                                     ).c_str()
                                   , ( "Jet Multiplicity - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; Jet multiplicity; Entries"
                                     ).c_str()
                                   , 10, -0.5, 9.5
                                   )
                         );
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::ObjectMultiplicity::Fill( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Electron*>& el_list
         , const std::vector<TruthNtuple::Muon*>& mu_list
         , const std::vector<TruthNtuple::Jet*>& jet_list
         , const TruthNtuple::Met&
         )
{
  size_t num_lep = el_list.size()+mu_list.size();
  size_t num_jet = jet_list.size();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {
      m_h_num_lep.at(fc)->Fill(num_lep);
      m_h_num_jet.at(fc)->Fill(num_jet);
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::ObjectMultiplicity::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
      m_h_num_lep.at(fc_it)->Write();
      m_h_num_jet.at(fc_it)->Write();
  }
}

// =============================================================================
// = LeptonKinematics
// =============================================================================
HistogramHandlers::LeptonKinematics::LeptonKinematics() : HistogramHandlers::Handle()
{
  const int    pt_bins = 150;
  const double pt_min  = 0.;
  const double pt_max  = 1500.;

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
                                    + "__lep_pt_all"
                                    ).c_str()
                                  , ( "p_{T} - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; p_{T} [GeV] ; Entries"
                                    ).c_str()
                                  , pt_bins, pt_min, pt_max
                                  )
                        );
    m_h_pt_0.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__lep_pt_0"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{0} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_1.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__lep_pt_1"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{1} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__lep_pt_diff"
                                     ).c_str()
                                   , ( "p_{T} diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; p_{T}^{0} - p_{T}^{1} [GeV] ; Entries"
                                     ).c_str()
                                   , pt_bins, pt_min, pt_max
                                   )
                         );
    m_h_pt_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__lep_pt_2d"
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
                                     + "__lep_eta_all"
                                     ).c_str()
                                   , ( "#eta - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #eta ; Entries"
                                     ).c_str()
                                   , eta_bins, eta_min, eta_max
                                   )
                         );
    m_h_eta_0.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__lep_eta_0"
                                   ).c_str()
                                 , ( "#eta - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #eta^{0} ; Entries"
                                   ).c_str()
                                 , eta_bins, eta_min, eta_max
                                 )
                       );
    m_h_eta_1.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__lep_eta_1"
                                   ).c_str()
                                 , ( "#eta - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #eta^{1} ; Entries"
                                   ).c_str()
                                 , eta_bins, eta_min, eta_max
                                 )
                       );
    m_h_eta_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__lep_eta_diff"
                                      ).c_str()
                                    , ( "#eta diff - "
                                      + TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "; #eta^{0} - #eta^{1} ; Entries"
                                      ).c_str()
                                    , eta_bins/2, 0, eta_max
                                    )
                          );
    m_h_eta_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__lep_eta_2d"
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
                                     + "__lep_phi_all"
                                     ).c_str()
                                   , ( "#phi - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #phi ; Entries"
                                     ).c_str()
                                   , phi_bins, phi_min, phi_max
                                   )
                         );
    m_h_phi_0.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__lep_phi_0"
                                   ).c_str()
                                 , ( "#phi - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #phi^{0} ; Entries"
                                   ).c_str()
                                 , phi_bins, phi_min, phi_max
                                 )
                       );
    m_h_phi_1.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__lep_phi_1"
                                   ).c_str()
                                 , ( "#phi - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #phi^{1} ; Entries"
                                   ).c_str()
                                 , phi_bins, phi_min, phi_max
                                 )
                       );
    m_h_phi_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__lep_phi_diff"
                                      ).c_str()
                                    , ( "#phi diff - "
                                      + TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "; #phi^{0} - #phi^{1} ; Entries"
                                      ).c_str()
                                    , phi_bins/2, 0, phi_max
                                    )
                          );
    m_h_phi_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__lep_phi_2d"
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
void HistogramHandlers::LeptonKinematics::Fill( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Electron*>& el_list
         , const std::vector<TruthNtuple::Muon*>& mu_list
         , const std::vector<TruthNtuple::Jet*>& /*jet_list*/
         , const TruthNtuple::Met&
         )
{
  double pt_0 = 0;
  double pt_1 = 0;

  double eta_0 = 0;
  double eta_1 = 0;

  double phi_0 = 0;
  double phi_1 = 0;

  // get lepton pt values delending on the channel
  if (flavor_channel == TruthNtuple::FLAVOR_EE) {
    pt_0 = el_list.at(0)->getPt()/1.e3;
    pt_1 = el_list.at(1)->getPt()/1.e3;

    eta_0 = el_list.at(0)->getEta();
    eta_1 = el_list.at(1)->getEta();

    phi_0 = el_list.at(0)->getPhi();
    phi_1 = el_list.at(1)->getPhi();
  }
  else if (flavor_channel == TruthNtuple::FLAVOR_MM) {
    pt_0 = mu_list.at(0)->getPt()/1.e3;
    pt_1 = mu_list.at(1)->getPt()/1.e3;

    eta_0 = mu_list.at(0)->getEta();
    eta_1 = mu_list.at(1)->getEta();

    phi_0 = mu_list.at(0)->getPhi();
    phi_1 = mu_list.at(1)->getPhi();
  }
  else if (  flavor_channel == TruthNtuple::FLAVOR_EM
          || flavor_channel == TruthNtuple::FLAVOR_ME
          ) {
    pt_0 = el_list.at(0)->getPt()/1.e3;
    pt_1 = mu_list.at(0)->getPt()/1.e3;

    eta_0 = el_list.at(0)->getEta();
    eta_1 = mu_list.at(0)->getEta();

    phi_0 = el_list.at(0)->getPhi();
    phi_1 = mu_list.at(0)->getPhi();
  }
  else {
    return;
  }

  // check pt ordering
  if (pt_0 < pt_1) {
    std::swap(pt_0 , pt_1 );
    std::swap(eta_0, eta_1);
    std::swap(phi_0, phi_1);
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
void HistogramHandlers::LeptonKinematics::write(TFile* f)
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
// = JetKinematics
// =============================================================================
HistogramHandlers::JetKinematics::JetKinematics() : HistogramHandlers::Handle()
{
  const int    pt_bins = 150;
  const double pt_min  = 0.;
  const double pt_max  = 1500.;

  const int    eta_bins = 50;
  const double eta_min  = -5.;
  const double eta_max  = +5.;

  const int    phi_bins   = 64;
  const double phi_min = -3.2;
  const double phi_max = +3.2;


  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    // initialize pt histograms
    m_h_pt_all.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__jet_pt_all"
                                    ).c_str()
                                  , ( "p_{T} - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; p_{T} [GeV] ; Entries"
                                    ).c_str()
                                  , pt_bins, pt_min, pt_max
                                  )
                        );
    m_h_pt_0.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__jet_pt_0"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{0} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_1.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__jet_pt_1"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{1} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__jet_pt_diff"
                                     ).c_str()
                                   , ( "p_{T} diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; p_{T}^{0} - p_{T}^{1} [GeV] ; Entries"
                                     ).c_str()
                                   , pt_bins, pt_min, pt_max
                                   )
                         );
    m_h_pt_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__jet_pt_2d"
                                   ).c_str()
                                 , ( "p_{T} map - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; p_{T}^{0} [GeV] ; p_{T}^{1} [GeV]"
                                   ).c_str()
                                 , pt_bins, pt_min, pt_max
                                 , pt_bins, pt_min, pt_max
                                 )
                       );

    // initialize eta histograms
    m_h_eta_all.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__jet_eta_all"
                                     ).c_str()
                                   , ( "#eta - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #eta ; Entries"
                                     ).c_str()
                                   , eta_bins, eta_min, eta_max
                                   )
                         );
    m_h_eta_0.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__jet_eta_0"
                                   ).c_str()
                                 , ( "#eta - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #eta^{0} ; Entries"
                                   ).c_str()
                                 , eta_bins, eta_min, eta_max
                                 )
                       );
    m_h_eta_1.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__jet_eta_1"
                                   ).c_str()
                                 , ( "#eta - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #eta^{1} ; Entries"
                                   ).c_str()
                                 , eta_bins, eta_min, eta_max
                                 )
                       );
    m_h_eta_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__jet_eta_diff"
                                      ).c_str()
                                    , ( "#eta diff - "
                                      + TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "; #eta^{0} - #eta^{1} ; Entries"
                                      ).c_str()
                                    , eta_bins/2, 0, eta_max
                                    )
                          );
    m_h_eta_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__jet_eta_2d"
                                    ).c_str()
                                  , ( "#eta map - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; #eta^{0} ; #eta^{1}"
                                    ).c_str()
                                  , eta_bins, eta_min, eta_max
                                  , eta_bins, eta_min, eta_max
                                  )
                        );

    // initialize phi hisogramas
    m_h_phi_all.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__jet_phi_all"
                                     ).c_str()
                                   , ( "#phi - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #phi ; Entries"
                                     ).c_str()
                                   , phi_bins, phi_min, phi_max
                                   )
                         );
    m_h_phi_0.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__jet_phi_0"
                                   ).c_str()
                                 , ( "#phi - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #phi^{0} ; Entries"
                                   ).c_str()
                                 , phi_bins, phi_min, phi_max
                                 )
                       );
    m_h_phi_1.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__jet_phi_1"
                                   ).c_str()
                                 , ( "#phi - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #phi^{1} ; Entries"
                                   ).c_str()
                                 , phi_bins, phi_min, phi_max
                                 )
                       );
    m_h_phi_diff.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__jet_phi_diff"
                                      ).c_str()
                                    , ( "#phi diff - "
                                      + TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "; #phi^{0} - #phi^{1} ; Entries"
                                      ).c_str()
                                    , phi_bins/2, 0, phi_max
                                    )
                          );
    m_h_phi_2d.push_back( new TH2F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__jet_phi_2d"
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
void HistogramHandlers::JetKinematics::Fill( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Electron*>& /*el_list*/
         , const std::vector<TruthNtuple::Muon*>& /*mu_list*/
         , const std::vector<TruthNtuple::Jet*>& jet_list
         , const TruthNtuple::Met&
         )
{
  size_t num_jet = jet_list.size();

  double pt_0 = -1;
  double pt_1 = -1;

  double eta_0 = 0;
  double eta_1 = 0;

  double phi_0 = 0;
  double phi_1 = 0;

  if (num_jet > 0) {
    pt_0 = jet_list.at(0)->getPt()/1.e3;
    eta_0 = jet_list.at(0)->getEta();
    phi_0 = jet_list.at(0)->getPhi();
  }
  else {
    return;
  }
  if (num_jet > 1) {
    pt_1 = jet_list.at(1)->getPt()/1.e3;
    eta_1 = jet_list.at(1)->getEta();
    phi_1 = jet_list.at(1)->getPhi();
  }

  // check pt ordering
  if (pt_0 < pt_1) {
    std::swap(pt_0 , pt_1 );
    std::swap(eta_0, eta_1);
    std::swap(phi_0, phi_1);
  }

  // fill histograms
  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {
      // always fill entries for leading jet
      m_h_pt_all.at(fc)->Fill(pt_0);
      m_h_pt_0.at(fc)->Fill(pt_0);

      m_h_eta_all.at(fc)->Fill(eta_0);
      m_h_eta_0.at(fc)->Fill(eta_0);

      m_h_phi_all.at(fc)->Fill(phi_0);
      m_h_phi_0.at(fc)->Fill(phi_0);

      // only fill entries for subleading jet if there are more than one
      if (num_jet > 1) {
        m_h_pt_all.at(fc)->Fill(pt_1);
        m_h_pt_1.at(fc)->Fill(pt_1);
        m_h_pt_diff.at(fc)->Fill(pt_0 - pt_1);
        m_h_pt_2d.at(fc)->Fill(pt_0, pt_1);

        m_h_eta_all.at(fc)->Fill(eta_1);
        m_h_eta_1.at(fc)->Fill(eta_1);
        m_h_eta_diff.at(fc)->Fill(TruthNtuple::deltaEta(eta_0, eta_1));
        m_h_eta_2d.at(fc)->Fill(eta_0, eta_1);

        m_h_phi_all.at(fc)->Fill(phi_1);
        m_h_phi_1.at(fc)->Fill(phi_1);
        m_h_phi_diff.at(fc)->Fill(TruthNtuple::deltaPhi(phi_0, phi_1));
        m_h_phi_2d.at(fc)->Fill(phi_0, phi_1);
      }
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::JetKinematics::write(TFile* f)
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
// = Mll
// =============================================================================
HistogramHandlers::Mll::Mll() : HistogramHandlers::Handle()
{
  const int bins   = 150;
  const double min = 0;
  const double max = 1500;

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    m_h_mll.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                 + "__mll"
                                 ).c_str()
                               , ( "m_{ll} - "
                                 + TruthNtuple::FlavorChannelStrings[fc_it]
                                 + "; m_{ll} [GeV] ; Entries"
                                 ).c_str()
                               , bins, min, max
                               )
                     );
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::Mll::Fill( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Electron*>& el_list
         , const std::vector<TruthNtuple::Muon*>& mu_list
         , const std::vector<TruthNtuple::Jet*>& /*jet_list*/
         , const TruthNtuple::Met&
         )
{
  if (flavor_channel == TruthNtuple::FLAVOR_NONE) return;

  double mll = 0;

  if (flavor_channel == TruthNtuple::FLAVOR_EE) {
    mll = TruthNtuple::invariantMass(el_list.at(0), el_list.at(1))/1.e3;
  }
  else if (flavor_channel == TruthNtuple::FLAVOR_MM) {
    mll = TruthNtuple::invariantMass(mu_list.at(0), mu_list.at(1))/1.e3;
  }
  else if (  flavor_channel == TruthNtuple::FLAVOR_EM
          || flavor_channel == TruthNtuple::FLAVOR_ME
          ) {
    mll = TruthNtuple::invariantMass(el_list.at(0), mu_list.at(0))/1.e3;
  }

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {
      m_h_mll.at(fc)->Fill(mll);
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::Mll::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
      m_h_mll.at(fc_it)->Write();
  }
}

// =============================================================================
// = Met
// =============================================================================
HistogramHandlers::Met::Met() : HistogramHandlers::Handle()
{
  const int bins   = 50;
  const double min = 0;
  const double max = 500;

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    m_h_met.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                 + "__met"
                                 ).c_str()
                               , ( "E_{T}^{miss} - "
                                 + TruthNtuple::FlavorChannelStrings[fc_it]
                                 + "; E_{T}^{miss} [GeV] ; Entries"
                                 ).c_str()
                               , bins, min, max
                               )
                     );
    m_h_metrel.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__metrel"
                                    ).c_str()
                                  , ( "E_{T}^{miss,rel} - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; E_{T}^{miss,rel} [GeV] ; Entries"
                                    ).c_str()
                                  , bins, min, max
                                  )
                        );
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::Met::Fill( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Electron*>& /*el_list*/
         , const std::vector<TruthNtuple::Muon*>& /*mu_list*/
         , const std::vector<TruthNtuple::Jet*>& /*jet_list*/
         , const TruthNtuple::Met& met
         )
{
  if (flavor_channel == TruthNtuple::FLAVOR_NONE) return;

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {
      m_h_met.at(fc)->Fill(met.getMetNoint()/1.e3);
      m_h_metrel.at(fc)->Fill(met.getMetRelNoint()/1.e3);
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::Met::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
      m_h_met.at(fc_it)->Write();
      m_h_metrel.at(fc_it)->Write();
  }
}

// =============================================================================
// = Mjl
// =============================================================================
HistogramHandlers::Mjl::Mjl() : HistogramHandlers::Handle()
{
  const int bins   = 150;
  const double min = 0;
  const double max = 1500;

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    m_h_mjl_truth.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__mjl_truth"
                                       ).c_str()
                                     , ( "m_{jl}^{truth} - "
                                       + TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "; m_{jl}^{truth} [GeV] ; Entries"
                                       ).c_str()
                                     , bins, min, max
                                     )
                           );
    m_h_mjl_dphi_matching.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__mjl_dphi_matching"
                                               ).c_str()
                                             , ( "m_{jl} - "
                                               + TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "; m_{jl} [GeV] ; Entries"
                                               ).c_str()
                                             , bins, min, max
                                             )
                                   );
    m_h_mjl_dr_matching.push_back( new TH1F( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                             + "__mjl_dr_matching"
                                             ).c_str()
                                           , ( "m_{jl} - "
                                             + TruthNtuple::FlavorChannelStrings[fc_it]
                                             + "; m_{jl} [GeV] ; Entries"
                                             ).c_str()
                                           , bins, min, max
                                           )
                                 );
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::Mjl::Fill( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Electron*>& el_list
         , const std::vector<TruthNtuple::Muon*>& mu_list
         , const std::vector<TruthNtuple::Jet*>& jet_list
         , const TruthNtuple::Met&
         )
{
  if (flavor_channel == TruthNtuple::FLAVOR_NONE) return;

  // merge el and mu lists to lepton list
  std::vector<TruthNtuple::Lepton*> lep_list;
  lep_list.reserve(el_list.size() + mu_list.size());
  lep_list.insert(lep_list.end(), el_list.begin(), el_list.end());
  lep_list.insert(lep_list.end(), mu_list.begin(), mu_list.end());

  // get lists of mjl values for different matching methods
  std::vector<double> mjl_list_truth         = TruthNtuple::getInvariantMassList(lep_list, jet_list, 0);
  std::vector<double> mjl_list_dphi_matching = TruthNtuple::getInvariantMassList(lep_list, jet_list, 1);
  std::vector<double> mjl_list_dr_matching   = TruthNtuple::getInvariantMassList(lep_list, jet_list, 2);

  // loop over flavor channels
  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {
      // fill mjl truth histogram
      for (size_t mjl_it = 0; mjl_it != mjl_list_truth.size(); ++mjl_it) {
        m_h_mjl_truth.at(fc)->Fill(mjl_list_truth.at(mjl_it)/1.e3);
      }
      // fill mjl dphi matching histogram
      for (size_t mjl_it = 0; mjl_it != mjl_list_dphi_matching.size(); ++mjl_it) {
        m_h_mjl_dphi_matching.at(fc)->Fill(mjl_list_dphi_matching.at(mjl_it)/1.e3);
      }
      // fill mjl dr matching histogram
      for (size_t mjl_it = 0; mjl_it != mjl_list_dr_matching.size(); ++mjl_it) {
        m_h_mjl_dr_matching.at(fc)->Fill(mjl_list_dr_matching.at(mjl_it)/1.e3);
      }
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::Mjl::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
      m_h_mjl_truth.at(fc_it)->Write();
      m_h_mjl_dphi_matching.at(fc_it)->Write();
      m_h_mjl_dr_matching.at(fc_it)->Write();
  }
}
