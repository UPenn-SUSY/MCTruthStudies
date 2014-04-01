#include "BMinusLCutflow/include/BMinusLHistogramHandlers.h"
#include "HistogramHandlers/include/HistogramHandlers.h"

#include "TruthNtupleLooper/include/Calculators.h"
#include "TruthNtupleLooper/include/ObjectDefs.h"
#include "BMinusLCutflow/include/BMinusLCalculators.h"

#include "TFile.h"

#include <iostream>
#include <math.h>
#include <algorithm>

// =============================================================================
// = StopKinematics
// =============================================================================
HistogramHandlers::StopKinematics::StopKinematics() : HistogramHandlers::Handle()
{
  const int e_bins   = 300;
  const double e_min = 0.;
  const double e_max = 3000.;

  const int pt_bins   = 150;
  const double pt_min = 0.;
  const double pt_max = 1500.;

  const int pt_ratio_bins   = 100;
  const double pt_ratio_min = 0.;
  const double pt_ratio_max = 1.;

  const int eta_bins   = 50;
  const double eta_min = -5.;
  const double eta_max = +5.;

  const int y_bins   = 50;
  const double y_min = -5.;
  const double y_max = +5.;

  const int phi_bins   = 64;
  const double phi_min = -3.2;
  const double phi_max = +3.2;

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_e_com.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__stop_e_com"
                                    ).c_str()
                                  , ( "E - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; E [GeV] ; Entries"
                                    ).c_str()
                                  , e_bins, e_min, e_max
                                  )
                        );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_e_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__stop_e_all"
                                    ).c_str()
                                  , ( "E - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; E [GeV] ; Entries"
                                    ).c_str()
                                  , e_bins, e_min, e_max
                                  )
                        );
    m_h_e_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_e_stop"
                                  ).c_str()
                                , ( "E - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; E^{#tilde{t}} [GeV] ; Entries"
                                  ).c_str()
                                , e_bins, e_min, e_max
                                )
                      );
    m_h_e_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_e_astp"
                                  ).c_str()
                                , ( "e - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; E^{#tilde{t}*} [GeV] ; Entries"
                                  ).c_str()
                                , e_bins, e_min, e_max
                                )
                      );
    m_h_e_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_e_diff"
                                     ).c_str()
                                   , ( "E diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; E^{#tilde{t}} - E^{#tilde{t}*} [GeV] ; Entries"
                                     ).c_str()
                                   , 2*e_bins, -e_max, e_max
                                   )
                         );
    m_h_e_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__stop_e_2d"
                                   ).c_str()
                                 , ( "p maE - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; E^{#tilde{t}} [GeV] ; E^{#tilde{t}*} [GeV]"
                                   ).c_str()
                                 , e_bins, e_min, e_max
                                 , e_bins, e_min, e_max
                                 )
                       );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_p_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__stop_p_all"
                                    ).c_str()
                                  , ( "p - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; p [GeV] ; Entries"
                                    ).c_str()
                                  , pt_bins, pt_min, pt_max
                                  )
                        );
    m_h_p_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_p_stop"
                                  ).c_str()
                                , ( "p - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p^{#tilde{t}} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_p_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_p_astp"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p^{#tilde{t}*} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_p_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_p_diff"
                                     ).c_str()
                                   , ( "p diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; p^{#tilde{t}} - p^{#tilde{t}*} [GeV] ; Entries"
                                     ).c_str()
                                   , 2*pt_bins, -pt_max, pt_max
                                   )
                         );
    m_h_p_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__stop_p_2d"
                                   ).c_str()
                                 , ( "p map - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; p^{#tilde{t}} [GeV] ; p^{#tilde{t}*} [GeV]"
                                   ).c_str()
                                 , pt_bins, pt_min, pt_max
                                 , pt_bins, pt_min, pt_max
                                 )
                       );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_pt_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__stop_pt_all"
                                    ).c_str()
                                  , ( "p_{T} - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; p_{T} [GeV] ; Entries"
                                    ).c_str()
                                  , pt_bins, pt_min, pt_max
                                  )
                        );
    m_h_pt_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_pt_stop"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{#tilde{t}} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_pt_astp"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{#tilde{t}*} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_pt_diff"
                                     ).c_str()
                                   , ( "p_{T} diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; p_{T}^{#tilde{t}} - p_{T}^{#tilde{t}*} [GeV] ; Entries"
                                     ).c_str()
                                   , 2*pt_bins, -pt_max, pt_max
                                   )
                         );
    m_h_pt_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_p_over_m_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__stop_p_over_m_all"
                                    ).c_str()
                                  , ( "p - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; p/m ; Entries"
                                    ).c_str()
                                  , pt_ratio_bins, pt_ratio_min, pt_ratio_max
                                  )
                        );
    m_h_p_over_m_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_p_over_m_stop"
                                  ).c_str()
                                , ( "p - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p^{#tilde{t}}/m ; Entries"
                                  ).c_str()
                                , pt_ratio_bins, pt_ratio_min, pt_ratio_max
                                )
                      );
    m_h_p_over_m_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_p_over_m_astp"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p^{#tilde{t}*}/m ; Entries"
                                  ).c_str()
                                , pt_ratio_bins, pt_ratio_min, pt_ratio_max
                                )
                      );
    m_h_p_over_m_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_p_over_m_diff"
                                     ).c_str()
                                   , ( "p diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; p^{#tilde{t}}/m - p^{#tilde{t}*}/m ; Entries"
                                     ).c_str()
                                   , 2*pt_ratio_bins, -pt_ratio_max, pt_ratio_max
                                   )
                         );
    m_h_p_over_m_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__stop_p_over_m_2d"
                                   ).c_str()
                                 , ( "p map - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; p^{#tilde{t}}/m ; p^{#tilde{t}*}/m"
                                   ).c_str()
                                 , pt_ratio_bins, pt_ratio_min, pt_ratio_max
                                 , pt_ratio_bins, pt_ratio_min, pt_ratio_max
                                 )
                       );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_pt_over_m_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__stop_pt_over_m_all"
                                    ).c_str()
                                  , ( "p_{T} - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; p_{T}/m ; Entries"
                                    ).c_str()
                                  , pt_ratio_bins, pt_ratio_min, pt_ratio_max
                                  )
                        );
    m_h_pt_over_m_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_pt_over_m_stop"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{#tilde{t}}/m ; Entries"
                                  ).c_str()
                                , pt_ratio_bins, pt_ratio_min, pt_ratio_max
                                )
                      );
    m_h_pt_over_m_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_pt_over_m_astp"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{#tilde{t}*}/m ; Entries"
                                  ).c_str()
                                , pt_ratio_bins, pt_ratio_min, pt_ratio_max
                                )
                      );
    m_h_pt_over_m_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_pt_over_m_diff"
                                     ).c_str()
                                   , ( "p_{T} diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; p_{T}^{#tilde{t}}/m - p_{T}^{#tilde{t}*}/m ; Entries"
                                     ).c_str()
                                   , 2*pt_ratio_bins, -pt_ratio_max, pt_ratio_max
                                   )
                         );
    m_h_pt_over_m_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__stop_pt_over_m_2d"
                                   ).c_str()
                                 , ( "p_{T} map - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; p_{T}^{#tilde{t}}/, ; p_{T}^{#tilde{t}*}/m"
                                   ).c_str()
                                 , pt_ratio_bins, pt_ratio_min, pt_ratio_max
                                 , pt_ratio_bins, pt_ratio_min, pt_ratio_max
                                 )
                       );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_eta_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_eta_all"
                                     ).c_str()
                                   , ( "#eta - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #eta ; Entries"
                                     ).c_str()
                                   , eta_bins, eta_min, eta_max
                                   )
                         );
    m_h_eta_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_eta_stop"
                                  ).c_str()
                                , ( "#eta - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; #eta^{#tilde{t}} ; Entries"
                                  ).c_str()
                                , eta_bins, eta_min, eta_max
                                )
                      );
    m_h_eta_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_eta_astp"
                                  ).c_str()
                                , ( "#eta - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; #eta^{#tilde{t}*} ; Entries"
                                  ).c_str()
                                , eta_bins, eta_min, eta_max
                                )
                      );
    m_h_eta_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_eta_diff"
                                     ).c_str()
                                   , ( "#eta diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #eta^{#tilde{t}} - #eta^{#tilde{t}*} ; Entries"
                                     ).c_str()
                                   , eta_bins/2, 0, eta_max
                                   )
                         );
    m_h_eta_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_y_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_y_all"
                                     ).c_str()
                                   , ( "y - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; y ; Entries"
                                     ).c_str()
                                   , y_bins, y_min, y_max
                                   )
                         );
    m_h_y_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_y_stop"
                                  ).c_str()
                                , ( "y - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; y^{#tilde{t}} ; Entries"
                                  ).c_str()
                                , y_bins, y_min, y_max
                                )
                      );
    m_h_y_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__stop_y_astp"
                                  ).c_str()
                                , ( "y - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; y^{#tilde{t}*} ; Entries"
                                  ).c_str()
                                , y_bins, y_min, y_max
                                )
                      );
    m_h_y_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_y_diff"
                                     ).c_str()
                                   , ( "y diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; y^{#tilde{t}} - y^{#tilde{t}*} ; Entries"
                                     ).c_str()
                                   , y_bins/2, 0, y_max
                                   )
                         );
    m_h_y_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__stop_y_2d"
                                   ).c_str()
                                 , ( "y map - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; y^{#tilde{t}} ; y^{#tilde{t}*}"
                                   ).c_str()
                                 , y_bins, y_min, y_max
                                 , y_bins, y_min, y_max
                                 )
                       );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_phi_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__stop_phi_all"
                                     ).c_str()
                                   , ( "#phi - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #phi ; Entries"
                                     ).c_str()
                                   , phi_bins, phi_min, phi_max
                                   )
                         );
    m_h_phi_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__stop_phi_stop"
                                      ).c_str()
                                    , ( "#phi - "
                                      + TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "; #phi^{#tilde{t}} ; Entries"
                                      ).c_str()
                                    , phi_bins, phi_min, phi_max
                                    )
                          );
    m_h_phi_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__stop_phi_astp"
                                      ).c_str()
                                    , ( "#phi - "
                                      + TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "; #phi^{#tilde{t}*} ; Entries"
                                      ).c_str()
                                    , phi_bins, phi_min, phi_max
                                    )
                          );
    m_h_phi_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__stop_phi_diff"
                                      ).c_str()
                                    , ( "#phi diff - "
                                      + TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "; #phi^{#tilde{t}} - #phi^{#tilde{t}*} ; Entries"
                                      ).c_str()
                                    , phi_bins/2, 0, phi_max
                                    )
                          );
    m_h_phi_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

  double e_com = TruthNtuple::comEnergy(stop_list.at(0), stop_list.at(1))/1.e3;

  double e_stop = 0;
  double e_astp = 0;

  double m_stop = 0;
  double m_astp = 0;

  double p_stop = 0;
  double p_astp = 0;

  double pt_stop = 0;
  double pt_astp = 0;

  double eta_stop = 0;
  double eta_astp = 0;

  double y_stop = 0;
  double y_astp = 0;

  double phi_stop = 0;
  double phi_astp = 0;

  if (  stop_list.at(0)->getPdgid() == +(1e6+6)
     && stop_list.at(1)->getPdgid() == -(1e6+6)
     ) {
    e_stop = stop_list.at(0)->getE()/1.e3;
    e_astp = stop_list.at(1)->getE()/1.e3;

    m_stop = stop_list.at(0)->getM()/1.e3;
    m_astp = stop_list.at(1)->getM()/1.e3;

    p_stop = stop_list.at(0)->getP()/1.e3;
    p_astp = stop_list.at(1)->getP()/1.e3;

    pt_stop = stop_list.at(0)->getPt()/1.e3;
    pt_astp = stop_list.at(1)->getPt()/1.e3;

    eta_stop = stop_list.at(0)->getEta();
    eta_astp = stop_list.at(1)->getEta();

    y_stop = TruthNtuple::rapidity(stop_list.at(0));
    y_astp = TruthNtuple::rapidity(stop_list.at(1));

    phi_stop = stop_list.at(0)->getPhi();
    phi_astp = stop_list.at(1)->getPhi();
  }
  else if (  stop_list.at(0)->getPdgid() == -(1e6+6)
          && stop_list.at(1)->getPdgid() == +(1e6+6)
          ) {
    e_stop = stop_list.at(1)->getE()/1.e3;
    e_astp = stop_list.at(0)->getE()/1.e3;

    m_stop = stop_list.at(1)->getM()/1.e3;
    m_astp = stop_list.at(0)->getM()/1.e3;

    p_stop = stop_list.at(1)->getP()/1.e3;
    p_astp = stop_list.at(0)->getP()/1.e3;

    pt_stop = stop_list.at(1)->getPt()/1.e3;
    pt_astp = stop_list.at(0)->getPt()/1.e3;

    eta_stop = stop_list.at(1)->getEta();
    eta_astp = stop_list.at(0)->getEta();

    y_stop = TruthNtuple::rapidity(stop_list.at(1));
    y_astp = TruthNtuple::rapidity(stop_list.at(0));

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
      // fill E com
      m_h_e_com.at(fc)->Fill(e_com);

      // fill E
      m_h_e_all.at(fc)->Fill(e_stop);
      m_h_e_all.at(fc)->Fill(e_astp);

      m_h_e_stop.at(fc)->Fill(e_stop);
      m_h_e_astp.at(fc)->Fill(e_astp);

      m_h_e_diff.at(fc)->Fill(e_stop - e_astp);
      m_h_e_2d.at(fc)->Fill(e_stop, e_astp);

      // fill p
      m_h_p_all.at(fc)->Fill(p_stop);
      m_h_p_all.at(fc)->Fill(p_astp);

      m_h_p_stop.at(fc)->Fill(p_stop);
      m_h_p_astp.at(fc)->Fill(p_astp);

      m_h_p_diff.at(fc)->Fill(p_stop - p_astp);
      m_h_p_2d.at(fc)->Fill(p_stop, p_astp);

      // fill pt/m
      m_h_pt_over_m_all.at(fc)->Fill(pt_stop/m_stop);
      m_h_pt_over_m_all.at(fc)->Fill(pt_astp/m_astp);

      m_h_pt_over_m_stop.at(fc)->Fill(pt_stop/m_stop);
      m_h_pt_over_m_astp.at(fc)->Fill(pt_astp/m_astp);

      m_h_pt_over_m_diff.at(fc)->Fill(pt_stop/m_stop - pt_astp/m_astp);
      m_h_pt_over_m_2d.at(fc)->Fill(pt_stop/m_stop, pt_astp/m_astp);

      // fill p/m
      m_h_p_over_m_all.at(fc)->Fill(p_stop/m_stop);
      m_h_p_over_m_all.at(fc)->Fill(p_astp/m_astp);

      m_h_p_over_m_stop.at(fc)->Fill(p_stop/m_stop);
      m_h_p_over_m_astp.at(fc)->Fill(p_astp/m_astp);

      m_h_p_over_m_diff.at(fc)->Fill(p_stop/m_stop - p_astp/m_astp);
      m_h_p_over_m_2d.at(fc)->Fill(p_stop/m_stop, p_astp/m_astp);

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

      m_h_eta_diff.at(fc)->Fill(TruthNtuple::deltaEta(eta_stop, eta_astp));
      m_h_eta_2d.at(fc)->Fill(eta_stop, eta_astp);

      // Fill y
      m_h_y_all.at(fc)->Fill(y_stop);
      m_h_y_all.at(fc)->Fill(y_astp);

      m_h_y_stop.at(fc)->Fill(y_stop);
      m_h_y_astp.at(fc)->Fill(y_astp);

      m_h_y_diff.at(fc)->Fill(TruthNtuple::deltaEta(y_stop, y_astp));
      m_h_y_2d.at(fc)->Fill(y_stop, y_astp);

      // Fill phi
      m_h_phi_all.at(fc)->Fill(phi_stop);
      m_h_phi_all.at(fc)->Fill(phi_astp);

      m_h_phi_stop.at(fc)->Fill(phi_stop);
      m_h_phi_astp.at(fc)->Fill(phi_astp);

      m_h_phi_diff.at(fc)->Fill(TruthNtuple::deltaPhi(phi_stop, phi_astp));
      m_h_phi_2d.at(fc)->Fill(phi_stop, phi_astp);
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::StopKinematics::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    m_h_e_com.at(fc_it)->Write();
    m_h_e_all.at(fc_it)->Write();
    m_h_e_stop.at(fc_it)->Write();
    m_h_e_astp.at(fc_it)->Write();
    m_h_e_diff.at(fc_it)->Write();
    m_h_e_2d.at(fc_it)->Write();

    m_h_p_all.at(fc_it)->Write();
    m_h_p_stop.at(fc_it)->Write();
    m_h_p_astp.at(fc_it)->Write();
    m_h_p_diff.at(fc_it)->Write();
    m_h_p_2d.at(fc_it)->Write();

    m_h_pt_all.at(fc_it)->Write();
    m_h_pt_stop.at(fc_it)->Write();
    m_h_pt_astp.at(fc_it)->Write();
    m_h_pt_diff.at(fc_it)->Write();
    m_h_pt_2d.at(fc_it)->Write();

    m_h_p_over_m_all.at(fc_it)->Write();
    m_h_p_over_m_stop.at(fc_it)->Write();
    m_h_p_over_m_astp.at(fc_it)->Write();
    m_h_p_over_m_diff.at(fc_it)->Write();
    m_h_p_over_m_2d.at(fc_it)->Write();

    m_h_pt_over_m_all.at(fc_it)->Write();
    m_h_pt_over_m_stop.at(fc_it)->Write();
    m_h_pt_over_m_astp.at(fc_it)->Write();
    m_h_pt_over_m_diff.at(fc_it)->Write();
    m_h_pt_over_m_2d.at(fc_it)->Write();

    m_h_eta_all.at(fc_it)->Write();
    m_h_eta_stop.at(fc_it)->Write();
    m_h_eta_astp.at(fc_it)->Write();
    m_h_eta_diff.at(fc_it)->Write();
    m_h_eta_2d.at(fc_it)->Write();

    m_h_y_all.at(fc_it)->Write();
    m_h_y_stop.at(fc_it)->Write();
    m_h_y_astp.at(fc_it)->Write();
    m_h_y_diff.at(fc_it)->Write();
    m_h_y_2d.at(fc_it)->Write();

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
    m_h_pt_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "__quark_pt_all"
                                    ).c_str()
                                  , ( "p_{T} - "
                                    + TruthNtuple::FlavorChannelStrings[fc_it]
                                    + "; p_{T} [GeV] ; Entries"
                                    ).c_str()
                                  , pt_bins, pt_min, pt_max
                                  )
                        );
    m_h_pt_0.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__quark_pt_0"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{0} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_1.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "__quark_pt_1"
                                  ).c_str()
                                , ( "p_{T} - "
                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                  + "; p_{T}^{1} [GeV] ; Entries"
                                  ).c_str()
                                , pt_bins, pt_min, pt_max
                                )
                      );
    m_h_pt_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__quark_pt_diff"
                                     ).c_str()
                                   , ( "p_{T} diff - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; p_{T}^{0} - p_{T}^{1} [GeV] ; Entries"
                                     ).c_str()
                                   , pt_bins, pt_min, pt_max
                                   )
                         );
    m_h_pt_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_eta_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__quark_eta_all"
                                     ).c_str()
                                   , ( "#eta - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #eta ; Entries"
                                     ).c_str()
                                   , eta_bins, eta_min, eta_max
                                   )
                         );
    m_h_eta_0.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__quark_eta_0"
                                   ).c_str()
                                 , ( "#eta - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #eta^{0} ; Entries"
                                   ).c_str()
                                 , eta_bins, eta_min, eta_max
                                 )
                       );
    m_h_eta_1.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__quark_eta_1"
                                   ).c_str()
                                 , ( "#eta - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #eta^{1} ; Entries"
                                   ).c_str()
                                 , eta_bins, eta_min, eta_max
                                 )
                       );
    m_h_eta_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__quark_eta_diff"
                                      ).c_str()
                                    , ( "#eta diff - "
                                      + TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "; #eta^{0} - #eta^{1} ; Entries"
                                      ).c_str()
                                    , eta_bins/2, 0, eta_max
                                    )
                          );
    m_h_eta_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_phi_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__quark_phi_all"
                                     ).c_str()
                                   , ( "#phi - "
                                     + TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "; #phi ; Entries"
                                     ).c_str()
                                   , phi_bins, phi_min, phi_max
                                   )
                         );
    m_h_phi_0.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__quark_phi_0"
                                   ).c_str()
                                 , ( "#phi - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #phi^{0} ; Entries"
                                   ).c_str()
                                 , phi_bins, phi_min, phi_max
                                 )
                       );
    m_h_phi_1.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "__quark_phi_1"
                                   ).c_str()
                                 , ( "#phi - "
                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                   + "; #phi^{1} ; Entries"
                                   ).c_str()
                                 , phi_bins, phi_min, phi_max
                                 )
                       );
    m_h_phi_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__quark_phi_diff"
                                      ).c_str()
                                    , ( "#phi diff - "
                                      + TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "; #phi^{0} - #phi^{1} ; Entries"
                                      ).c_str()
                                    , phi_bins/2, 0, phi_max
                                    )
                          );
    m_h_phi_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
  const int    pt_bins = 150;
  const double pt_min  = 0.;
  const double pt_max  = 1500.;

  const int    eta_bins = 50;
  const double eta_min  = -5.;
  const double eta_max  = +5.;

  const int    phi_bins = 64;
  const double phi_min  = -3.2;
  const double phi_max  = +3.2;

  const int    dr_bins = 50;
  const double dr_min  = 0.;
  const double dr_max  = +5.;

  const int    mbl_bins = 150;
  const double mbl_min  = 0.;
  const double mbl_max  = 1500.;

  const int    mbl_ratio_bins = 110;
  const double mbl_ratio_min  = 0.;
  const double mbl_ratio_max  = 1.1;

  const int    mbl_sq_sum_bins = 300;
  const double mbl_sq_sum_min  = 0.;
  const double mbl_sq_sum_max  = 3000.;

  const int    ptbl_bins = 150;
  const double ptbl_min  = 0.;
  const double ptbl_max  = 1500.;

  const int    ptbl_sq_sum_bins = 300;
  const double ptbl_sq_sum_min  = 0.;
  const double ptbl_sq_sum_max  = 3000.;

  const int    ptbl_ratio_bins = 110;
  const double ptbl_ratio_min  = 0.;
  const double ptbl_ratio_max  = 1.1;

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_l_pt_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__l_pt_all" // name suffix
                                      ).c_str()
                                    , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + " - p_{T}" // title suffix
                                      + " ; p_{T} [GeV]" // x-axis label
                                      + " ; Entries" // y-axis label
                                      ).c_str()
                                    , pt_bins, pt_min, pt_max
                                    )
                          );
    m_h_l_pt_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_l_pt_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_l_pt_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__l_pt_diff" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - p_{T}" // title suffix
                                       + " ; |p_{T}^{from #tilde{t}} - p_{T}^{from #tilde{t}*}| [GeV]" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , pt_bins, pt_min, pt_max
                                     )
                           );
    m_h_l_pt_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__l_pt_2d" // name suffix
                                     ).c_str()
                                   , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + " - p_{T}" // title suffix
                                     + " ; p_{T}^{from #tilde{t}} [GeV]" // x-axis label
                                     + " ; p_{T}^{from #tilde{t}*} [GeV]" // y-axis label
                                     ).c_str()
                                   , pt_bins, pt_min, pt_max
                                   , pt_bins, pt_min, pt_max
                                   )
                         );

    m_h_l_eta_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__l_eta_all" // name suffix
                                      ).c_str()
                                    , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + " - #eta" // title suffix
                                      + " ; #eta" // x-axis label
                                      + " ; Entries" // y-axis label
                                      ).c_str()
                                    , eta_bins, eta_min, eta_max
                                    )
                          );
    m_h_l_eta_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_l_eta_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_l_eta_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__l_eta_diff" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #eta" // title suffix
                                       + " ; |#eta^{from #tilde{t}} - #eta^{from #tilde{t}*}|" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , eta_bins/2, 0, eta_max
                                     )
                           );
    m_h_l_eta_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__l_eta_2d" // name suffix
                                     ).c_str()
                                   , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + " - #eta" // title suffix
                                     + " ; #eta^{from #tilde{t}}" // x-axis label
                                     + " ; #eta^{from #tilde{t}*}" // y-axis label
                                     ).c_str()
                                   , eta_bins, eta_min, eta_max
                                   , eta_bins, eta_min, eta_max
                                   )
                         );

    m_h_l_phi_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__l_phi_all" // name suffix
                                      ).c_str()
                                    , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + " - #phi" // title suffix
                                      + " ; #phi" // x-axis label
                                      + " ; Entries" // y-axis label
                                      ).c_str()
                                    , phi_bins, phi_min, phi_max
                                    )
                          );
    m_h_l_phi_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_l_phi_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_l_phi_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__l_phi_diff" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #phi" // title suffix
                                       + " ; |#phi^{from #tilde{t}} - #phi^{from #tilde{t}*}|" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , phi_bins/2, 0, phi_max
                                     )
                           );
    m_h_l_phi_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__l_phi_2d" // name suffix
                                     ).c_str()
                                   , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + " - #phi" // title suffix
                                     + " ; #phi^{from #tilde{t}}" // x-axis label
                                     + " ; #phi^{from #tilde{t}*}" // y-axis label
                                     ).c_str()
                                   , phi_bins, phi_min, phi_max
                                   , phi_bins, phi_min, phi_max
                                   )
                         );


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_b_pt_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__b_pt_all" // name suffix
                                      ).c_str()
                                    , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + " - p_{T}" // title suffix
                                      + " ; p_{T} [GeV]" // x-axis label
                                      + " ; Entries" // y-axis label
                                      ).c_str()
                                    , pt_bins, pt_min, pt_max
                                    )
                          );
    m_h_b_pt_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_b_pt_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_b_pt_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__b_pt_diff" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - p_{T}" // title suffix
                                       + " ; |p_{T}^{from #tilde{t}} - p_{T}^{from #tilde{t}*}| [GeV]" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , pt_bins, pt_min, pt_max
                                     )
                           );
    m_h_b_pt_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__b_pt_2d" // name suffix
                                     ).c_str()
                                   , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + " - p_{T}" // title suffix
                                     + " ; p_{T}^{from #tilde{t}} [GeV]" // x-axis label
                                     + " ; p_{T}^{from #tilde{t}*} [GeV]" // y-axis label
                                     ).c_str()
                                   , pt_bins, pt_min, pt_max
                                   , pt_bins, pt_min, pt_max
                                   )
                         );

    m_h_b_eta_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + "__b_eta_all" // name suffix
                                      ).c_str()
                                    , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                      + " - #eta" // title suffix
                                      + " ; #eta" // x-axis label
                                      + " ; Entries" // y-axis label
                                      ).c_str()
                                    , eta_bins, eta_min, eta_max
                                    )
                          );
    m_h_b_eta_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_b_eta_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_b_eta_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__b_eta_diff" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #eta" // title suffix
                                       + " ; |#eta^{#tilde{t}} - #eta^{#tilde{t}*}|" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , eta_bins/2, 0, eta_max
                                     )
                           );
    m_h_b_eta_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__b_eta_2d" // name suffix
                                     ).c_str()
                                   , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + " - #eta" // title suffix
                                     + " ; #eta^{#tilde{t}}" // x-axis label
                                     + " ; #eta^{#tilde{t}*}" // y-axis label
                                     ).c_str()
                                   , eta_bins, eta_min, eta_max
                                   , eta_bins, eta_min, eta_max
                                   )
                         );

    m_h_b_phi_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__b_phi_all" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #phi" // title suffix
                                       + " ; #phi" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , phi_bins, phi_min, phi_max
                                     )
                           );
    m_h_b_phi_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_b_phi_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_b_phi_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + "__b_phi_diff" // name suffix
                                       ).c_str()
                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                       + " - #phi" // title suffix
                                       + " ; |#phi^{from #tilde{t}} - #phi^{from #tilde{t}*}|" // x-axis label
                                       + " ; Entries" // y-axis label
                                       ).c_str()
                                     , phi_bins/2, 0, phi_max
                                     )
                           );
    m_h_b_phi_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + "__b_phi_2d" // name suffix
                                     ).c_str()
                                   , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                     + " - #phi" // title suffix
                                     + " ; #phi^{from #tilde{t}}" // x-axis label
                                     + " ; #phi^{from #tilde{t}*}" // y-axis label
                                     ).c_str()
                                   , phi_bins, phi_min, phi_max
                                   , phi_bins, phi_min, phi_max
                                   )
                         );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_right_pair_bl_dpt.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__right_pair_bl_dpt" // name suffix
                                               ).c_str()
                                             , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + " - p_{T}" // title suffix
                                               + " ; p_{T}^{l} - p_{T}^{b} [GeV]" // x-axis label
                                               + " ; Entries" // y-axis label
                                               ).c_str()
                                             , 2*pt_bins, -pt_max, pt_max
                                             )
                                   );
    m_h_right_pair_bl_deta.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__right_pair_bl_deta" // name suffix
                                                ).c_str()
                                              , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + " - #eta" // title suffix
                                                + " ; |#eta^{l} - #eta^{b}|" // x-axis label
                                                + " ; Entries" // y-axis label
                                                ).c_str()
                                              , eta_bins/2, 0, eta_max
                                              )
                                    );
    m_h_right_pair_bl_dphi.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__right_pair_bl_dphi" // name suffix
                                                ).c_str()
                                              , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + " - #phi" // title suffix
                                                + " ; |#phi^{l} - #phi^{b}|" // x-axis label
                                                + " ; Entries" // y-axis label
                                                ).c_str()
                                              , phi_bins/2, 0, phi_max
                                              )
                                    );
    m_h_right_pair_bl_dr.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + "__right_pair_bl_dr" // name suffix
                                              ).c_str()
                                            , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + " - #delta(R)" // title suffix
                                              + " ; #delta R(l,b)" // x-axis label
                                              + " ; Entries" // y-axis label
                                              ).c_str()
                                            , dr_bins, dr_min, dr_max
                                            )
                                  );

    m_h_right_pair_bl_pt_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_right_pair_bl_eta_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_right_pair_bl_phi_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_wrong_pair_bl_dpt.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__wrong_pair_bl_dpt" // name suffix
                                               ).c_str()
                                             , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + " - p_{T}" // title suffix
                                               + " ; p_{T}^{l} - p_{T}^{b} [GeV]" // x-axis label
                                               + " ; Entries" // y-axis label
                                               ).c_str()
                                             , 2*pt_bins, -pt_max, pt_max
                                             )
                                   );
    m_h_wrong_pair_bl_deta.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__wrong_pair_bl_deta" // name suffix
                                                ).c_str()
                                              , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + " - #eta" // title suffix
                                                + " ; |#eta^{l} - #eta^{b}|" // x-axis label
                                                + " ; Entries" // y-axis label
                                                ).c_str()
                                              , eta_bins/2, 0, eta_max
                                              )
                                    );
    m_h_wrong_pair_bl_dphi.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__wrong_pair_bl_dphi" // name suffix
                                                ).c_str()
                                              , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + " - #phi" // title suffix
                                                + " ; |#phi^{l} - #phi^{b}|" // x-axis label
                                                + " ; Entries" // y-axis label
                                                ).c_str()
                                              , phi_bins/2, 0, phi_max
                                              )
                                    );
    m_h_wrong_pair_bl_dr.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + "__wrong_pair_bl_dr" // name suffix
                                              ).c_str()
                                            , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + " - #delta(R)" // title suffix
                                              + " ; #delta R(l,b)" // x-axis label
                                              + " ; Entries" // y-axis label
                                              ).c_str()
                                            , dr_bins, dr_min, dr_max
                                            )
                                  );

    m_h_wrong_pair_bl_pt_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    m_h_wrong_pair_bl_eta_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    m_h_wrong_pair_bl_phi_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_right_pair_mbl_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    m_h_right_pair_mbl_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    m_h_right_pair_mbl_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    m_h_right_pair_mbl_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    m_h_right_pair_mbl_ratio.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__right_pair_mbl_ratio" // name suffix
                                                  ).c_str()
                                                , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - mbl" // title suffix
                                                  + "; m_{bl}^{1}/m_{bl}^{0}" // x-axis label
                                                  + "; Entries" // y-axis label
                                                  ).c_str()
                                                , mbl_ratio_bins, mbl_ratio_min, mbl_ratio_max
                                                )
                                      );

    m_h_right_pair_mbl_sq_sum.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + "__right_pair_mbl_sq_sum" // name suffix
                                                   ).c_str()
                                                 , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + " - mbl" // title suffix
                                                   + "; #sqrt{m_{bl}^{#tilde{t} 2} + m_{bl}^{#tilde{t}* 2}} [GeV]" // x-axis label
                                                   + "; Entries" // y-axis label
                                                   ).c_str()
                                                 , mbl_sq_sum_bins, mbl_sq_sum_min, mbl_sq_sum_max
                                                 )
                                       );

    m_h_right_pair_mbl_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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
    m_h_wrong_pair_mbl_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    m_h_wrong_pair_mbl_0.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    m_h_wrong_pair_mbl_1.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    m_h_wrong_pair_mbl_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    m_h_wrong_pair_mbl_ratio.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__wrong_pair_mbl_ratio" // name suffix
                                                  ).c_str()
                                                , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - mbl" // title suffix
                                                  + "; m_{bl}^{1}/m_{bl}^{0}" // x-axis label
                                                  + "; Entries" // y-axis label
                                                  ).c_str()
                                                , mbl_ratio_bins, mbl_ratio_min, mbl_ratio_max
                                                )
                                      );
    m_h_wrong_pair_mbl_sq_sum.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + "__wrong_pair_mbl_sq_sum" // name suffix
                                                   ).c_str()
                                                 , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + " - mbl" // title suffix
                                                   + "; #sqrt{m_{bl}^{0 2} + m_{bl}^{1 2}} [GeV]" // x-axis label
                                                   + "; Entries" // y-axis label
                                                   ).c_str()
                                                 , mbl_sq_sum_bins, mbl_sq_sum_min, mbl_sq_sum_max
                                                 )
                                       );


    m_h_wrong_pair_mbl_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_diff_pair_mbl_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__diff_pair_mbl_all" // name suffix
                                               ).c_str()
                                             , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + " - mbl" // title suffix
                                               + "; m_{bl} [GeV]" // x-axis label
                                               + "; Entries" // y-axis label
                                               ).c_str()
                                             , mbl_bins, mbl_min, mbl_max
                                             )
                                   );

    m_h_diff_pair_mbl_0.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                             + "__diff_pair_mbl_0" // name suffix
                                             ).c_str()
                                           , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                             + " - mbl" // title suffix
                                             + "; m_{bl}^{0} [GeV]" // x-axis label
                                             + "; Entries" // y-axis label
                                             ).c_str()
                                           , mbl_bins, mbl_min, mbl_max
                                           )
                                 );

    m_h_diff_pair_mbl_1.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                             + "__diff_pair_mbl_1" // name suffix
                                             ).c_str()
                                           , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                             + " - mbl" // title suffix
                                             + "; m_{bl}^{1} [GeV]" // x-axis label
                                             + "; Entries" // y-axis label
                                             ).c_str()
                                           , mbl_bins, mbl_min, mbl_max
                                           )
                                 );

    m_h_diff_pair_mbl_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__diff_pair_mbl_diff" // name suffix
                                                ).c_str()
                                              , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + " - mbl" // title suffix
                                                + "; |m_{bl}^{0} - m_{bl}^{1}| [GeV]" // x-axis label
                                                + "; Entries" // y-axis label
                                                ).c_str()
                                              , mbl_bins, mbl_min, mbl_max
                                              )
                                    );

    m_h_diff_pair_mbl_ratio.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__diff_pair_mbl_ratio" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - mbl" // title suffix
                                                 + "; m_{bl}^{1}/m_{bl}^{0}" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , mbl_ratio_bins, mbl_ratio_min, mbl_ratio_max
                                               )
                                     );
    m_h_diff_pair_mbl_sq_sum.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__diff_pair_mbl_sq_sum" // name suffix
                                                  ).c_str()
                                                , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - mbl" // title suffix
                                                  + "; #sqrt{m_{bl}^{0 2} + m_{bl}^{1 2}} [GeV]" // x-axis label
                                                  + "; Entries" // y-axis label
                                                  ).c_str()
                                                , mbl_sq_sum_bins, mbl_sq_sum_min, mbl_sq_sum_max
                                                )
                                      );

    m_h_diff_pair_mbl_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + "__diff_pair_mbl_2d" // name suffix
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

    m_h_diff_pair_cor_pairing.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + "__diff_pair_cor_pairing" // name suffix
                                                   ).c_str()
                                                 , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + " - pairing eff" // title suffix
                                                   + "; pairing" // x-axis label
                                                   + "; Entries" // y-axis label
                                                   ).c_str()
                                                 , 2, -0.5, 1.5
                                                 )
                                       );
    m_h_diff_pair_cor_pairing.at(m_h_diff_pair_cor_pairing.size()-1)->GetXaxis()->SetBinLabel(1, "Incorrect");
    m_h_diff_pair_cor_pairing.at(m_h_diff_pair_cor_pairing.size()-1)->GetXaxis()->SetBinLabel(2, "Correct");

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_ratio_pair_mbl_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__ratio_pair_mbl_all" // name suffix
                                                ).c_str()
                                              , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + " - mbl" // title suffix
                                                + "; m_{bl} [GeV]" // x-axis label
                                                + "; Entries" // y-axis label
                                                ).c_str()
                                              , mbl_bins, mbl_min, mbl_max
                                              )
                                    );

    m_h_ratio_pair_mbl_0.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + "__ratio_pair_mbl_0" // name suffix
                                              ).c_str()
                                            , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + " - mbl" // title suffix
                                              + "; m_{bl}^{0} [GeV]" // x-axis label
                                              + "; Entries" // y-axis label
                                              ).c_str()
                                            , mbl_bins, mbl_min, mbl_max
                                            )
                                  );

    m_h_ratio_pair_mbl_1.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + "__ratio_pair_mbl_1" // name suffix
                                              ).c_str()
                                            , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + " - mbl" // title suffix
                                              + "; m_{bl}^{1} [GeV]" // x-axis label
                                              + "; Entries" // y-axis label
                                              ).c_str()
                                            , mbl_bins, mbl_min, mbl_max
                                            )
                                  );

    m_h_ratio_pair_mbl_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__ratio_pair_mbl_diff" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - mbl" // title suffix
                                                 + "; |m_{bl}^{0} - m_{bl}^{1}| [GeV]" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , mbl_bins, mbl_min, mbl_max
                                               )
                                     );

    m_h_ratio_pair_mbl_ratio.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__ratio_pair_mbl_ratio" // name suffix
                                                  ).c_str()
                                                , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - mbl" // title suffix
                                                  + "; m_{bl}^{1}/m_{bl}^{0}" // x-axis label
                                                  + "; Entries" // y-axis label
                                                  ).c_str()
                                                , mbl_ratio_bins, mbl_ratio_min, mbl_ratio_max
                                                )
                                      );
    m_h_ratio_pair_mbl_sq_sum.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + "__ratio_pair_mbl_sq_sum" // name suffix
                                                   ).c_str()
                                                 , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + " - mbl" // title suffix
                                                   + "; #sqrt{m_{bl}^{0 2} + m_{bl}^{1 2}} [GeV]" // x-axis label
                                                   + "; Entries" // y-axis label
                                                   ).c_str()
                                                 , mbl_sq_sum_bins, mbl_sq_sum_min, mbl_sq_sum_max
                                                 )
                                       );


    m_h_ratio_pair_mbl_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__ratio_pair_mbl_2d" // name suffix
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

    m_h_ratio_pair_cor_pairing.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                    + "__ratio_pair_cor_pairing" // name suffix
                                                    ).c_str()
                                                  , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                    + " - pairing eff" // title suffix
                                                    + "; pairing" // x-axis label
                                                    + "; Entries" // y-axis label
                                                    ).c_str()
                                                  , 2, -0.5, 1.5
                                                  )
                                        );
    m_h_ratio_pair_cor_pairing.at(m_h_ratio_pair_cor_pairing.size()-1)->GetXaxis()->SetBinLabel(1, "Incorrect");
    m_h_ratio_pair_cor_pairing.at(m_h_ratio_pair_cor_pairing.size()-1)->GetXaxis()->SetBinLabel(2, "Correct");

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_sq_sum_pair_mbl_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__sq_sum_pair_mbl_all" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - mbl" // title suffix
                                                 + "; m_{bl} [GeV]" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , mbl_bins, mbl_min, mbl_max
                                               )
                                     );

    m_h_sq_sum_pair_mbl_0.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__sq_sum_pair_mbl_0" // name suffix
                                               ).c_str()
                                             , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + " - mbl" // title suffix
                                               + "; m_{bl}^{0} [GeV]" // x-axis label
                                               + "; Entries" // y-axis label
                                               ).c_str()
                                             , mbl_bins, mbl_min, mbl_max
                                             )
                                   );

    m_h_sq_sum_pair_mbl_1.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__sq_sum_pair_mbl_1" // name suffix
                                               ).c_str()
                                             , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + " - mbl" // title suffix
                                               + "; m_{bl}^{1} [GeV]" // x-axis label
                                               + "; Entries" // y-axis label
                                               ).c_str()
                                             , mbl_bins, mbl_min, mbl_max
                                             )
                                   );

    m_h_sq_sum_pair_mbl_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__sq_sum_pair_mbl_diff" // name suffix
                                                  ).c_str()
                                                , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - mbl" // title suffix
                                                  + "; |m_{bl}^{0} - m_{bl}^{1}| [GeV]" // x-axis label
                                                  + "; Entries" // y-axis label
                                                  ).c_str()
                                                , mbl_bins, mbl_min, mbl_max
                                                )
                                      );

    m_h_sq_sum_pair_mbl_ratio.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + "__sq_sum_pair_mbl_ratio" // name suffix
                                                   ).c_str()
                                                 , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + " - mbl" // title suffix
                                                   + "; m_{bl}^{1}/m_{bl}^{0}" // x-axis label
                                                   + "; Entries" // y-axis label
                                                   ).c_str()
                                                 , mbl_ratio_bins, mbl_ratio_min, mbl_ratio_max
                                                 )
                                       );
    m_h_sq_sum_pair_mbl_sq_sum.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                    + "__sq_sum_pair_mbl_sq_sum" // name suffix
                                                    ).c_str()
                                                  , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                    + " - mbl" // title suffix
                                                    + "; #sqrt{m_{bl}^{0 2} + m_{bl}^{1 2}} [GeV]" // x-axis label
                                                    + "; Entries" // y-axis label
                                                    ).c_str()
                                                  , mbl_sq_sum_bins, mbl_sq_sum_min, mbl_sq_sum_max
                                                  )
                                        );


    m_h_sq_sum_pair_mbl_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__sq_sum_pair_mbl_2d" // name suffix
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

    m_h_sq_sum_pair_cor_pairing.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                     + "__sq_sum_pair_cor_pairing" // name suffix
                                                     ).c_str()
                                                   , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                     + " - pairing eff" // title suffix
                                                     + "; pairing" // x-axis label
                                                     + "; Entries" // y-axis label
                                                     ).c_str()
                                                   , 2, -0.5, 1.5
                                                   )
                                         );
    m_h_sq_sum_pair_cor_pairing.at(m_h_ratio_pair_cor_pairing.size()-1)->GetXaxis()->SetBinLabel(1, "Incorrect");
    m_h_sq_sum_pair_cor_pairing.at(m_h_ratio_pair_cor_pairing.size()-1)->GetXaxis()->SetBinLabel(2, "Correct");

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_vs_right_mbl_diff.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                     + "__wrong_vs_right_mbl_diff" // name suffix
                                                     ).c_str()
                                                   , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                     + " - mbl" // title suffix
                                                     + "; |m_{bl}^{#tilde{t}} - m_{bl}^{#tilde{t}*}|(right pairing) [GeV]" // x-axis label
                                                     + "; |m_{bl}^{0} - m_{bl}^{1}|(wrong pairing) [GeV]" // y-axis label
                                                     ).c_str()
                                                   , mbl_bins, mbl_min, mbl_max
                                                   , mbl_bins, mbl_min, mbl_max
                                                   )
                                         );

    m_h_wrong_vs_right_mbl_ratio.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                      + "__wrong_vs_right_mbl_ratio" // name suffix
                                                      ).c_str()
                                                    , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                      + " - mbl" // title suffix
                                                      + "; m_{bl}^{1}/m_{bl}^{0} (right pairing)" // x-axis label
                                                      + "; m_{bl}^{1}/m_{bl}^{0} (wrong pairing)" // y-axis label
                                                      ).c_str()
                                                    , mbl_ratio_bins, mbl_ratio_min, mbl_ratio_max
                                                    , mbl_ratio_bins, mbl_ratio_min, mbl_ratio_max
                                                    )
                                          );

    m_h_wrong_vs_right_mbl_sq_sum.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                       + "__wrong_vs_right_mbl_sq_sum" // name suffix
                                                       ).c_str()
                                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                       + " - mbl" // title suffix
                                                       + "; #sqrt{m_{bl}^{#tilde{t} 2} + m_{bl}^{#tilde{t}* 2}}(right pairing) [GeV]" // x-axis label
                                                       + "; #sqrt{m_{bl}^{0 2} + m_{bl}^{1 2}}(wrong pairing) [GeV]" // y-axis label
                                                       ).c_str()
                                                     , mbl_sq_sum_bins, mbl_sq_sum_min, mbl_sq_sum_max
                                                     , mbl_sq_sum_bins, mbl_sq_sum_min, mbl_sq_sum_max
                                                     )
                                           );


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_right_pair_ptbl_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__right_pair_ptbl_all" // name suffix
                                                ).c_str()
                                              , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + " - ptbl" // title suffix
                                                + "; p_{T}^{bl} [GeV]" // x-axis label
                                                + "; Entries" // y-axis label
                                                ).c_str()
                                              , ptbl_bins, ptbl_min, ptbl_max
                                              )
                                    );

    m_h_right_pair_ptbl_stop.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__right_pair_ptbl_stop" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - ptbl" // title suffix
                                                 + "; p_{T}^{bl}^{#tilde{t}} [GeV]" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , ptbl_bins, ptbl_min, ptbl_max
                                               )
                                     );

    m_h_right_pair_ptbl_astp.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__right_pair_ptbl_astp" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - ptbl" // title suffix
                                                 + "; p_{T}^{bl}^{#tilde{t}*} [GeV]" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , ptbl_bins, ptbl_min, ptbl_max
                                               )
                                     );

    m_h_right_pair_ptbl_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__right_pair_ptbl_diff" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - ptbl" // title suffix
                                                 + "; |p_{T}^{bl}^{#tilde{t}} - p_{T}^{bl}^{#tilde{t}*}| [GeV]" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , ptbl_bins, ptbl_min, ptbl_max
                                               )
                                     );

    m_h_right_pair_ptbl_ratio.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__right_pair_ptbl_ratio" // name suffix
                                                  ).c_str()
                                                , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - ptbl" // title suffix
                                                  + "; p_{T}^{bl}^{1}/p_{T}^{bl}^{0}" // x-axis label
                                                  + "; Entries" // y-axis label
                                                  ).c_str()
                                                , ptbl_ratio_bins, ptbl_ratio_min, ptbl_ratio_max
                                                )
                                      );

    m_h_right_pair_ptbl_sq_sum.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + "__right_pair_ptbl_sq_sum" // name suffix
                                                   ).c_str()
                                                 , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + " - ptbl" // title suffix
                                                   + "; #sqrt{p_{T}^{bl}^{#tilde{t} 2} + p_{T}^{bl}^{#tilde{t}* 2}} [GeV]" // x-axis label
                                                   + "; Entries" // y-axis label
                                                   ).c_str()
                                                 , ptbl_sq_sum_bins, ptbl_sq_sum_min, ptbl_sq_sum_max
                                                 )
                                       );

    m_h_right_pair_ptbl_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__right_pair_ptbl_2d" // name suffix
                                               ).c_str()
                                             , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + " - ptbl" // title suffix
                                               + "; p_{T}^{bl}^{#tilde{t}} [GeV]" // x-axis label
                                               + "; #p_{T}^{bl}^{#tilde{t}*} [GeV]" // y-axis label
                                               ).c_str()
                                             , ptbl_bins, ptbl_min, ptbl_max
                                             , ptbl_bins, ptbl_min, ptbl_max
                                             )
                                   );


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_pair_ptbl_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__wrong_pair_ptbl_all" // name suffix
                                                ).c_str()
                                              , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + " - ptbl" // title suffix
                                                + "; p_{T}^{bl} [GeV]" // x-axis label
                                                + "; Entries" // y-axis label
                                                ).c_str()
                                              , ptbl_bins, ptbl_min, ptbl_max
                                              )
                                    );

    m_h_wrong_pair_ptbl_0.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + "__wrong_pair_ptbl_0" // name suffix
                                              ).c_str()
                                            , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + " - ptbl" // title suffix
                                              + "; p_{T}^{bl}^{0} [GeV]" // x-axis label
                                              + "; Entries" // y-axis label
                                              ).c_str()
                                            , ptbl_bins, ptbl_min, ptbl_max
                                            )
                                  );

    m_h_wrong_pair_ptbl_1.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + "__wrong_pair_ptbl_1" // name suffix
                                              ).c_str()
                                            , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + " - ptbl" // title suffix
                                              + "; p_{T}^{bl}^{1} [GeV]" // x-axis label
                                              + "; Entries" // y-axis label
                                              ).c_str()
                                            , ptbl_bins, ptbl_min, ptbl_max
                                            )
                                  );

    m_h_wrong_pair_ptbl_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__wrong_pair_ptbl_diff" // name suffix
                                                 ).c_str()
                                               , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - ptbl" // title suffix
                                                 + "; |p_{T}^{bl}^{0} - p_{T}^{bl}^{1}| [GeV]" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , ptbl_bins, ptbl_min, ptbl_max
                                               )
                                     );

    m_h_wrong_pair_ptbl_ratio.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__wrong_pair_ptbl_ratio" // name suffix
                                                  ).c_str()
                                                , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - ptbl" // title suffix
                                                  + "; p_{T}^{bl}^{1}/p_{T}^{bl}^{0}" // x-axis label
                                                  + "; Entries" // y-axis label
                                                  ).c_str()
                                                , ptbl_ratio_bins, ptbl_ratio_min, ptbl_ratio_max
                                                )
                                      );
    m_h_wrong_pair_ptbl_sq_sum.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + "__wrong_pair_ptbl_sq_sum" // name suffix
                                                   ).c_str()
                                                 , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + " - ptbl" // title suffix
                                                   + "; #sqrt{p_{T}^{bl}^{0 2} + p_{T}^{bl}^{1 2}} [GeV]" // x-axis label
                                                   + "; Entries" // y-axis label
                                                   ).c_str()
                                                 , ptbl_sq_sum_bins, ptbl_sq_sum_min, ptbl_sq_sum_max
                                                 )
                                       );


    m_h_wrong_pair_ptbl_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__wrong_pair_ptbl_2d" // name suffix
                                               ).c_str()
                                             , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + " - ptbl" // title suffix
                                               + "; p_{T}^{bl}^{0} [GeV]" // x-axis label
                                               + "; #p_{T}^{bl}^{1} [GeV]" // y-axis label
                                               ).c_str()
                                             , ptbl_bins, ptbl_min, ptbl_max
                                             , ptbl_bins, ptbl_min, ptbl_max
                                             )
                                   );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_vs_right_ptbl_diff.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                     + "__wrong_vs_right_ptbl_diff" // name suffix
                                                     ).c_str()
                                                   , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                     + " - ptbl" // title suffix
                                                     + "; |p_{T}^{bl}^{#tilde{t}} - p_{T}^{bl}^{#tilde{t}*}|(right pairing) [GeV]" // x-axis label
                                                     + "; |p_{T}^{bl}^{0} - p_{T}^{bl}^{1}|(wrong pairing) [GeV]" // y-axis label
                                                     ).c_str()
                                                   , ptbl_bins, ptbl_min, ptbl_max
                                                   , ptbl_bins, ptbl_min, ptbl_max
                                                   )
                                         );

    m_h_wrong_vs_right_ptbl_ratio.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                      + "__wrong_vs_right_ptbl_ratio" // name suffix
                                                      ).c_str()
                                                    , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                      + " - ptbl" // title suffix
                                                      + "; p_{T}^{bl}^{1}/p_{T}^{bl}^{0} (right pairing)" // x-axis label
                                                      + "; p_{T}^{bl}^{1}/p_{T}^{bl}^{0} (wrong pairing)" // y-axis label
                                                      ).c_str()
                                                    , ptbl_ratio_bins, ptbl_ratio_min, ptbl_ratio_max
                                                    , ptbl_ratio_bins, ptbl_ratio_min, ptbl_ratio_max
                                                    )
                                          );

    m_h_wrong_vs_right_ptbl_sq_sum.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                       + "__wrong_vs_right_ptbl_sq_sum" // name suffix
                                                       ).c_str()
                                                     , ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                       + " - ptbl" // title suffix
                                                       + "; #sqrt{p_{T}^{bl}^{#tilde{t} 2} + p_{T}^{bl}^{#tilde{t}* 2}}(right pairing) [GeV]" // x-axis label
                                                       + "; #sqrt{p_{T}^{bl}^{0 2} + p_{T}^{bl}^{1 2}}(wrong pairing) [GeV]" // y-axis label
                                                       ).c_str()
                                                     , ptbl_sq_sum_bins, ptbl_sq_sum_min, ptbl_sq_sum_max
                                                     , ptbl_sq_sum_bins, ptbl_sq_sum_min, ptbl_sq_sum_max
                                                     )
                                           );


  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::BLPairKinematics::FillSpecial( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Particle*>& el_list
         , const std::vector<TruthNtuple::Particle*>& mu_list
         , const std::vector<TruthNtuple::Particle*>& quark_list
         )
{
  if (flavor_channel == TruthNtuple::FLAVOR_NONE) return;

  // sort objects based on parent particle - if sorting fails, exis the function
  if (sortObjects(el_list, mu_list, quark_list) == false) return;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // compute mbl for correct pairing
  double mbl_stop = TruthNtuple::invariantMass(m_l_from_stop, m_b_from_stop)/1.e3;
  double mbl_astp = TruthNtuple::invariantMass(m_l_from_astp, m_b_from_astp)/1.e3;

  double mbl_right_diff   = fabs(mbl_stop - mbl_astp);
  double mbl_right_ratio  = ( mbl_stop <= mbl_astp
                            ? mbl_stop/mbl_astp
                            : mbl_astp/mbl_stop
                            );
  double mbl_right_sq_sum = sqrt(mbl_stop*mbl_stop + mbl_astp*mbl_astp);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // compute mbl for incorrect pairing
  double mbl_wrong_0 = TruthNtuple::invariantMass( m_l_from_stop
                                                 , m_b_from_astp
                                                 )/1.e3;
  double mbl_wrong_1 = TruthNtuple::invariantMass( m_l_from_astp
                                                 , m_b_from_stop
                                                 )/1.e3;
  if (mbl_wrong_0 < mbl_wrong_1) std::swap(mbl_wrong_0, mbl_wrong_1);

  double mbl_wrong_diff   = fabs(mbl_wrong_0 - mbl_wrong_1);
  double mbl_wrong_ratio  = mbl_wrong_1 / mbl_wrong_0;
  double mbl_wrong_sq_sum = sqrt(mbl_wrong_0*mbl_wrong_0 + mbl_wrong_1*mbl_wrong_1);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // compute mbl for difference pairing
  std::vector<std::pair<TruthNtuple::Particle*, TruthNtuple::Particle*> >
      diff_pairs = BMinusL::doMblDiffPairing(m_b_list, m_l_list);
  double mbl_diff_0 = TruthNtuple::invariantMass( diff_pairs.at(0).first
                                                , diff_pairs.at(0).second
                                                )/1.e3;
  double mbl_diff_1 = TruthNtuple::invariantMass( diff_pairs.at(1).first
                                                , diff_pairs.at(1).second
                                                )/1.e3;
  if (mbl_diff_0 < mbl_diff_1) std::swap(mbl_diff_0, mbl_diff_1);

  double mbl_diff_diff   = fabs(mbl_diff_0 - mbl_diff_1);
  double mbl_diff_ratio  = mbl_diff_1 / mbl_diff_0;
  double mbl_diff_sq_sum = sqrt(mbl_diff_0*mbl_diff_0 + mbl_diff_1*mbl_diff_1);
  bool diff_pariing = (  (  diff_pairs.at(0).first  == m_b_from_stop
                         && diff_pairs.at(0).second == m_l_from_stop
                         && diff_pairs.at(1).first  == m_b_from_astp
                         && diff_pairs.at(1).second == m_l_from_astp
                         )
                      || (  diff_pairs.at(0).first  == m_b_from_astp
                         && diff_pairs.at(0).second == m_l_from_astp
                         && diff_pairs.at(1).first  == m_b_from_stop
                         && diff_pairs.at(1).second == m_l_from_stop
                         )
                      );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // compute mbl for ratio pairing
  std::vector<std::pair<TruthNtuple::Particle*, TruthNtuple::Particle*> >
      ratio_pairs = BMinusL::doMblRatioPairing(m_b_list, m_l_list);
  double mbl_ratio_0 = TruthNtuple::invariantMass( ratio_pairs.at(0).first
                                                 , ratio_pairs.at(0).second
                                                 )/1.e3;
  double mbl_ratio_1 = TruthNtuple::invariantMass( ratio_pairs.at(1).first
                                                 , ratio_pairs.at(1).second
                                                 )/1.e3;
  if (mbl_ratio_0 < mbl_ratio_1) std::swap(mbl_ratio_0, mbl_ratio_1);

  double mbl_ratio_diff   = fabs(mbl_ratio_0 - mbl_ratio_1);
  double mbl_ratio_ratio  = mbl_ratio_1 / mbl_ratio_0;
  double mbl_ratio_sq_sum = sqrt(mbl_ratio_0*mbl_ratio_0 + mbl_ratio_1*mbl_ratio_1);
  bool ratio_pariing = (  (  ratio_pairs.at(0).first  == m_b_from_stop
                          && ratio_pairs.at(0).second == m_l_from_stop
                          && ratio_pairs.at(1).first  == m_b_from_astp
                          && ratio_pairs.at(1).second == m_l_from_astp
                          )
                       || (  ratio_pairs.at(0).first  == m_b_from_astp
                          && ratio_pairs.at(0).second == m_l_from_astp
                          && ratio_pairs.at(1).first  == m_b_from_stop
                          && ratio_pairs.at(1).second == m_l_from_stop
                          )
                       );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // compute mbl for squared sum pairing
  std::vector<std::pair<TruthNtuple::Particle*, TruthNtuple::Particle*> >
      sq_sum_pairs = BMinusL::doMblSqSumPairing(m_b_list, m_l_list);
  double mbl_sq_sum_0 = TruthNtuple::invariantMass( sq_sum_pairs.at(0).first
                                                  , sq_sum_pairs.at(0).second
                                                  )/1.e3;
  double mbl_sq_sum_1 = TruthNtuple::invariantMass( sq_sum_pairs.at(1).first
                                                  , sq_sum_pairs.at(1).second
                                                  )/1.e3;
  if (mbl_sq_sum_0 < mbl_sq_sum_1) std::swap(mbl_sq_sum_0, mbl_sq_sum_1);

  double mbl_sq_sum_diff   = fabs(mbl_sq_sum_0 - mbl_sq_sum_1);
  double mbl_sq_sum_ratio  = mbl_sq_sum_1 / mbl_sq_sum_0;
  double mbl_sq_sum_sq_sum = sqrt(mbl_sq_sum_0*mbl_sq_sum_0 + mbl_sq_sum_1*mbl_sq_sum_1);
  bool sq_sum_pariing = (  (  sq_sum_pairs.at(0).first  == m_b_from_stop
                           && sq_sum_pairs.at(0).second == m_l_from_stop
                           && sq_sum_pairs.at(1).first  == m_b_from_astp
                           && sq_sum_pairs.at(1).second == m_l_from_astp
                           )
                        || (  sq_sum_pairs.at(0).first  == m_b_from_astp
                           && sq_sum_pairs.at(0).second == m_l_from_astp
                           && sq_sum_pairs.at(1).first  == m_b_from_stop
                           && sq_sum_pairs.at(1).second == m_l_from_stop
                           )
                        );


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // compute ptbl for correct pairing
  double ptbl_stop = TruthNtuple::ptDiObject(m_l_from_stop, m_b_from_stop)/1.e3;
  double ptbl_astp = TruthNtuple::ptDiObject(m_l_from_astp, m_b_from_astp)/1.e3;

  double right_ptbl_diff   = fabs(ptbl_stop - ptbl_astp);
  double right_ptbl_ratio  = ( ptbl_stop <= ptbl_astp
                            ? ptbl_stop/ptbl_astp
                            : ptbl_astp/ptbl_stop
                            );
  double right_ptbl_sq_sum = sqrt(ptbl_stop*ptbl_stop + ptbl_astp*ptbl_astp);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // compute ptbl for incorrect pairing
  double ptbl_wrong_0 = TruthNtuple::ptDiObject(m_l_from_stop, m_b_from_astp)/1.e3;
  double ptbl_wrong_1 = TruthNtuple::ptDiObject(m_l_from_astp, m_b_from_stop)/1.e3;
  if (ptbl_wrong_0 < ptbl_wrong_1) std::swap(ptbl_wrong_0, ptbl_wrong_1);

  double wrong_ptbl_diff   = fabs(ptbl_wrong_0 - ptbl_wrong_1);
  double wrong_ptbl_ratio  = ptbl_wrong_1 / ptbl_wrong_0;
  double wrong_ptbl_sq_sum = sqrt(ptbl_wrong_0*ptbl_wrong_0 - ptbl_wrong_1*ptbl_wrong_1);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // loop through flavor channels
  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // fill lepton kinematic plots
      m_h_l_pt_all.at(fc)->Fill(m_l_from_stop->getPt()/1.e3);
      m_h_l_pt_all.at(fc)->Fill(m_l_from_astp->getPt()/1.e3);
      m_h_l_pt_stop.at(fc)->Fill(m_l_from_stop->getPt()/1.e3);
      m_h_l_pt_astp.at(fc)->Fill(m_l_from_astp->getPt()/1.e3);
      m_h_l_pt_diff.at(fc)->Fill(fabs( m_l_from_stop->getPt()
                                     - m_l_from_astp->getPt()
                                     )/1.e3
                                );
      m_h_l_pt_2d.at(fc)->Fill( m_l_from_stop->getPt()/1.e3
                              , m_l_from_astp->getPt()/1.e3
                              );

      m_h_l_eta_all.at(fc)->Fill(m_l_from_stop->getEta());
      m_h_l_eta_all.at(fc)->Fill(m_l_from_astp->getEta());
      m_h_l_eta_stop.at(fc)->Fill(m_l_from_stop->getEta());
      m_h_l_eta_astp.at(fc)->Fill(m_l_from_astp->getEta());
      m_h_l_eta_diff.at(fc)->Fill(TruthNtuple::deltaEta( m_l_from_stop
                                                       , m_l_from_astp
                                                       )
                                 );
      m_h_l_eta_2d.at(fc)->Fill( m_l_from_stop->getEta()
                               , m_l_from_astp->getEta()
                               );

      m_h_l_phi_all.at(fc)->Fill(m_l_from_stop->getPhi());
      m_h_l_phi_all.at(fc)->Fill(m_l_from_astp->getPhi());
      m_h_l_phi_stop.at(fc)->Fill(m_l_from_stop->getPhi());
      m_h_l_phi_astp.at(fc)->Fill(m_l_from_astp->getPhi());
      m_h_l_phi_diff.at(fc)->Fill(TruthNtuple::deltaPhi( m_l_from_stop
                                                       , m_l_from_astp
                                                       )
                                 );
      m_h_l_phi_2d.at(fc)->Fill( m_l_from_stop->getPhi()
                               , m_l_from_astp->getPhi()
                               );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // fill b quark kinematic plots
      m_h_b_pt_all.at(fc)->Fill(m_b_from_stop->getPt()/1.e3);
      m_h_b_pt_all.at(fc)->Fill(m_b_from_astp->getPt()/1.e3);
      m_h_b_pt_stop.at(fc)->Fill(m_b_from_stop->getPt()/1.e3);
      m_h_b_pt_astp.at(fc)->Fill(m_b_from_astp->getPt()/1.e3);
      m_h_b_pt_diff.at(fc)->Fill(fabs( m_b_from_stop->getPt()
                                     - m_b_from_astp->getPt()
                                     )/1.e3
                                );
      m_h_b_pt_2d.at(fc)->Fill( m_b_from_stop->getPt()/1.e3
                              , m_b_from_astp->getPt()/1.e3
                              );

      m_h_b_eta_all.at(fc)->Fill(m_b_from_stop->getEta());
      m_h_b_eta_all.at(fc)->Fill(m_b_from_astp->getEta());
      m_h_b_eta_stop.at(fc)->Fill(m_b_from_stop->getEta());
      m_h_b_eta_astp.at(fc)->Fill(m_b_from_astp->getEta());
      m_h_b_eta_diff.at(fc)->Fill(TruthNtuple::deltaEta( m_b_from_stop
                                                       , m_b_from_astp
                                                       )
                                 );
      m_h_b_eta_2d.at(fc)->Fill( m_b_from_stop->getEta()
                               , m_b_from_astp->getEta()
                               );

      m_h_b_phi_all.at(fc)->Fill(m_b_from_stop->getPhi());
      m_h_b_phi_all.at(fc)->Fill(m_b_from_astp->getPhi());
      m_h_b_phi_stop.at(fc)->Fill(m_b_from_stop->getPhi());
      m_h_b_phi_astp.at(fc)->Fill(m_b_from_astp->getPhi());
      m_h_b_phi_diff.at(fc)->Fill(TruthNtuple::deltaPhi( m_b_from_stop
                                                       , m_b_from_astp
                                                       )
                                 );
      m_h_b_phi_2d.at(fc)->Fill( m_b_from_stop->getPhi()
                               , m_b_from_astp->getPhi()
                               );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // fill b-vs-l histograms for right pairing
      m_h_right_pair_bl_dpt.at(fc)->Fill( m_l_from_stop->getPt()/1.e3
                                        - m_b_from_stop->getPt()/1.e3
                                        );
      m_h_right_pair_bl_deta.at(fc)->Fill( TruthNtuple::deltaEta( m_l_from_stop
                                                                , m_b_from_stop
                                                                )
                                         );
      m_h_right_pair_bl_dphi.at(fc)->Fill( TruthNtuple::deltaPhi( m_l_from_stop
                                                                , m_b_from_stop
                                                                )
                                         );
      m_h_right_pair_bl_dr.at(fc)->Fill( TruthNtuple::deltaR( m_l_from_stop
                                                            , m_b_from_stop
                                                            )
                                       );

      m_h_right_pair_bl_dpt.at(fc)->Fill( m_l_from_astp->getPt()/1.e3
                                        - m_b_from_astp->getPt()/1.e3
                                        );
      m_h_right_pair_bl_deta.at(fc)->Fill( TruthNtuple::deltaEta( m_l_from_astp
                                                                , m_b_from_astp
                                                                )
                                         );
      m_h_right_pair_bl_dphi.at(fc)->Fill( TruthNtuple::deltaPhi( m_l_from_astp
                                                                , m_b_from_astp
                                                                )
                                         );
      m_h_right_pair_bl_dr.at(fc)->Fill( TruthNtuple::deltaR( m_l_from_astp
                                                            , m_b_from_astp
                                                            )
                                       );

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
      m_h_wrong_pair_bl_dpt.at(fc)->Fill( m_l_from_stop->getPt()/1.e3
                                        - m_b_from_astp->getPt()/1.e3
                                        );
      m_h_wrong_pair_bl_deta.at(fc)->Fill( TruthNtuple::deltaEta( m_l_from_stop
                                                                , m_b_from_astp
                                                                )
                                         );
      m_h_wrong_pair_bl_dphi.at(fc)->Fill( TruthNtuple::deltaPhi( m_l_from_stop
                                                                , m_b_from_astp
                                                                )
                                         );
      m_h_wrong_pair_bl_dr.at(fc)->Fill( TruthNtuple::deltaR( m_l_from_stop
                                                            , m_b_from_astp
                                                            )
                                       );

      m_h_wrong_pair_bl_dpt.at(fc)->Fill( m_l_from_astp->getPt()/1.e3
                                        - m_b_from_stop->getPt()/1.e3
                                        );
      m_h_wrong_pair_bl_deta.at(fc)->Fill( TruthNtuple::deltaEta( m_l_from_astp
                                                                , m_b_from_stop
                                                                )
                                         );
      m_h_wrong_pair_bl_dphi.at(fc)->Fill( TruthNtuple::deltaPhi( m_l_from_astp
                                                                , m_b_from_stop
                                                                )
                                         );
      m_h_wrong_pair_bl_dr.at(fc)->Fill( TruthNtuple::deltaR( m_l_from_astp
                                                            , m_b_from_stop
                                                            )
                                       );

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
      m_h_right_pair_mbl_diff.at(fc)->Fill(mbl_right_diff);
      m_h_right_pair_mbl_ratio.at(fc)->Fill(mbl_right_ratio);
      m_h_right_pair_mbl_sq_sum.at(fc)->Fill(mbl_right_sq_sum);
      m_h_right_pair_mbl_2d.at(fc)->Fill(mbl_stop, mbl_astp);

      // fill mbl histograms for wrong pairing
      m_h_wrong_pair_mbl_all.at(fc)->Fill(mbl_wrong_0);
      m_h_wrong_pair_mbl_all.at(fc)->Fill(mbl_wrong_1);
      m_h_wrong_pair_mbl_0.at(fc)->Fill(mbl_wrong_0);
      m_h_wrong_pair_mbl_1.at(fc)->Fill(mbl_wrong_1);
      m_h_wrong_pair_mbl_diff.at(fc)->Fill(mbl_wrong_diff);
      m_h_wrong_pair_mbl_ratio.at(fc)->Fill(mbl_wrong_ratio);
      m_h_wrong_pair_mbl_sq_sum.at(fc)->Fill(mbl_wrong_sq_sum);
      m_h_wrong_pair_mbl_2d.at(fc)->Fill(mbl_wrong_0, mbl_wrong_1);

      m_h_diff_pair_mbl_all.at(fc)->Fill(mbl_diff_0);
      m_h_diff_pair_mbl_all.at(fc)->Fill(mbl_diff_1);
      m_h_diff_pair_mbl_0.at(fc)->Fill(mbl_diff_0);
      m_h_diff_pair_mbl_1.at(fc)->Fill(mbl_diff_1);
      m_h_diff_pair_mbl_diff.at(fc)->Fill(mbl_diff_diff);
      m_h_diff_pair_mbl_ratio.at(fc)->Fill(mbl_diff_ratio);
      m_h_diff_pair_mbl_sq_sum.at(fc)->Fill(mbl_diff_sq_sum);
      m_h_diff_pair_mbl_2d.at(fc)->Fill(mbl_diff_0, mbl_diff_1);
      m_h_diff_pair_cor_pairing.at(fc)->Fill(diff_pariing);

      m_h_ratio_pair_mbl_all.at(fc)->Fill(mbl_ratio_0);
      m_h_ratio_pair_mbl_all.at(fc)->Fill(mbl_ratio_1);
      m_h_ratio_pair_mbl_0.at(fc)->Fill(mbl_ratio_0);
      m_h_ratio_pair_mbl_1.at(fc)->Fill(mbl_ratio_1);
      m_h_ratio_pair_mbl_diff.at(fc)->Fill(mbl_ratio_diff);
      m_h_ratio_pair_mbl_ratio.at(fc)->Fill(mbl_ratio_ratio);
      m_h_ratio_pair_mbl_sq_sum.at(fc)->Fill(mbl_ratio_sq_sum);
      m_h_ratio_pair_mbl_2d.at(fc)->Fill(mbl_ratio_0, mbl_ratio_1);
      m_h_ratio_pair_cor_pairing.at(fc)->Fill(ratio_pariing);

      m_h_sq_sum_pair_mbl_all.at(fc)->Fill(mbl_sq_sum_0);
      m_h_sq_sum_pair_mbl_all.at(fc)->Fill(mbl_sq_sum_1);
      m_h_sq_sum_pair_mbl_0.at(fc)->Fill(mbl_sq_sum_0);
      m_h_sq_sum_pair_mbl_1.at(fc)->Fill(mbl_sq_sum_1);
      m_h_sq_sum_pair_mbl_diff.at(fc)->Fill(mbl_sq_sum_diff);
      m_h_sq_sum_pair_mbl_ratio.at(fc)->Fill(mbl_sq_sum_ratio);
      m_h_sq_sum_pair_mbl_sq_sum.at(fc)->Fill(mbl_sq_sum_sq_sum);
      m_h_sq_sum_pair_mbl_2d.at(fc)->Fill(mbl_sq_sum_0, mbl_sq_sum_1);
      m_h_sq_sum_pair_cor_pairing.at(fc)->Fill(sq_sum_pariing);

      // fill mbl histograms wrong vs right pairing
      m_h_wrong_vs_right_mbl_diff.at(fc)->Fill(  mbl_right_diff  , mbl_wrong_diff  );
      m_h_wrong_vs_right_mbl_ratio.at(fc)->Fill( mbl_right_ratio , mbl_wrong_ratio );
      m_h_wrong_vs_right_mbl_sq_sum.at(fc)->Fill(mbl_right_sq_sum, mbl_wrong_sq_sum);

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // fill ptbl histograms for right pairing
      m_h_right_pair_ptbl_all.at(fc)->Fill(ptbl_stop);
      m_h_right_pair_ptbl_all.at(fc)->Fill(ptbl_astp);
      m_h_right_pair_ptbl_stop.at(fc)->Fill(ptbl_stop);
      m_h_right_pair_ptbl_astp.at(fc)->Fill(ptbl_astp);
      m_h_right_pair_ptbl_diff.at(fc)->Fill(right_ptbl_diff);
      m_h_right_pair_ptbl_ratio.at(fc)->Fill(right_ptbl_ratio);
      m_h_right_pair_ptbl_sq_sum.at(fc)->Fill(right_ptbl_sq_sum);
      m_h_right_pair_ptbl_2d.at(fc)->Fill(ptbl_stop, ptbl_astp);

      // fill ptbl histograms for wrong pairing
      m_h_wrong_pair_ptbl_all.at(fc)->Fill(ptbl_wrong_0);
      m_h_wrong_pair_ptbl_all.at(fc)->Fill(ptbl_wrong_1);
      m_h_wrong_pair_ptbl_0.at(fc)->Fill(ptbl_wrong_0);
      m_h_wrong_pair_ptbl_1.at(fc)->Fill(ptbl_wrong_1);
      m_h_wrong_pair_ptbl_diff.at(fc)->Fill(wrong_ptbl_diff);
      m_h_wrong_pair_ptbl_ratio.at(fc)->Fill(wrong_ptbl_ratio);
      m_h_wrong_pair_ptbl_sq_sum.at(fc)->Fill(wrong_ptbl_sq_sum);
      m_h_wrong_pair_ptbl_2d.at(fc)->Fill(ptbl_wrong_0, ptbl_wrong_1);

      // fill ptbl histograms wrong vs right pairing
      m_h_wrong_vs_right_ptbl_diff.at(fc)->Fill(  right_ptbl_diff  , wrong_ptbl_diff  );
      m_h_wrong_vs_right_ptbl_ratio.at(fc)->Fill( right_ptbl_ratio , wrong_ptbl_ratio );
      m_h_wrong_vs_right_ptbl_sq_sum.at(fc)->Fill(right_ptbl_sq_sum, wrong_ptbl_sq_sum);
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::BLPairKinematics::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    m_h_l_pt_all.at(fc_it)->Write();
    m_h_l_pt_stop.at(fc_it)->Write();
    m_h_l_pt_astp.at(fc_it)->Write();
    m_h_l_pt_diff.at(fc_it)->Write();
    m_h_l_pt_2d.at(fc_it)->Write();

    m_h_l_eta_all.at(fc_it)->Write();
    m_h_l_eta_stop.at(fc_it)->Write();
    m_h_l_eta_astp.at(fc_it)->Write();
    m_h_l_eta_diff.at(fc_it)->Write();
    m_h_l_eta_2d.at(fc_it)->Write();

    m_h_l_phi_all.at(fc_it)->Write();
    m_h_l_phi_stop.at(fc_it)->Write();
    m_h_l_phi_astp.at(fc_it)->Write();
    m_h_l_phi_diff.at(fc_it)->Write();
    m_h_l_phi_2d.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_b_pt_all.at(fc_it)->Write();
    m_h_b_pt_stop.at(fc_it)->Write();
    m_h_b_pt_astp.at(fc_it)->Write();
    m_h_b_pt_diff.at(fc_it)->Write();
    m_h_b_pt_2d.at(fc_it)->Write();

    m_h_b_eta_all.at(fc_it)->Write();
    m_h_b_eta_stop.at(fc_it)->Write();
    m_h_b_eta_astp.at(fc_it)->Write();
    m_h_b_eta_diff.at(fc_it)->Write();
    m_h_b_eta_2d.at(fc_it)->Write();

    m_h_b_phi_all.at(fc_it)->Write();
    m_h_b_phi_stop.at(fc_it)->Write();
    m_h_b_phi_astp.at(fc_it)->Write();
    m_h_b_phi_diff.at(fc_it)->Write();
    m_h_b_phi_2d.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_right_pair_bl_dpt.at(fc_it)->Write();
    m_h_right_pair_bl_deta.at(fc_it)->Write();
    m_h_right_pair_bl_dphi.at(fc_it)->Write();
    m_h_right_pair_bl_dr.at(fc_it)->Write();

    m_h_right_pair_bl_pt_2d.at(fc_it)->Write();
    m_h_right_pair_bl_eta_2d.at(fc_it)->Write();
    m_h_right_pair_bl_phi_2d.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_pair_bl_dpt.at(fc_it)->Write();
    m_h_wrong_pair_bl_deta.at(fc_it)->Write();
    m_h_wrong_pair_bl_dphi.at(fc_it)->Write();
    m_h_wrong_pair_bl_dr.at(fc_it)->Write();

    m_h_wrong_pair_bl_pt_2d.at(fc_it)->Write();
    m_h_wrong_pair_bl_eta_2d.at(fc_it)->Write();
    m_h_wrong_pair_bl_phi_2d.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_right_pair_mbl_all.at(fc_it)->Write();
    m_h_right_pair_mbl_stop.at(fc_it)->Write();
    m_h_right_pair_mbl_astp.at(fc_it)->Write();
    m_h_right_pair_mbl_diff.at(fc_it)->Write();
    m_h_right_pair_mbl_ratio.at(fc_it)->Write();
    m_h_right_pair_mbl_sq_sum.at(fc_it)->Write();
    m_h_right_pair_mbl_2d.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_pair_mbl_all.at(fc_it)->Write();
    m_h_wrong_pair_mbl_0.at(fc_it)->Write();
    m_h_wrong_pair_mbl_1.at(fc_it)->Write();
    m_h_wrong_pair_mbl_diff.at(fc_it)->Write();
    m_h_wrong_pair_mbl_ratio.at(fc_it)->Write();
    m_h_wrong_pair_mbl_sq_sum.at(fc_it)->Write();
    m_h_wrong_pair_mbl_2d.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_diff_pair_mbl_all.at(fc_it)->Write();
    m_h_diff_pair_mbl_0.at(fc_it)->Write();
    m_h_diff_pair_mbl_1.at(fc_it)->Write();
    m_h_diff_pair_mbl_diff.at(fc_it)->Write();
    m_h_diff_pair_mbl_ratio.at(fc_it)->Write();
    m_h_diff_pair_mbl_sq_sum.at(fc_it)->Write();
    m_h_diff_pair_mbl_2d.at(fc_it)->Write();
    m_h_diff_pair_cor_pairing.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_ratio_pair_mbl_all.at(fc_it)->Write();
    m_h_ratio_pair_mbl_0.at(fc_it)->Write();
    m_h_ratio_pair_mbl_1.at(fc_it)->Write();
    m_h_ratio_pair_mbl_diff.at(fc_it)->Write();
    m_h_ratio_pair_mbl_ratio.at(fc_it)->Write();
    m_h_ratio_pair_mbl_sq_sum.at(fc_it)->Write();
    m_h_ratio_pair_mbl_2d.at(fc_it)->Write();
    m_h_ratio_pair_cor_pairing.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_sq_sum_pair_mbl_all.at(fc_it)->Write();
    m_h_sq_sum_pair_mbl_0.at(fc_it)->Write();
    m_h_sq_sum_pair_mbl_1.at(fc_it)->Write();
    m_h_sq_sum_pair_mbl_diff.at(fc_it)->Write();
    m_h_sq_sum_pair_mbl_ratio.at(fc_it)->Write();
    m_h_sq_sum_pair_mbl_sq_sum.at(fc_it)->Write();
    m_h_sq_sum_pair_mbl_2d.at(fc_it)->Write();
    m_h_sq_sum_pair_cor_pairing.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_vs_right_mbl_diff.at(fc_it)->Write();
    m_h_wrong_vs_right_mbl_ratio.at(fc_it)->Write();
    m_h_wrong_vs_right_mbl_sq_sum.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_right_pair_ptbl_all.at(fc_it)->Write();
    m_h_right_pair_ptbl_stop.at(fc_it)->Write();
    m_h_right_pair_ptbl_astp.at(fc_it)->Write();
    m_h_right_pair_ptbl_diff.at(fc_it)->Write();
    m_h_right_pair_ptbl_ratio.at(fc_it)->Write();
    m_h_right_pair_ptbl_sq_sum.at(fc_it)->Write();
    m_h_right_pair_ptbl_2d.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_pair_ptbl_all.at(fc_it)->Write();
    m_h_wrong_pair_ptbl_0.at(fc_it)->Write();
    m_h_wrong_pair_ptbl_1.at(fc_it)->Write();
    m_h_wrong_pair_ptbl_diff.at(fc_it)->Write();
    m_h_wrong_pair_ptbl_ratio.at(fc_it)->Write();
    m_h_wrong_pair_ptbl_sq_sum.at(fc_it)->Write();
    m_h_wrong_pair_ptbl_2d.at(fc_it)->Write();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    m_h_wrong_vs_right_ptbl_diff.at(fc_it)->Write();
    m_h_wrong_vs_right_ptbl_ratio.at(fc_it)->Write();
    m_h_wrong_vs_right_ptbl_sq_sum.at(fc_it)->Write();
  }
}

// -----------------------------------------------------------------------------
bool HistogramHandlers::BLPairKinematics::sortObjects( const std::vector<TruthNtuple::Particle*>& el_list
                                                     , const std::vector<TruthNtuple::Particle*>& mu_list
                                                     , const std::vector<TruthNtuple::Particle*>& quark_list
                                                     )
{
  // merge el and mu lists to lepton list
  m_l_list.clear();
  m_l_list.reserve(el_list.size() + mu_list.size());
  m_l_list.insert(m_l_list.end(), el_list.begin(), el_list.end());
  m_l_list.insert(m_l_list.end(), mu_list.begin(), mu_list.end());

  // Find lepton from stop and anti-stop
  m_l_from_stop = 0;
  m_l_from_astp = 0;

  // I'm probably being overly careful here
  for (size_t lep_it = 0; lep_it != m_l_list.size(); ++lep_it) {
    // if parent is stop (or anti-tau)
    // if (m_l_list.at(lep_it)->getParentPdgid() == +(1e6+6))  {
    if (  m_l_list.at(lep_it)->getParentPdgid() == +(1e6+6)
       || m_l_list.at(lep_it)->getParentPdgid() == -15
       )  {
      if (m_l_from_stop == 0)
        m_l_from_stop = m_l_list.at(lep_it);
      else
        std::cout << "WARNING! Found multiple leptons paired to stop!\n";
    }
    // if parent is anti-stop (or tau)
    // if (m_l_list.at(lep_it)->getParentPdgid() == -(1e6+6)) {
    if (  m_l_list.at(lep_it)->getParentPdgid() == -(1e6+6)
       || m_l_list.at(lep_it)->getParentPdgid() == +15
       ) {
      if (m_l_from_astp == 0)
        m_l_from_astp = m_l_list.at(lep_it);
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
  m_b_list = quark_list;
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
  if (  m_l_from_stop->getParentBarcode() != m_b_from_stop->getParentBarcode()
     && m_l_from_stop->getParentPdgid() != -15
     ) {
    std::cout << "ERROR! Lepton and b from stop have different parent barcodes\n";
    std::cout << "lep parent: " << m_l_from_stop->getParentPdgid() << "\n";
    return false;
  }
  // check that the objects from the anti-stop have the same parent barcode
  if (  m_l_from_astp->getParentBarcode() != m_b_from_astp->getParentBarcode()
     && m_l_from_astp->getParentPdgid() != +15
     ) {
    std::cout << "ERROR! Lepton and b from anti-stop have different parent barcodes\n";
    std::cout << "lep parent: " << m_l_from_astp->getParentPdgid() << "\n";
    return false;
  }

  return true;
}


// =============================================================================
// = Mbl
// =============================================================================
HistogramHandlers::Mbl::Mbl() : HistogramHandlers::Handle()
{
  const int mbl_bins   = 150;
  const double mbl_min = 0;
  const double mbl_max = 1500;

  const int    mbl_ratio_bins = 110;
  const double mbl_ratio_min  = 0.;
  const double mbl_ratio_max  = 1.1;

  const int    mbl_sq_sum_bins = 300;
  const double mbl_sq_sum_min  = 0.;
  const double mbl_sq_sum_max  = 3000.;

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    m_h_ratio_pair_mbl_all.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                + "__mjl__ratio_pair_mbl_all"
                                                ).c_str()
                                              , ( "m_{bl}^{truth} - "
                                                + TruthNtuple::FlavorChannelStrings[fc_it]
                                                + " - mbl" // title suffix
                                                + "; m_{bl} [GeV]" // x-axis label
                                                + "; Entries" // y-axis label
                                                ).c_str()
                                              , mbl_bins, mbl_min, mbl_max
                                              )
                                    );
    m_h_ratio_pair_mbl_0.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + "__mjl__ratio_pair_mbl_0"
                                              ).c_str()
                                            , ( "m_{bl}^{truth} - "
                                              + TruthNtuple::FlavorChannelStrings[fc_it]
                                              + " - mbl" // title suffix
                                              + "; m_{bl}^{0} [GeV]" // x-axis label
                                              + "; Entries" // y-axis label
                                              ).c_str()
                                            , mbl_bins, mbl_min, mbl_max
                                            )
                                  );
    m_h_ratio_pair_mbl_1.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                              + "__mjl__ratio_pair_mbl_1"
                                              ).c_str()
                                            , ( "m_{bl}^{truth} - "
                                              + TruthNtuple::FlavorChannelStrings[fc_it]
                                              + " - mbl" // title suffix
                                              + "; m_{bl}^{1} [GeV]" // x-axis label
                                              + "; Entries" // y-axis label
                                              ).c_str()
                                            , mbl_bins, mbl_min, mbl_max
                                            )
                                  );
    m_h_ratio_pair_mbl_diff.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + "__mjl__ratio_pair_mbl_diff"
                                                 ).c_str()
                                               , ( "m_{bl}^{truth} - "
                                                 + TruthNtuple::FlavorChannelStrings[fc_it]
                                                 + " - mbl" // title suffix
                                                 + "; |m_{bl}^{0} - m_{bl}^{1}| [GeV]" // x-axis label
                                                 + "; Entries" // y-axis label
                                                 ).c_str()
                                               , mbl_bins, mbl_min, mbl_max
                                               )
                                     );
    m_h_ratio_pair_mbl_ratio.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + "__mjl__ratio_pair_mbl_ratio"
                                                  ).c_str()
                                                , ( "m_{bl}^{truth} - "
                                                  + TruthNtuple::FlavorChannelStrings[fc_it]
                                                  + " - mbl" // title suffix
                                                  + "; m_{bl}^{1}/m_{bl}^{0}" // x-axis label
                                                  + "; Entries" // y-axis label
                                                  ).c_str()
                                                , mbl_ratio_bins, mbl_ratio_min, mbl_ratio_max
                                                )
                                      );
    m_h_ratio_pair_mbl_sq_sum.push_back( new TH1D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + "__mjl__ratio_pair_mbl_sq_sum"
                                                   ).c_str()
                                                 , ( "m_{bl}^{truth} - "
                                                   + TruthNtuple::FlavorChannelStrings[fc_it]
                                                   + " - mbl" // title suffix
                                                   + "; #sqrt{m_{bl}^{0 2} + m_{bl}^{1 2}} [GeV]" // x-axis label
                                                   + "; Entries" // y-axis label
                                                   ).c_str()
                                                 , mbl_sq_sum_bins, mbl_sq_sum_min, mbl_sq_sum_max
                                                 )
                                       );
    m_h_ratio_pair_mbl_2d.push_back( new TH2D( ( TruthNtuple::FlavorChannelStrings[fc_it]
                                               + "__mjl__ratio_pair_mbl_2d"
                                               ).c_str()
                                             , ( "m_{bl}^{truth} - "
                                               + TruthNtuple::FlavorChannelStrings[fc_it]
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
void HistogramHandlers::Mbl::FillSpecial( const TruthNtuple::FLAVOR_CHANNEL flavor_channel
         , const std::vector<TruthNtuple::Particle*>& el_list
         , const std::vector<TruthNtuple::Particle*>& mu_list
         , const std::vector<TruthNtuple::Particle*>& b_jet_list
         )
{
  if (flavor_channel == TruthNtuple::FLAVOR_NONE) return;
  if (b_jet_list.size() != 2) return;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // merge el and mu lists to lepton list
  std::vector<TruthNtuple::Particle*> lep_list;
  lep_list.reserve(el_list.size() + mu_list.size());
  lep_list.insert(lep_list.end(), el_list.begin(), el_list.end());
  lep_list.insert(lep_list.end(), mu_list.begin(), mu_list.end());

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // compute mbl for ratio pairing
  std::vector<std::pair<TruthNtuple::Particle*, TruthNtuple::Particle*> >
      ratio_pairs = BMinusL::doMblRatioPairing(b_jet_list, lep_list);
  double mbl_ratio_0 = TruthNtuple::invariantMass( ratio_pairs.at(0).first
                                                 , ratio_pairs.at(0).second
                                                 )/1.e3;
  double mbl_ratio_1 = TruthNtuple::invariantMass( ratio_pairs.at(1).first
                                                 , ratio_pairs.at(1).second
                                                 )/1.e3;
  if (mbl_ratio_0 < mbl_ratio_1) std::swap(mbl_ratio_0, mbl_ratio_1);

  double mbl_ratio_diff   = fabs(mbl_ratio_0 - mbl_ratio_1);
  double mbl_ratio_ratio  = mbl_ratio_1 / mbl_ratio_0;
  double mbl_ratio_sq_sum = sqrt(mbl_ratio_0*mbl_ratio_0 + mbl_ratio_1*mbl_ratio_1);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // loop over flavor channels
  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
    TruthNtuple::FLAVOR_CHANNEL fc = TruthNtuple::FLAVOR_CHANNEL(fc_it);
    if (fc == TruthNtuple::FLAVOR_ALL || fc == flavor_channel) {
      m_h_ratio_pair_mbl_all.at(fc)->Fill(mbl_ratio_0);
      m_h_ratio_pair_mbl_all.at(fc)->Fill(mbl_ratio_1);
      m_h_ratio_pair_mbl_0.at(fc)->Fill(mbl_ratio_0);
      m_h_ratio_pair_mbl_1.at(fc)->Fill(mbl_ratio_1);
      m_h_ratio_pair_mbl_diff.at(fc)->Fill(mbl_ratio_diff);
      m_h_ratio_pair_mbl_ratio.at(fc)->Fill(mbl_ratio_ratio);
      m_h_ratio_pair_mbl_sq_sum.at(fc)->Fill(mbl_ratio_sq_sum);
      m_h_ratio_pair_mbl_2d.at(fc)->Fill(mbl_ratio_0, mbl_ratio_1);
    }
  }
}

// -----------------------------------------------------------------------------
void HistogramHandlers::Mbl::write(TFile* f)
{
  f->cd();

  for (unsigned int fc_it = 0; fc_it != TruthNtuple::FLAVOR_N; ++fc_it) {
      m_h_ratio_pair_mbl_all.at(fc_it)->Write();
      m_h_ratio_pair_mbl_0.at(fc_it)->Write();
      m_h_ratio_pair_mbl_1.at(fc_it)->Write();
      m_h_ratio_pair_mbl_diff.at(fc_it)->Write();
      m_h_ratio_pair_mbl_ratio.at(fc_it)->Write();
      m_h_ratio_pair_mbl_sq_sum.at(fc_it)->Write();
      m_h_ratio_pair_mbl_2d.at(fc_it)->Write();
  }
}
