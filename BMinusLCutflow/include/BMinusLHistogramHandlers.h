#ifndef BMINUSL_HISTOGRAMHANDLERS_H
#define BMINUSL_HISTOGRAMHANDLERS_H

#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TruthNtupleLooper/include/TruthNtupleEnums.h"
#include "HistogramHandlers/include/HistogramHandlers.h"

// =============================================================================
class TFile;

namespace TruthNtuple
{
  class Particle;
  class Lepton;
  class Electron;
  class Muon;
  class Jet;
  class Met;
}

// =============================================================================
namespace HistogramHandlers
{
  // =============================================================================
  // = Stop kinematics
  // =============================================================================
  class StopKinematics : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      StopKinematics();

      virtual void FillSpecial( const TruthNtuple::FLAVOR_CHANNEL
                              , const std::vector<TruthNtuple::Particle*>&
                              );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_pt_all;
      std::vector<TH1F*> m_h_pt_stop; // stop
      std::vector<TH1F*> m_h_pt_astp; // anti-stop
      std::vector<TH1F*> m_h_pt_diff;
      std::vector<TH2F*> m_h_pt_2d;

      std::vector<TH1F*> m_h_eta_all;
      std::vector<TH1F*> m_h_eta_stop; // stop
      std::vector<TH1F*> m_h_eta_astp; // anti-stop
      std::vector<TH1F*> m_h_eta_diff;
      std::vector<TH2F*> m_h_eta_2d;

      std::vector<TH1F*> m_h_phi_all;
      std::vector<TH1F*> m_h_phi_stop; // stop
      std::vector<TH1F*> m_h_phi_astp; // anti-stop
      std::vector<TH1F*> m_h_phi_diff;
      std::vector<TH2F*> m_h_phi_2d;
  };

  // =============================================================================
  // = QuarkKinematics
  // =============================================================================
  class QuarkKinematics : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      QuarkKinematics();

      virtual void FillSpecial( const TruthNtuple::FLAVOR_CHANNEL
                              , const std::vector<TruthNtuple::Particle*>&
                              );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_pt_all;
      std::vector<TH1F*> m_h_pt_0;
      std::vector<TH1F*> m_h_pt_1;
      std::vector<TH1F*> m_h_pt_diff;
      std::vector<TH2F*> m_h_pt_2d;

      std::vector<TH1F*> m_h_eta_all;
      std::vector<TH1F*> m_h_eta_0;
      std::vector<TH1F*> m_h_eta_1;
      std::vector<TH1F*> m_h_eta_diff;
      std::vector<TH2F*> m_h_eta_2d;

      std::vector<TH1F*> m_h_phi_all;
      std::vector<TH1F*> m_h_phi_0;
      std::vector<TH1F*> m_h_phi_1;
      std::vector<TH1F*> m_h_phi_diff;
      std::vector<TH2F*> m_h_phi_2d;

  };

  // =============================================================================
  // = B-l pair kinematics
  // =============================================================================
  class BLPairKinematics : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      BLPairKinematics();

      virtual void FillSpecial( const TruthNtuple::FLAVOR_CHANNEL
                              , const std::vector<TruthNtuple::Electron*>&
                              , const std::vector<TruthNtuple::Muon*>&
                              , const std::vector<TruthNtuple::Particle*>&
                              );
      virtual void write(TFile*);

      bool sortObjects( const std::vector<TruthNtuple::Electron*>&
                      , const std::vector<TruthNtuple::Muon*>&
                      , const std::vector<TruthNtuple::Particle*>&
                      );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      // helper pointers to keep track of leptons from stop and anti-stop
      TruthNtuple::Lepton* m_l_from_stop;
      TruthNtuple::Lepton* m_l_from_astp;

      // helper pointers to keep track of b quarks from stop and anti-stop
      TruthNtuple::Particle* m_b_from_stop;
      TruthNtuple::Particle* m_b_from_astp;

      // 1d plots of lepton kinematics. separate plots for leptons from stop
      // and anti-stop
      std::vector<TH1F*> m_h_l_pt_stop;
      std::vector<TH1F*> m_h_l_pt_astp;
      std::vector<TH1F*> m_h_l_pt_diff;
      std::vector<TH2F*> m_h_l_pt_2d;

      std::vector<TH1F*> m_h_l_eta_stop;
      std::vector<TH1F*> m_h_l_eta_astp;
      std::vector<TH1F*> m_h_l_eta_diff;
      std::vector<TH2F*> m_h_l_eta_2d;

      std::vector<TH1F*> m_h_l_phi_stop;
      std::vector<TH1F*> m_h_l_phi_astp;
      std::vector<TH1F*> m_h_l_phi_diff;
      std::vector<TH2F*> m_h_l_phi_2d;

      // 1d plots of b quark kinematics. separate plots for leptons from stop
      // and anti-stop
      std::vector<TH1F*> m_h_b_pt_stop;
      std::vector<TH1F*> m_h_b_pt_astp;
      std::vector<TH1F*> m_h_b_pt_diff;
      std::vector<TH2F*> m_h_b_pt_2d;

      std::vector<TH1F*> m_h_b_eta_stop;
      std::vector<TH1F*> m_h_b_eta_astp;
      std::vector<TH1F*> m_h_b_eta_diff;
      std::vector<TH2F*> m_h_b_eta_2d;

      std::vector<TH1F*> m_h_b_phi_stop;
      std::vector<TH1F*> m_h_b_phi_astp;
      std::vector<TH1F*> m_h_b_phi_diff;
      std::vector<TH2F*> m_h_b_phi_2d;

      // differrence and 2d plots of b quark kinematics vs lepton kinematics
      // do for right and wrong pairing
      std::vector<TH1F*> m_h_right_pair_bl_dpt;
      std::vector<TH1F*> m_h_right_pair_bl_deta;
      std::vector<TH1F*> m_h_right_pair_bl_dphi;
      std::vector<TH1F*> m_h_right_pair_bl_dr;

      std::vector<TH2F*> m_h_right_pair_bl_pt_2d;
      std::vector<TH2F*> m_h_right_pair_bl_eta_2d;
      std::vector<TH2F*> m_h_right_pair_bl_phi_2d;

      std::vector<TH1F*> m_h_wrong_pair_bl_dpt;
      std::vector<TH1F*> m_h_wrong_pair_bl_deta;
      std::vector<TH1F*> m_h_wrong_pair_bl_dphi;
      std::vector<TH1F*> m_h_wrong_pair_bl_dr;

      std::vector<TH2F*> m_h_wrong_pair_bl_pt_2d;
      std::vector<TH2F*> m_h_wrong_pair_bl_eta_2d;
      std::vector<TH2F*> m_h_wrong_pair_bl_phi_2d;

      // m_bl invariant mass plots - do for right and wrong pairing
      std::vector<TH1F*> m_h_right_pair_mbl_all;
      std::vector<TH1F*> m_h_right_pair_mbl_stop;
      std::vector<TH1F*> m_h_right_pair_mbl_astp;
      std::vector<TH1F*> m_h_right_pair_mbl_diff;
      std::vector<TH2F*> m_h_right_pair_mbl_2d;

      std::vector<TH1F*> m_h_wrong_pair_mbl_all;
      std::vector<TH1F*> m_h_wrong_pair_mbl_0;
      std::vector<TH1F*> m_h_wrong_pair_mbl_1;
      std::vector<TH1F*> m_h_wrong_pair_mbl_diff;
      std::vector<TH2F*> m_h_wrong_pair_mbl_2d;
  };

  // =============================================================================
  // = mbl
  // =============================================================================
  class Mbl : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Mbl();

      virtual void FillSpecial( const TruthNtuple::FLAVOR_CHANNEL
                              , const std::vector<TruthNtuple::Electron*>&
                              , const std::vector<TruthNtuple::Muon*>&
                              , const std::vector<TruthNtuple::Particle*>&
                              );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_mbl_truth;
      std::vector<TH1F*> m_h_mbl_dphi_matching;
      std::vector<TH1F*> m_h_mbl_dr_matching;
  };
}

#endif
