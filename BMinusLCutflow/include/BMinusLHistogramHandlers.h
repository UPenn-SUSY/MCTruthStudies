#ifndef BMINUSL_HISTOGRAMHANDLERS_H
#define BMINUSL_HISTOGRAMHANDLERS_H

#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TruthNtupleLooper/include/TruthNtupleEnums.h"
#include "HistogramHandlers/include/HistogramHandlers.h"

// =============================================================================
class TFile;

namespace TruthNtuple
{
  class Particle;
  // class Lepton;
  // class Electron;
  // class Muon;
  // class Jet;
  class Met;
}

// =============================================================================
namespace HistogramHandlers
{
  // ===========================================================================
  // = Stop kinematics
  // ===========================================================================
  class StopKinematics : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      StopKinematics();

      virtual void FillSpecial( const TruthNtuple::FLAVOR_CHANNEL
                              , const std::vector<TruthNtuple::Particle*>&
                              );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1D*> m_h_e_com;
      std::vector<TH1D*> m_h_m;

      std::vector<TH1D*> m_h_e_all;
      std::vector<TH1D*> m_h_e_stop;
      std::vector<TH1D*> m_h_e_astp;
      std::vector<TH1D*> m_h_e_diff;
      std::vector<TH2D*> m_h_e_2d;

      std::vector<TH1D*> m_h_p_all;
      std::vector<TH1D*> m_h_p_stop;
      std::vector<TH1D*> m_h_p_astp;
      std::vector<TH1D*> m_h_p_diff;
      std::vector<TH2D*> m_h_p_2d;

      std::vector<TH1D*> m_h_pt_all;
      std::vector<TH1D*> m_h_pt_stop; // stop
      std::vector<TH1D*> m_h_pt_astp; // anti-stop
      std::vector<TH1D*> m_h_pt_diff;
      std::vector<TH2D*> m_h_pt_2d;

      std::vector<TH1D*> m_h_pz_all;
      std::vector<TH1D*> m_h_pz_stop; // stop
      std::vector<TH1D*> m_h_pz_astp; // anti-stop
      std::vector<TH1D*> m_h_pz_diff;
      std::vector<TH2D*> m_h_pz_2d;

      std::vector<TH1D*> m_h_p_over_m_all;
      std::vector<TH1D*> m_h_p_over_m_stop;
      std::vector<TH1D*> m_h_p_over_m_astp;
      std::vector<TH1D*> m_h_p_over_m_diff;
      std::vector<TH2D*> m_h_p_over_m_2d;

      std::vector<TH1D*> m_h_pt_over_m_all;
      std::vector<TH1D*> m_h_pt_over_m_stop; // stop
      std::vector<TH1D*> m_h_pt_over_m_astp; // anti-stop
      std::vector<TH1D*> m_h_pt_over_m_diff;
      std::vector<TH2D*> m_h_pt_over_m_2d;

      std::vector<TH1D*> m_h_eta_all;
      std::vector<TH1D*> m_h_eta_stop; // stop
      std::vector<TH1D*> m_h_eta_astp; // anti-stop
      std::vector<TH1D*> m_h_eta_diff;
      std::vector<TH2D*> m_h_eta_2d;

      std::vector<TH1D*> m_h_y_all;
      std::vector<TH1D*> m_h_y_stop; // stop
      std::vector<TH1D*> m_h_y_astp; // anti-stop
      std::vector<TH1D*> m_h_y_diff;
      std::vector<TH2D*> m_h_y_2d;

      std::vector<TH1D*> m_h_phi_all;
      std::vector<TH1D*> m_h_phi_stop; // stop
      std::vector<TH1D*> m_h_phi_astp; // anti-stop
      std::vector<TH1D*> m_h_phi_diff;
      std::vector<TH2D*> m_h_phi_2d;
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
      std::vector<TH1D*> m_h_pt_all;
      std::vector<TH1D*> m_h_pt_0;
      std::vector<TH1D*> m_h_pt_1;
      std::vector<TH1D*> m_h_pt_diff;
      std::vector<TH2D*> m_h_pt_2d;

      std::vector<TH1D*> m_h_eta_all;
      std::vector<TH1D*> m_h_eta_0;
      std::vector<TH1D*> m_h_eta_1;
      std::vector<TH1D*> m_h_eta_diff;
      std::vector<TH2D*> m_h_eta_2d;

      std::vector<TH1D*> m_h_phi_all;
      std::vector<TH1D*> m_h_phi_0;
      std::vector<TH1D*> m_h_phi_1;
      std::vector<TH1D*> m_h_phi_diff;
      std::vector<TH2D*> m_h_phi_2d;

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
                              , const std::vector<TruthNtuple::Particle*>& el
                              , const std::vector<TruthNtuple::Particle*>& mu
                              , const std::vector<TruthNtuple::Particle*>& b
                              , const std::vector<TruthNtuple::Particle*>& stop_list
                              );
      virtual void write(TFile*);

      bool sortObjects( const std::vector<TruthNtuple::Particle*>& el
                      , const std::vector<TruthNtuple::Particle*>& mu
                      , const std::vector<TruthNtuple::Particle*>& b
                      );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      // helper pointers to keep track of leptons from stop and anti-stop
      std::vector<TruthNtuple::Particle*> m_l_list;
      TruthNtuple::Particle* m_l_from_stop;
      TruthNtuple::Particle* m_l_from_astp;

      // helper pointers to keep track of b quarks from stop and anti-stop
      std::vector<TruthNtuple::Particle*> m_b_list;
      TruthNtuple::Particle* m_b_from_stop;
      TruthNtuple::Particle* m_b_from_astp;

      // 1d plots of lepton kinematics. separate plots for leptons from stop
      // and anti-stop
      std::vector<TH1D*> m_h_l_pt_all;
      std::vector<TH1D*> m_h_l_pt_stop;
      std::vector<TH1D*> m_h_l_pt_astp;
      std::vector<TH1D*> m_h_l_pt_diff;
      std::vector<TH2D*> m_h_l_pt_2d;

      std::vector<TH1D*> m_h_l_eta_all;
      std::vector<TH1D*> m_h_l_eta_stop;
      std::vector<TH1D*> m_h_l_eta_astp;
      std::vector<TH1D*> m_h_l_eta_diff;
      std::vector<TH2D*> m_h_l_eta_2d;

      std::vector<TH1D*> m_h_l_phi_all;
      std::vector<TH1D*> m_h_l_phi_stop;
      std::vector<TH1D*> m_h_l_phi_astp;
      std::vector<TH1D*> m_h_l_phi_diff;
      std::vector<TH2D*> m_h_l_phi_2d;

      // 1d plots of b quark kinematics. separate plots for leptons from stop
      // and anti-stop
      std::vector<TH1D*> m_h_b_pt_all;
      std::vector<TH1D*> m_h_b_pt_stop;
      std::vector<TH1D*> m_h_b_pt_astp;
      std::vector<TH1D*> m_h_b_pt_diff;
      std::vector<TH2D*> m_h_b_pt_2d;

      std::vector<TH1D*> m_h_b_eta_all;
      std::vector<TH1D*> m_h_b_eta_stop;
      std::vector<TH1D*> m_h_b_eta_astp;
      std::vector<TH1D*> m_h_b_eta_diff;
      std::vector<TH2D*> m_h_b_eta_2d;

      std::vector<TH1D*> m_h_b_phi_all;
      std::vector<TH1D*> m_h_b_phi_stop;
      std::vector<TH1D*> m_h_b_phi_astp;
      std::vector<TH1D*> m_h_b_phi_diff;
      std::vector<TH2D*> m_h_b_phi_2d;

      // differrence and 2d plots of b quark kinematics vs lepton kinematics
      // do for right and wrong pairing
      std::vector<TH1D*> m_h_right_pair_bl_dpt;
      std::vector<TH1D*> m_h_right_pair_bl_deta;
      std::vector<TH1D*> m_h_right_pair_bl_dphi;
      std::vector<TH1D*> m_h_right_pair_bl_dr;

      std::vector<TH2D*> m_h_right_pair_bl_pt_2d;
      std::vector<TH2D*> m_h_right_pair_bl_eta_2d;
      std::vector<TH2D*> m_h_right_pair_bl_phi_2d;

      std::vector<TH1D*> m_h_wrong_pair_bl_dpt;
      std::vector<TH1D*> m_h_wrong_pair_bl_deta;
      std::vector<TH1D*> m_h_wrong_pair_bl_dphi;
      std::vector<TH1D*> m_h_wrong_pair_bl_dr;

      std::vector<TH2D*> m_h_wrong_pair_bl_pt_2d;
      std::vector<TH2D*> m_h_wrong_pair_bl_eta_2d;
      std::vector<TH2D*> m_h_wrong_pair_bl_phi_2d;

      // m_bl invariant mass plots - do for right and wrong and chosen pairings
      std::vector<TH1D*> m_h_right_pair_mbl_all;
      std::vector<TH1D*> m_h_right_pair_mbl_stop;
      std::vector<TH1D*> m_h_right_pair_mbl_astp;
      std::vector<TH1D*> m_h_right_pair_mbl_diff;
      std::vector<TH1D*> m_h_right_pair_mbl_ratio;
      std::vector<TH1D*> m_h_right_pair_mbl_sq_sum;
      std::vector<TH2D*> m_h_right_pair_mbl_2d;

      std::vector<TH1D*> m_h_wrong_pair_mbl_all;
      std::vector<TH1D*> m_h_wrong_pair_mbl_0;
      std::vector<TH1D*> m_h_wrong_pair_mbl_1;
      std::vector<TH1D*> m_h_wrong_pair_mbl_diff;
      std::vector<TH1D*> m_h_wrong_pair_mbl_ratio;
      std::vector<TH1D*> m_h_wrong_pair_mbl_sq_sum;
      std::vector<TH2D*> m_h_wrong_pair_mbl_2d;

      std::vector<TH2D*> m_h_wrong_vs_right_mbl_diff;
      std::vector<TH2D*> m_h_wrong_vs_right_mbl_ratio;
      std::vector<TH2D*> m_h_wrong_vs_right_mbl_sq_sum;

      std::vector<TH1D*> m_h_diff_pair_mbl_all;
      std::vector<TH1D*> m_h_diff_pair_mbl_0;
      std::vector<TH1D*> m_h_diff_pair_mbl_1;
      std::vector<TH1D*> m_h_diff_pair_mbl_diff;
      std::vector<TH1D*> m_h_diff_pair_mbl_ratio;
      std::vector<TH1D*> m_h_diff_pair_mbl_sq_sum;
      std::vector<TH2D*> m_h_diff_pair_mbl_2d;
      std::vector<TH1D*> m_h_diff_pair_cor_pairing;

      std::vector<TH1D*> m_h_ratio_pair_mbl_all;
      std::vector<TH1D*> m_h_ratio_pair_mbl_0;
      std::vector<TH1D*> m_h_ratio_pair_mbl_1;
      std::vector<TH1D*> m_h_ratio_pair_mbl_diff;
      std::vector<TH1D*> m_h_ratio_pair_mbl_ratio;
      std::vector<TH1D*> m_h_ratio_pair_mbl_sq_sum;
      std::vector<TH2D*> m_h_ratio_pair_mbl_2d;
      std::vector<TH1D*> m_h_ratio_pair_cor_pairing;

      std::vector<TH1D*> m_h_sq_sum_pair_mbl_all;
      std::vector<TH1D*> m_h_sq_sum_pair_mbl_0;
      std::vector<TH1D*> m_h_sq_sum_pair_mbl_1;
      std::vector<TH1D*> m_h_sq_sum_pair_mbl_diff;
      std::vector<TH1D*> m_h_sq_sum_pair_mbl_ratio;
      std::vector<TH1D*> m_h_sq_sum_pair_mbl_sq_sum;
      std::vector<TH2D*> m_h_sq_sum_pair_mbl_2d;
      std::vector<TH1D*> m_h_sq_sum_pair_cor_pairing;

      // pt_bl invariant mass plots - do for right and wrong pairing
      std::vector<TH1D*> m_h_right_pair_ptbl_all;
      std::vector<TH1D*> m_h_right_pair_ptbl_stop;
      std::vector<TH1D*> m_h_right_pair_ptbl_astp;
      std::vector<TH1D*> m_h_right_pair_ptbl_diff;
      std::vector<TH1D*> m_h_right_pair_ptbl_ratio;
      std::vector<TH1D*> m_h_right_pair_ptbl_sq_sum;
      std::vector<TH2D*> m_h_right_pair_ptbl_2d;

      std::vector<TH1D*> m_h_wrong_pair_ptbl_all;
      std::vector<TH1D*> m_h_wrong_pair_ptbl_0;
      std::vector<TH1D*> m_h_wrong_pair_ptbl_1;
      std::vector<TH1D*> m_h_wrong_pair_ptbl_diff;
      std::vector<TH1D*> m_h_wrong_pair_ptbl_ratio;
      std::vector<TH1D*> m_h_wrong_pair_ptbl_sq_sum;
      std::vector<TH2D*> m_h_wrong_pair_ptbl_2d;

      std::vector<TH2D*> m_h_wrong_vs_right_ptbl_diff;
      std::vector<TH2D*> m_h_wrong_vs_right_ptbl_ratio;
      std::vector<TH2D*> m_h_wrong_vs_right_ptbl_sq_sum;

      // pt of stop vs bl pair
      std::vector<TH1D*> m_h_stop_vs_bl_px_all;
      std::vector<TH1D*> m_h_stop_vs_bl_px_stop;
      std::vector<TH1D*> m_h_stop_vs_bl_px_astp;
      std::vector<TH2D*> m_h_stop_vs_bl_px_2d;
      std::vector<TH1D*> m_h_stop_vs_bl_px_event;

      std::vector<TH1D*> m_h_stop_vs_bl_py_all;
      std::vector<TH1D*> m_h_stop_vs_bl_py_stop;
      std::vector<TH1D*> m_h_stop_vs_bl_py_astp;
      std::vector<TH2D*> m_h_stop_vs_bl_py_2d;
      std::vector<TH1D*> m_h_stop_vs_bl_py_event;

      std::vector<TH1D*> m_h_stop_vs_bl_pz_all;
      std::vector<TH1D*> m_h_stop_vs_bl_pz_stop;
      std::vector<TH1D*> m_h_stop_vs_bl_pz_astp;
      std::vector<TH2D*> m_h_stop_vs_bl_pz_2d;
      std::vector<TH1D*> m_h_stop_vs_bl_pz_event;

      std::vector<TH1D*> m_h_stop_vs_bl_pt_all;
      std::vector<TH1D*> m_h_stop_vs_bl_pt_stop;
      std::vector<TH1D*> m_h_stop_vs_bl_pt_astp;
      std::vector<TH2D*> m_h_stop_vs_bl_pt_2d;
      std::vector<TH1D*> m_h_stop_vs_bl_pt_event;
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
                              , const std::vector<TruthNtuple::Particle*>& el
                              , const std::vector<TruthNtuple::Particle*>& mu
                              , const std::vector<TruthNtuple::Particle*>& b_jets
                              );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1D*> m_h_ratio_pair_mbl_all;
      std::vector<TH1D*> m_h_ratio_pair_mbl_0;
      std::vector<TH1D*> m_h_ratio_pair_mbl_1;
      std::vector<TH1D*> m_h_ratio_pair_mbl_diff;
      std::vector<TH1D*> m_h_ratio_pair_mbl_ratio;
      std::vector<TH1D*> m_h_ratio_pair_mbl_sq_sum;
      std::vector<TH2D*> m_h_ratio_pair_mbl_2d;
      // std::vector<TH1D*> m_h_ratio_pair_cor_pairing;

  };
}

#endif
