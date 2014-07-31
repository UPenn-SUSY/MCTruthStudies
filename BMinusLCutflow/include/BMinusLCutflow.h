#ifndef BMINUSLCUTFLOW_H
#define BMINUSLCUTFLOW_H

#include "TruthNtupleLooper/include/TruthNtupleLooper.h"
#include "TruthNtupleLooper/include/TruthNtupleEnums.h"
// #include "include/PDFTool.h"
// =============================================================================
class TTree;
class TH1D;
class TH2D;
class TRandom;

namespace TruthNtuple
{
  class Electron;
  class Muon;
  class Jet;
  // class Met;
}

namespace HistogramHandlers
{
  class Handle;
  class Mbl;
  class BLPairKinematics;
  class QuarkKinematics;
  class StopKinematics;
}

namespace BMinusL
{
  class BMinusLStandAloneHistograms;
}

// =============================================================================
namespace BMinusL
{
  class Cutflow : public TruthNtuple::TruthNtupleLooper
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Cutflow(TTree *tree=0, bool isSignal = true);
      ~Cutflow();

      virtual void clearObjects();
      virtual void processEvent();

      void writeToFile();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      bool m_is_signal;
      void doObjectSelection();
      void print();

      bool isLeptonFromTauFromStop(const TruthNtuple::Particle*);
      double scalePDF(double m_event_weight);
      void broadenResolution();
      double findQuarkSigma(double x);
      double findElSigma(double x);
      double findMuSigma(double x);

      TruthNtuple::FLAVOR_CHANNEL m_flavor_channel;

      // keep track of stops and ALL b quarks
      std::vector<TruthNtuple::Particle*> m_truth_stops;
      std::vector<TruthNtuple::Particle*> m_truth_electrons;
      std::vector<TruthNtuple::Particle*> m_truth_muons;
      std::vector<TruthNtuple::Particle*> m_truth_taus;
      std::vector<TruthNtuple::Particle*> m_truth_b_quarks;
      std::vector<TruthNtuple::Jet*> m_b_jets;

      // objects matched to SUSY mother
      std::vector<TruthNtuple::Particle*> m_daughter_el;
      std::vector<TruthNtuple::Particle*> m_daughter_mu;
      std::vector<TruthNtuple::Particle*> m_daughter_tau;
      std::vector<TruthNtuple::Particle*> m_daughter_b_quarks;

      std::vector<TruthNtuple::Particle*> m_leading_b_jets;

      TruthNtuple::Met m_met;

      std::vector<HistogramHandlers::Handle*> m_histograms;
      HistogramHandlers::Mbl*                 m_h_mbl;
      HistogramHandlers::BLPairKinematics*    m_h_bl_pair_kinematics;
      HistogramHandlers::QuarkKinematics*     m_h_quark_kinematics;
      HistogramHandlers::StopKinematics*      m_h_stop_kinematics;

      BMinusL::BMinusLStandAloneHistograms* m_sa_hists;

      friend class BMinusLStandAloneHistograms;
  };

  // ===========================================================================
  // = Stand alone histograms
  // ===========================================================================
  class BMinusLStandAloneHistograms
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      BMinusLStandAloneHistograms();

      void Fill( const BMinusL::Cutflow* );

      void deltaRCalc(const TruthNtuple::Particle*
		      ,const TruthNtuple::Particle*
		      ,const TruthNtuple::Particle*
		      ,const TruthNtuple::Particle*
		      ,double&
		      ,double&
		      ,double&
		      ,double&
		      );

      void write(TFile*);

      TH2D* calcEffPt(TH2D*, std::string tag="");
      TH2D* calcEffFid(TH2D*, TH2D*);
      TH1D* m_h_meff;
      TH2D* m_h_fc_all__pt_event_b1vl1;
      TH2D* m_h_fc_all__pt_event_b1ve1;
      TH2D* m_h_fc_all__pt_event_b1vm1;
      TH2D* m_h_fc_all__pt_event_b1vl1_eff;
      TH2D* m_h_fc_all__pt_event_b1ve1_eff;
      TH2D* m_h_fc_all__pt_event_b1vm1_eff;
      TH1D* m_h_fc_all__eta_event_eff;
      TH1D* m_h_fc_ee__eta_event_eff;
      TH1D* m_h_fc_em__eta_event_eff;
      TH1D* m_h_fc_me__eta_event_eff;
      TH1D* m_h_fc_mm__eta_event_eff;
      TH1D* m_h_fc_all__lep_eta_event_eff;
      TH1D* m_h_fc_all__quark_eta_event_eff;
      TH1D* m_h_fc_all__lep_fiducial_eventall_pass;
      TH1D* m_h_fc_all__lep_fiducial_eventall_fail;
      TH1D* m_h_fc_all__lep_fiducial_eventlep_pass;
      TH1D* m_h_fc_all__lep_fiducial_eventlep_fail;
      TH1D* m_h_fc_all__lep_fiducial_single_pass;
      TH1D* m_h_fc_all__lep_fiducial_single_fail;
      TH1D* m_h_fc_all__quark_fiducial_eventall_pass;
      TH1D* m_h_fc_all__quark_fiducial_eventall_fail;
      TH1D* m_h_fc_all__quark_fiducial_eventquark_pass;
      TH1D* m_h_fc_all__quark_fiducial_eventquark_fail;
      TH1D* m_h_fc_all__quark_fiducial_single_pass;
      TH1D* m_h_fc_all__quark_fiducial_single_fail;
      //TH1D* m_h_fc_all__fiducial_event_e_pass;
      //TH1D* m_h_fc_all__fiducial_event_e_fail;
      //TH1D* m_h_fc_all__fiducial_event_m_pass;
      //TH1D* m_h_fc_all__fiducial_event_m_fail;
      TH2D* m_h_fc_all__fiducial_event_b1vl1_pass;
      TH2D* m_h_fc_all__fiducial_event_b1vl1_fail;
      TH2D* m_h_fc_all__fiducial_event_b1vl1_eff;
      TH1D* m_h_fc_all__lep_deltaRq;
      TH2D* m_h_fc_all__lep_deltaRq_l0vl1;
      TH1D* m_h_fc_all__lep_deltaRsameq0;
      TH1D* m_h_fc_all__lep_deltaRsameq1;
      TH1D* m_h_fc_all__lep_deltaRl;
      TH1D* m_h_fc_all__quark_deltaRq;
  };

}

#endif
