#ifndef BMINUSLCUTFLOW_H
#define BMINUSLCUTFLOW_H

#include "TruthNtupleLooper/include/TruthNtupleLooper.h"
#include "TruthNtupleLooper/include/TruthNtupleEnums.h"

// =============================================================================
class TTree;
class TH1D;
class TH2D;

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

      void write(TFile*);

      TH2D* calcEff(TH2D*, std::string tag="");

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      TH1D* m_h_meff;
      TH2D* m_h_pt_b1vsl1;
      TH2D* m_h_pt_b1vse1;
      TH2D* m_h_pt_b1vsm1;
      TH2D* m_h_pt_b1vsl1_eff;
      TH2D* m_h_pt_b1vse1_eff;
      TH2D* m_h_pt_b1vsm1_eff;
  };

}

#endif
