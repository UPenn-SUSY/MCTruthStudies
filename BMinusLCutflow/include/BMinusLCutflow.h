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

// =============================================================================
namespace BMinusL
{
  class Cutflow : public TruthNtuple::TruthNtupleLooper
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Cutflow(TTree *tree=0);
      ~Cutflow();

      virtual void clearObjects();
      virtual void processEvent();

      void writeToFile();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      void doObjectSelection();

      TruthNtuple::FLAVOR_CHANNEL m_flavor_channel;

      // keep track of stops and ALL b quarks
      std::vector<TruthNtuple::Particle*> m_stops;
      std::vector<TruthNtuple::Particle*> m_b_quarks;

      // objects matched to SUSY mother
      std::vector<TruthNtuple::Electron*> m_daughter_el;
      std::vector<TruthNtuple::Muon*>     m_daughter_mu;
      std::vector<TruthNtuple::Particle*> m_daughter_b_quarks;
      std::vector<TruthNtuple::Jet*>      m_daughter_jet;

      TruthNtuple::Met m_met;

      std::vector<HistogramHandlers::Handle*> m_histograms;
      HistogramHandlers::Mbl*    m_h_mbl;
      HistogramHandlers::BLPairKinematics* m_h_bl_pair_kinematics;
      HistogramHandlers::QuarkKinematics* m_h_quark_kinematics;
      HistogramHandlers::StopKinematics* m_h_stop_kinematics;
  };
}

#endif
