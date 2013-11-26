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

      // // objects passing selection criteria
      // std::vector<TruthNtuple::Electron*> m_selected_el;
      // std::vector<TruthNtuple::Muon*>     m_selected_mu;
      // std::vector<TruthNtuple::Jet*>      m_selected_jet;

      // objects matched to SUSY mother
      std::vector<TruthNtuple::Electron*> m_daughter_el;
      std::vector<TruthNtuple::Muon*>     m_daughter_mu;
      std::vector<TruthNtuple::Jet*>      m_daughter_jet;

      TruthNtuple::Met m_met;

      std::vector<HistogramHandlers::Handle*> m_histograms;
  };
}

#endif
