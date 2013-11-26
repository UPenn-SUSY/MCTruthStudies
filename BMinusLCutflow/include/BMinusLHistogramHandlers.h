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
