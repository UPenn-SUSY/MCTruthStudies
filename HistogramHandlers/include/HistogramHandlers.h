#ifndef HISTOGRAMHANDLERS_H
#define HISTOGRAMHANDLERS_H

#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TruthNtupleLooper/include/TruthNtupleEnums.h"

// =============================================================================
class TFile;

namespace TruthNtuple
{
  class Electron;
  class Muon;
  class Jet;
  class Met;
}

// =============================================================================
namespace HistogramHandlers
{
  // =============================================================================
  class Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Handle();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                      , const std::vector<TruthNtuple::Electron*>&
                      , const std::vector<TruthNtuple::Muon*>&
                      , const std::vector<TruthNtuple::Jet*>&
                      , const TruthNtuple::Met&
                      );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
  };

  // =============================================================================
  class ObjectMultiplicity : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      ObjectMultiplicity();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Electron*>&
                       , const std::vector<TruthNtuple::Muon*>&
                       , const std::vector<TruthNtuple::Jet*>&
                       , const TruthNtuple::Met&
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_num_lep;
      std::vector<TH1F*> m_h_num_jet;
  };

  // =============================================================================
  class LeptonPt : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      LeptonPt();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Electron*>&
                       , const std::vector<TruthNtuple::Muon*>&
                       , const std::vector<TruthNtuple::Jet*>&
                       , const TruthNtuple::Met&
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_pt_all;
      std::vector<TH1F*> m_h_pt_0;
      std::vector<TH1F*> m_h_pt_1;
      std::vector<TH1F*> m_h_pt_diff;
      std::vector<TH2F*> m_h_pt_2d;
  };

  // =============================================================================
  class LeptonEta : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      LeptonEta();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Electron*>&
                       , const std::vector<TruthNtuple::Muon*>&
                       , const std::vector<TruthNtuple::Jet*>&
                       , const TruthNtuple::Met&
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_eta_all;
      std::vector<TH1F*> m_h_eta_0;
      std::vector<TH1F*> m_h_eta_1;
      std::vector<TH1F*> m_h_eta_diff;
      std::vector<TH2F*> m_h_eta_2d;
  };

  // =============================================================================
  class LeptonPhi : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      LeptonPhi();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Electron*>&
                       , const std::vector<TruthNtuple::Muon*>&
                       , const std::vector<TruthNtuple::Jet*>&
                       , const TruthNtuple::Met&
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_phi_all;
      std::vector<TH1F*> m_h_phi_0;
      std::vector<TH1F*> m_h_phi_1;
      std::vector<TH1F*> m_h_phi_diff;
      std::vector<TH2F*> m_h_phi_2d;
  };

  // =============================================================================
  class JetPt : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      JetPt();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Electron*>&
                       , const std::vector<TruthNtuple::Muon*>&
                       , const std::vector<TruthNtuple::Jet*>&
                       , const TruthNtuple::Met&
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_pt_all;
      std::vector<TH1F*> m_h_pt_0;
      std::vector<TH1F*> m_h_pt_1;
      std::vector<TH1F*> m_h_pt_diff;
      std::vector<TH2F*> m_h_pt_2d;
  };

  // =============================================================================
  class JetEta : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      JetEta();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Electron*>&
                       , const std::vector<TruthNtuple::Muon*>&
                       , const std::vector<TruthNtuple::Jet*>&
                       , const TruthNtuple::Met&
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_eta_all;
      std::vector<TH1F*> m_h_eta_0;
      std::vector<TH1F*> m_h_eta_1;
      std::vector<TH1F*> m_h_eta_diff;
      std::vector<TH2F*> m_h_eta_2d;
  };

  // =============================================================================
  class JetPhi : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      JetPhi();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Electron*>&
                       , const std::vector<TruthNtuple::Muon*>&
                       , const std::vector<TruthNtuple::Jet*>&
                       , const TruthNtuple::Met&
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_phi_all;
      std::vector<TH1F*> m_h_phi_0;
      std::vector<TH1F*> m_h_phi_1;
      std::vector<TH1F*> m_h_phi_diff;
      std::vector<TH2F*> m_h_phi_2d;
  };

  // =============================================================================
  class Met : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Met();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Electron*>&
                       , const std::vector<TruthNtuple::Muon*>&
                       , const std::vector<TruthNtuple::Jet*>&
                       , const TruthNtuple::Met&
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_met;
      std::vector<TH1F*> m_h_metrel;
  };

  // =============================================================================
  class Mll : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Mll();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Electron*>&
                       , const std::vector<TruthNtuple::Muon*>&
                       , const std::vector<TruthNtuple::Jet*>&
                       , const TruthNtuple::Met&
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_mll;
  };

  // =============================================================================
  class Mbl : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Mbl();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Electron*>&
                       , const std::vector<TruthNtuple::Muon*>&
                       , const std::vector<TruthNtuple::Jet*>&
                       , const TruthNtuple::Met&
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
