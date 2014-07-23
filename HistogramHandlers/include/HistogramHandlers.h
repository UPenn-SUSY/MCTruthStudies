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
  class Particle;
  // class Electron;
  // class Muon;
  // class Jet;
  class Met;
}

// =============================================================================
namespace HistogramHandlers
{
  // =============================================================================
  // = handle - used as parent class for other histograms
  // =============================================================================
  class Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Handle();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Particle*>& el
                       , const std::vector<TruthNtuple::Particle*>& mu
                       , const std::vector<TruthNtuple::Particle*>& jet
                       , const std::vector<TruthNtuple::Particle*>& quark
                       , const TruthNtuple::Met&
		       , double m_event_weight
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
  };

  // =============================================================================
  // = flavor channel
  // =============================================================================
  class FlavorChannel : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      FlavorChannel();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Particle*>& el
                       , const std::vector<TruthNtuple::Particle*>& mu
                       , const std::vector<TruthNtuple::Particle*>& jet
                       , const std::vector<TruthNtuple::Particle*>& quark
                       , const TruthNtuple::Met&
		       , double m_event_weight
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_flavor_channel;
  };

  // =============================================================================
  // = Object multiplicity
  // =============================================================================
  class ObjectMultiplicity : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      ObjectMultiplicity();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Particle*>& el
                       , const std::vector<TruthNtuple::Particle*>& mu
                       , const std::vector<TruthNtuple::Particle*>& jet
                       , const std::vector<TruthNtuple::Particle*>& quark
                       , const TruthNtuple::Met&
		       , double m_event_weight
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_num_lep;
      std::vector<TH1F*> m_h_num_jet;
  };

  // =============================================================================
  // = Lepton kinematics
  // =============================================================================
  class LeptonKinematics : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      LeptonKinematics();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Particle*>& el
                       , const std::vector<TruthNtuple::Particle*>& mu
                       , const std::vector<TruthNtuple::Particle*>& jet
                       , const std::vector<TruthNtuple::Particle*>& quark
                       , const TruthNtuple::Met&
		       , double m_event_weight
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
  // = JetKinematics
  // =============================================================================
  class JetKinematics : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      JetKinematics();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Particle*>& el
                       , const std::vector<TruthNtuple::Particle*>& mu
                       , const std::vector<TruthNtuple::Particle*>& jet
                       , const std::vector<TruthNtuple::Particle*>& quark
                       , const TruthNtuple::Met&
		       , double m_event_weight
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
  // = met
  // =============================================================================
  class Met : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Met();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Particle*>& el
                       , const std::vector<TruthNtuple::Particle*>& mu
                       , const std::vector<TruthNtuple::Particle*>& jet
                       , const std::vector<TruthNtuple::Particle*>& quark
                       , const TruthNtuple::Met&
		       , double m_event_weight
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_met;
      std::vector<TH1F*> m_h_metrel;
      std::vector<TH1F*> m_h_met_sig;
  };

  // =============================================================================
  // = mll
  // =============================================================================
  class Mll : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Mll();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Particle*>& el
                       , const std::vector<TruthNtuple::Particle*>& mu
                       , const std::vector<TruthNtuple::Particle*>& jet
                       , const std::vector<TruthNtuple::Particle*>& quark
                       , const TruthNtuple::Met&
		       , double m_event_weight
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_mll;
  };

  // =============================================================================
  // = mjl
  // =============================================================================
  class Mjl : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  public:
    Mjl();

      virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Particle*>& el
                       , const std::vector<TruthNtuple::Particle*>& mu
                       , const std::vector<TruthNtuple::Particle*>& jet
                       , const std::vector<TruthNtuple::Particle*>& quark
                       , const TruthNtuple::Met&
		       , double m_event_weight
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_mjl_truth;
      std::vector<TH1F*> m_h_mjl_dphi_matching;
      std::vector<TH1F*> m_h_mjl_dr_matching;
  };

  // =============================================================================
  // = Ht
  // =============================================================================
  class Ht : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  public:
    Ht();
    
    virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Particle*>& el
                       , const std::vector<TruthNtuple::Particle*>& mu
                       , const std::vector<TruthNtuple::Particle*>& jet
                       , const std::vector<TruthNtuple::Particle*>& quark
                       , const TruthNtuple::Met&
		       , double m_event_weight
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_ht;
  };

  // =============================================================================
  // = Dr
  // =============================================================================
  class Dr : public Handle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  public:
    Dr();
    
    virtual void Fill( const TruthNtuple::FLAVOR_CHANNEL
                       , const std::vector<TruthNtuple::Particle*>& el
                       , const std::vector<TruthNtuple::Particle*>& mu
                       , const std::vector<TruthNtuple::Particle*>& jet
                       , const std::vector<TruthNtuple::Particle*>& quark
                       , const TruthNtuple::Met&
		       , double m_event_weight
                       );
      virtual void write(TFile*);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      std::vector<TH1F*> m_h_dr_ll;
      std::vector<TH1F*> m_h_dr_qq;
      std::vector<TH1F*> m_h_dr_l0q0;
      std::vector<TH1F*> m_h_dr_l0q1;
      std::vector<TH1F*> m_h_dr_l1q0;
      std::vector<TH1F*> m_h_dr_l1q1;
/*       std::vector<TH2F*> m_h_dr_lq0vlq1; */
/*       std::vector<TH1F*> m_h_dr_lsameq0; */
/*       std::vector<TH1F*> m_h_dr_lsameq1; */
  };
}
#endif
