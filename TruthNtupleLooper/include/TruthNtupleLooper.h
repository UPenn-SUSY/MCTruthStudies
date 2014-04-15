#ifndef TRUTHNTUPLELOOPER_H
#define TRUTHNTUPLELOOPER_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

#include "TruthNtupleLooper/include/ObjectDefs.h"

// =============================================================================
class TBranch;

// =============================================================================
namespace TruthNtuple
{
  // ===========================================================================
  class TruthNtupleLooper {
    public :
      TruthNtupleLooper(TTree *tree=0);
      virtual ~TruthNtupleLooper();
      virtual Int_t    GetEntry(Long64_t entry);
      virtual Long64_t LoadTree(Long64_t entry);
      virtual void     Init(TTree *tree);
      virtual void     Loop();
      virtual Bool_t   Notify();

      virtual void clearObjects();
      virtual void constructObjects();
      virtual void processEvent();


      // Clean multiple instances of same particle in std::vector
      // protect us from times when particle "decays" to itself many times
      virtual void cleanParticleList( std::vector<TruthNtuple::Particle*>& );

    protected:
      std::vector<TruthNtuple::Particle> m_particle_list;
      std::vector<TruthNtuple::Electron> m_el_list;
      std::vector<TruthNtuple::Muon>     m_mu_list;
      std::vector<TruthNtuple::Jet>      m_jet_list;
      TruthNtuple::Met                   m_met;

      TTree          *fChain;   //!pointer to the analyzed TTree or TChain
      Int_t           fCurrent; //!current Tree number in a TChain

      // Declaration of leaf types
      UInt_t          RunNumber;
      UInt_t          EventNumber;
      UInt_t          mc_channel_number;
      UInt_t          mc_event_number;
      Float_t         mc_event_weight;
      Int_t           mcevt_n;
      std::vector<int>     *mcevt_signal_process_id;
      std::vector<int>     *mcevt_event_number;
      std::vector<double>  *mcevt_event_scale;
      std::vector<double>  *mcevt_alphaQCD;
      std::vector<double>  *mcevt_alphaQED;
      std::vector<int>     *mcevt_pdf_id1;
      std::vector<int>     *mcevt_pdf_id2;
      std::vector<double>  *mcevt_pdf_x1;
      std::vector<double>  *mcevt_pdf_x2;
      std::vector<double>  *mcevt_pdf_scale;
      std::vector<double>  *mcevt_pdf1;
      std::vector<double>  *mcevt_pdf2;
      std::vector<std::vector<double> > *mcevt_weight;
      Int_t           jet_AntiKt4TruthJets_n;
      std::vector<float>   *jet_AntiKt4TruthJets_E;
      std::vector<float>   *jet_AntiKt4TruthJets_pt;
      std::vector<float>   *jet_AntiKt4TruthJets_m;
      std::vector<float>   *jet_AntiKt4TruthJets_eta;
      std::vector<float>   *jet_AntiKt4TruthJets_phi;
      std::vector<float>   *jet_AntiKt4TruthJets_flavor_partonDR;
      std::vector<int>     *jet_AntiKt4TruthJets_flavor_partonFlavor;
      std::vector<int>     *jet_AntiKt4TruthJets_flavor_hadronFlavor;
      std::vector<int>     *jet_AntiKt4TruthJets_flavor_hadronPDGID;
      Int_t           jet_AntiKt4TruthJets_WZ_n;
      std::vector<float>   *jet_AntiKt4TruthJets_WZ_E;
      std::vector<float>   *jet_AntiKt4TruthJets_WZ_pt;
      std::vector<float>   *jet_AntiKt4TruthJets_WZ_m;
      std::vector<float>   *jet_AntiKt4TruthJets_WZ_eta;
      std::vector<float>   *jet_AntiKt4TruthJets_WZ_phi;
      std::vector<float>   *jet_AntiKt4TruthJets_WZ_flavor_partonDR;
      std::vector<int>     *jet_AntiKt4TruthJets_WZ_flavor_partonFlavor;
      std::vector<int>     *jet_AntiKt4TruthJets_WZ_flavor_hadronFlavor;
      std::vector<int>     *jet_AntiKt4TruthJets_WZ_flavor_hadronPDGID;
      Int_t           jet_AntiKt6TruthJets_WZ_n;
      std::vector<float>   *jet_AntiKt6TruthJets_WZ_E;
      std::vector<float>   *jet_AntiKt6TruthJets_WZ_pt;
      std::vector<float>   *jet_AntiKt6TruthJets_WZ_m;
      std::vector<float>   *jet_AntiKt6TruthJets_WZ_eta;
      std::vector<float>   *jet_AntiKt6TruthJets_WZ_phi;
      std::vector<float>   *jet_AntiKt6TruthJets_WZ_flavor_partonDR;
      std::vector<int>     *jet_AntiKt6TruthJets_WZ_flavor_partonFlavor;
      std::vector<int>     *jet_AntiKt6TruthJets_WZ_flavor_hadronFlavor;
      std::vector<int>     *jet_AntiKt6TruthJets_WZ_flavor_hadronPDGID;
      Int_t           mc_n;
      std::vector<float>   *mc_pt;
      std::vector<float>   *mc_m;
      std::vector<float>   *mc_eta;
      std::vector<float>   *mc_phi;
      std::vector<int>     *mc_status;
      std::vector<int>     *mc_barcode;
      std::vector<int>     *mc_pdgId;
      std::vector<float>   *mc_charge;
      std::vector<std::vector<int> > *mc_parents;
      std::vector<std::vector<int> > *mc_children;
      std::vector<float>   *mc_vx_x;
      std::vector<float>   *mc_vx_y;
      std::vector<float>   *mc_vx_z;
      std::vector<int>     *mc_vx_barcode;
      std::vector<std::vector<int> > *mc_child_index;
      std::vector<std::vector<int> > *mc_parent_index;
      Float_t         MET_Truth_NonInt_etx;
      Float_t         MET_Truth_NonInt_ety;
      Float_t         MET_Truth_Int_etx;
      Float_t         MET_Truth_Int_ety;
      Float_t         MET_Truth_IntCentral_etx;
      Float_t         MET_Truth_IntCentral_ety;
      Float_t         MET_Truth_IntFwd_etx;
      Float_t         MET_Truth_IntFwd_ety;
      Float_t         MET_Truth_IntOutCover_etx;
      Float_t         MET_Truth_IntOutCover_ety;
      Float_t         MET_Truth_IntMuons_etx;
      Float_t         MET_Truth_IntMuons_ety;
      Int_t           el_n;
      std::vector<float>   *el_pt;
      std::vector<float>   *el_m;
      std::vector<float>   *el_eta;
      std::vector<float>   *el_phi;
      std::vector<int>     *el_status;
      std::vector<int>     *el_barcode;
      std::vector<float>   *el_charge;
      std::vector<std::vector<int> > *el_parent_index;
      Int_t           mu_n;
      std::vector<float>   *mu_pt;
      std::vector<float>   *mu_m;
      std::vector<float>   *mu_eta;
      std::vector<float>   *mu_phi;
      std::vector<int>     *mu_status;
      std::vector<int>     *mu_barcode;
      std::vector<float>   *mu_charge;
      std::vector<std::vector<int> > *mu_parent_index;
      Int_t           tau_n;
      std::vector<float>   *tau_pt;
      std::vector<float>   *tau_m;
      std::vector<float>   *tau_eta;
      std::vector<float>   *tau_phi;
      std::vector<int>     *tau_status;
      std::vector<int>     *tau_barcode;
      std::vector<float>   *tau_charge;
      std::vector<std::vector<int> > *tau_parent_index;
      std::vector<std::vector<int> > *tau_decay_index;
      Int_t           ph_n;
      std::vector<float>   *ph_pt;
      std::vector<float>   *ph_m;
      std::vector<float>   *ph_eta;
      std::vector<float>   *ph_phi;
      std::vector<int>     *ph_status;
      std::vector<int>     *ph_barcode;

      // List of branches
      TBranch        *b_RunNumber;   //!
      TBranch        *b_EventNumber;   //!
      TBranch        *b_mc_channel_number;   //!
      TBranch        *b_mc_event_number;   //!
      TBranch        *b_mc_event_weight;   //!
      TBranch        *b_mcevt_n;   //!
      TBranch        *b_mcevt_signal_process_id;   //!
      TBranch        *b_mcevt_event_number;   //!
      TBranch        *b_mcevt_event_scale;   //!
      TBranch        *b_mcevt_alphaQCD;   //!
      TBranch        *b_mcevt_alphaQED;   //!
      TBranch        *b_mcevt_pdf_id1;   //!
      TBranch        *b_mcevt_pdf_id2;   //!
      TBranch        *b_mcevt_pdf_x1;   //!
      TBranch        *b_mcevt_pdf_x2;   //!
      TBranch        *b_mcevt_pdf_scale;   //!
      TBranch        *b_mcevt_pdf1;   //!
      TBranch        *b_mcevt_pdf2;   //!
      TBranch        *b_mcevt_weight;   //!
      TBranch        *b_jet_AntiKt4TruthJets_n;   //!
      TBranch        *b_jet_AntiKt4TruthJets_E;   //!
      TBranch        *b_jet_AntiKt4TruthJets_pt;   //!
      TBranch        *b_jet_AntiKt4TruthJets_m;   //!
      TBranch        *b_jet_AntiKt4TruthJets_eta;   //!
      TBranch        *b_jet_AntiKt4TruthJets_phi;   //!
      TBranch        *b_jet_AntiKt4TruthJets_flavor_partonDR;   //!
      TBranch        *b_jet_AntiKt4TruthJets_flavor_partonFlavor;   //!
      TBranch        *b_jet_AntiKt4TruthJets_flavor_hadronFlavor;   //!
      TBranch        *b_jet_AntiKt4TruthJets_flavor_hadronPDGID;   //!
      TBranch        *b_jet_AntiKt4TruthJets_WZ_n;   //!
      TBranch        *b_jet_AntiKt4TruthJets_WZ_E;   //!
      TBranch        *b_jet_AntiKt4TruthJets_WZ_pt;   //!
      TBranch        *b_jet_AntiKt4TruthJets_WZ_m;   //!
      TBranch        *b_jet_AntiKt4TruthJets_WZ_eta;   //!
      TBranch        *b_jet_AntiKt4TruthJets_WZ_phi;   //!
      TBranch        *b_jet_AntiKt4TruthJets_WZ_flavor_partonDR;   //!
      TBranch        *b_jet_AntiKt4TruthJets_WZ_flavor_partonFlavor;   //!
      TBranch        *b_jet_AntiKt4TruthJets_WZ_flavor_hadronFlavor;   //!
      TBranch        *b_jet_AntiKt4TruthJets_WZ_flavor_hadronPDGID;   //!
      TBranch        *b_jet_AntiKt6TruthJets_WZ_n;   //!
      TBranch        *b_jet_AntiKt6TruthJets_WZ_E;   //!
      TBranch        *b_jet_AntiKt6TruthJets_WZ_pt;   //!
      TBranch        *b_jet_AntiKt6TruthJets_WZ_m;   //!
      TBranch        *b_jet_AntiKt6TruthJets_WZ_eta;   //!
      TBranch        *b_jet_AntiKt6TruthJets_WZ_phi;   //!
      TBranch        *b_jet_AntiKt6TruthJets_WZ_flavor_partonDR;   //!
      TBranch        *b_jet_AntiKt6TruthJets_WZ_flavor_partonFlavor;   //!
      TBranch        *b_jet_AntiKt6TruthJets_WZ_flavor_hadronFlavor;   //!
      TBranch        *b_jet_AntiKt6TruthJets_WZ_flavor_hadronPDGID;   //!
      TBranch        *b_mc_n;   //!
      TBranch        *b_mc_pt;   //!
      TBranch        *b_mc_m;   //!
      TBranch        *b_mc_eta;   //!
      TBranch        *b_mc_phi;   //!
      TBranch        *b_mc_status;   //!
      TBranch        *b_mc_barcode;   //!
      TBranch        *b_mc_pdgId;   //!
      TBranch        *b_mc_charge;   //!
      TBranch        *b_mc_parents;   //!
      TBranch        *b_mc_children;   //!
      TBranch        *b_mc_vx_x;   //!
      TBranch        *b_mc_vx_y;   //!
      TBranch        *b_mc_vx_z;   //!
      TBranch        *b_mc_vx_barcode;   //!
      TBranch        *b_mc_child_index;   //!
      TBranch        *b_mc_parent_index;   //!
      TBranch        *b_MET_Truth_NonInt_etx;   //!
      TBranch        *b_MET_Truth_NonInt_ety;   //!
      TBranch        *b_MET_Truth_Int_etx;   //!
      TBranch        *b_MET_Truth_Int_ety;   //!
      TBranch        *b_MET_Truth_IntCentral_etx;   //!
      TBranch        *b_MET_Truth_IntCentral_ety;   //!
      TBranch        *b_MET_Truth_IntFwd_etx;   //!
      TBranch        *b_MET_Truth_IntFwd_ety;   //!
      TBranch        *b_MET_Truth_IntOutCover_etx;   //!
      TBranch        *b_MET_Truth_IntOutCover_ety;   //!
      TBranch        *b_MET_Truth_IntMuons_etx;   //!
      TBranch        *b_MET_Truth_IntMuons_ety;   //!
      TBranch        *b_el_n;   //!
      TBranch        *b_el_pt;   //!
      TBranch        *b_el_m;   //!
      TBranch        *b_el_eta;   //!
      TBranch        *b_el_phi;   //!
      TBranch        *b_el_status;   //!
      TBranch        *b_el_barcode;   //!
      TBranch        *b_el_charge;   //!
      TBranch        *b_el_parent_index;   //!
      TBranch        *b_mu_n;   //!
      TBranch        *b_mu_pt;   //!
      TBranch        *b_mu_m;   //!
      TBranch        *b_mu_eta;   //!
      TBranch        *b_mu_phi;   //!
      TBranch        *b_mu_status;   //!
      TBranch        *b_mu_barcode;   //!
      TBranch        *b_mu_charge;   //!
      TBranch        *b_mu_parent_index;   //!
      TBranch        *b_tau_n;   //!
      TBranch        *b_tau_pt;   //!
      TBranch        *b_tau_m;   //!
      TBranch        *b_tau_eta;   //!
      TBranch        *b_tau_phi;   //!
      TBranch        *b_tau_status;   //!
      TBranch        *b_tau_barcode;   //!
      TBranch        *b_tau_charge;   //!
      TBranch        *b_tau_parent_index;   //!
      TBranch        *b_tau_decay_index;   //!
      TBranch        *b_ph_n;   //!
      TBranch        *b_ph_pt;   //!
      TBranch        *b_ph_m;   //!
      TBranch        *b_ph_eta;   //!
      TBranch        *b_ph_phi;   //!
      TBranch        *b_ph_status;   //!
      TBranch        *b_ph_barcode;   //!

    private:

      friend class Particle;
      friend class Electron;
      friend class Muon;
      friend class Jet;
  };
}

#endif
