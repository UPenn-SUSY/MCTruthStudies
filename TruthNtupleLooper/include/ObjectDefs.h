#ifndef ObjectDefs_H
#define ObjectDefs_H

#include <vector>

// =============================================================================
namespace TruthNtuple
{
  class TruthNtupleLooper;
}

// =============================================================================
namespace TruthNtuple
{
  // ===========================================================================
  class Particle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Particle();
      Particle( const TruthNtuple::TruthNtupleLooper*
              , int mc_index
              , bool verbose = false
              );

      void setMCIndex(int);
      void setPdgid(int);
      void setStatus(int);
      void setBarcode(int);
      void setPt(double);
      void setEta(double);
      void setTheta(double);
      void setPhi(double);
      void setE(double);
      void setM(double);
      void setPx(double);
      void setPy(double);
      void setPz(double);
      void setParentPdgid(int);
      void setParentMCIndex(int);
      void setParentBarcode(int);
      void setImmediateParentPdgid(int);
      void setImmediateParentMCIndex(int);
      void setImmediateParentBarcode(int);

      int getMCIndex() const;
      int getPdgid() const;
      int getStatus() const;
      int getBarcode() const;
      double getPt() const;
      double getP() const;
      double getEta() const;
      double getTheta() const;
      double getPhi() const;
      double getE() const;
      double getM() const;
      double getPx() const;
      double getPy() const;
      double getPz() const;
      int getParentPdgid() const;
      int getParentMCIndex() const;
      int getParentBarcode() const;
      int getImmediateParentPdgid() const;
      int getImmediateParentMCIndex() const;
      int getImmediateParentBarcode() const;

      virtual void printGeneralInfo() const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      int          m_mc_index;
      double       m_pdgid;
      int          m_status;
      int          m_barcode;
      double       m_pt;
      double       m_eta;
      double       m_theta;
      double       m_phi;
      double       m_e;
      double       m_m;
      double       m_px;
      double       m_py;
      double       m_pz;
      int m_parent_pdgid;
      int m_parent_index;
      int m_parent_barcode;
      int m_immediate_parent_pdgid;
      int m_immediate_parent_index;
      int m_immediate_parent_barcode;
  };

  // ===========================================================================
  class Lepton : public Particle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Lepton();

      void setIsElectron(bool);
      void setCharge(double);

      bool isElectron() const;
      bool isMuon() const;
      double getCharge() const;

      virtual void print(TruthNtuple::TruthNtupleLooper* tnl = 0) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      bool m_is_electron;
      double m_charge;
  };

  // ===========================================================================
  class Electron : public Lepton
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Electron();
      Electron( const TruthNtuple::TruthNtupleLooper*
              , int el_index
              , bool get_final_state = true
              , bool verbose = false
              );

      void setElIndex(int);

      int getElIndex() const;

      virtual void print(TruthNtuple::TruthNtupleLooper* tnl = 0) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      int m_el_index;
  };

  // ===========================================================================
  class Muon : public Lepton
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Muon();
      Muon( const TruthNtuple::TruthNtupleLooper*
          , int mu_index
          , bool get_final_state = true
          , bool verbose = false
          );

      void setMuIndex(int);

      int getMuIndex() const;

      virtual void print(TruthNtuple::TruthNtupleLooper* tnl = 0) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      int m_mu_index;
  };

  // ===========================================================================
  class Jet : public Particle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Jet();
      Jet( const TruthNtuple::TruthNtupleLooper*
         , int jet_index
         , bool verbose = false
         );

      void setJetIndex(int);
      void setIsBJet(bool);
      void setBQuarkIndex(int);

      int getJetIndex() const;
      bool getIsBJet() const;
      int getBQuarkIndex() const;

      virtual void print(TruthNtuple::TruthNtupleLooper* tnl = 0) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      int m_jet_index;
      bool m_is_b_jet;
      int m_b_quark_index;
  };

  // ===========================================================================
  class Met
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Met();
      Met( double met_etx_noint
         , double met_ety_noint
         , bool verbose = false
         );
      void clear();

      void setMetNoint(double met_etx, double met_ety);
      void calculateMetRelNoint(const std::vector<TruthNtuple::Particle*>&);
      void calculateMetRelNoint( const std::vector<TruthNtuple::Electron*>&
                               , const std::vector<TruthNtuple::Muon*>&
                               , const std::vector<TruthNtuple::Jet*>&
                               );

      double getMetNoint() const;
      double getMetPhiNoint() const;
      double getMetRelNoint() const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      double m_met_etx_noint;
      double m_met_ety_noint;
      double m_met_et_noint;
      double m_met_phi_noint;

      double m_met_rel_noint;
  };
}

#endif
