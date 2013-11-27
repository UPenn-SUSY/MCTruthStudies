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
  // =============================================================================
  class Particle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Particle();
      Particle(const TruthNtuple::TruthNtupleLooper*, unsigned int mc_index);

      void setMCIndex(unsigned int);
      void setPdgid(int);
      void setPt(double);
      void setEta(double);
      void setPhi(double);
      void setE(double);
      void setPx(double);
      void setPy(double);
      void setPz(double);
      void setParentPdgid(int);
      void setParentMCIndex(int);
      void setParentBarcode(int);

      unsigned int getMCIndex() const;
      int getPdgid() const;
      double getPt() const;
      double getEta() const;
      double getPhi() const;
      double getE() const;
      double getPx() const;
      double getPy() const;
      double getPz() const;
      int getParentPdgid() const;
      int getParentMCIndex() const;
      int getParentBarcode() const;

      virtual void printGeneralInfo() const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      unsigned int m_mc_index;
      double       m_pdgid;
      double       m_pt;
      double       m_eta;
      double       m_phi;
      double       m_e;
      double       m_px;
      double       m_py;
      double       m_pz;
      int m_parent_pdgid;
      int m_parent_index;
      int m_parent_barcode;
  };

  // =============================================================================
  class Lepton : public Particle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Lepton();

      void setIsElectron(bool);
      void setCharge(double);

      bool isElectron() const;
      bool isMuon() const;
      double getCharge() const;

      virtual void print(TruthNtuple::TruthNtupleLooper* tnl = 0) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      bool m_is_electron;
      double m_charge;
  };

  // =============================================================================
  class Electron : public Lepton
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Electron();
      Electron(const TruthNtuple::TruthNtupleLooper*, unsigned int el_index);

      void setElIndex(unsigned int);

      unsigned int getElIndex() const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      unsigned int m_el_index;
  };

  // =============================================================================
  class Muon : public Lepton
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Muon();
      Muon(const TruthNtuple::TruthNtupleLooper*, unsigned int mu_index);

      void setMuIndex(unsigned int);

      unsigned int getMuIndex() const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      unsigned int m_mu_index;
  };

  // =============================================================================
  class Jet : public Particle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Jet();
      Jet(const TruthNtuple::TruthNtupleLooper*, unsigned int jet_index);

      void setJetIndex(unsigned int);
      void setTheta(double);
      void setIsBJet(bool);
      void setBQuarkIndex(int);

      unsigned int getJetIndex() const;
      double getTheta() const;
      bool getIsBJet() const;
      int getBQuarkIndex() const;

      virtual void print(TruthNtuple::TruthNtupleLooper* tnl = 0) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      unsigned int m_jet_index;
      double m_theta;
      bool m_is_b_jet;
      int m_b_quark_index;
  };

  // =============================================================================
  class Met
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Met();
      Met(double met_etx_noint, double met_ety_noint);
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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      double m_met_etx_noint;
      double m_met_ety_noint;
      double m_met_et_noint;
      double m_met_phi_noint;

      double m_met_rel_noint;
  };
}

#endif
