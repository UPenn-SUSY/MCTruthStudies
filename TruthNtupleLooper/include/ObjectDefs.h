#ifndef ObjectDefs_H
#define ObjectDefs_H

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

      void setIndex(unsigned int);
      void setPt(double);
      void setEta(double);
      void setPhi(double);
      void setE(double);
      void setPx(double);
      void setPy(double);
      void setPz(double);
      void setParentPdgid(int);
      void setParentIndex(int);
      void setParentBarcode(int);

      unsigned int getIndex() const;
      double getPt() const;
      double getEta() const;
      double getPhi() const;
      double getE() const;
      double getPx() const;
      double getPy() const;
      double getPz() const;
      int getParentPdgid() const;
      int getParentIndex() const;
      int getParentBarcode() const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      unsigned int m_index;
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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
  };

  // =============================================================================
  class Muon : public Lepton
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Muon();
      Muon(const TruthNtuple::TruthNtupleLooper*, unsigned int mu_index);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
  };

  // =============================================================================
  class Jet : public Particle
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      Jet();
      Jet(const TruthNtuple::TruthNtupleLooper*, unsigned int jet_index);

      void setTheta(double);
      void setIsBJet(bool);

      double getTheta() const;
      bool getIsBJet() const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      double m_theta;
      bool m_is_b_jet;

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

      double getMetNoint() const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    protected:
      double m_met_etx_noint;
      double m_met_ety_noint;
      double m_met_et_noint;
  };
}

#endif
