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

      unsigned int getIndex() const;
      double getPt() const;
      double getEta() const;
      double getPhi() const;
      double getE() const;
      double getPx() const;
      double getPy() const;
      double getPz() const;

    protected:
      unsigned int m_index;
      double       m_pt;
      double       m_eta;
      double       m_phi;
      double       m_e;
      double       m_px;
      double       m_py;
      double       m_pz;
  };

  // =============================================================================
  class Electron : public Particle
  {
    public:
      Electron();
      Electron(const TruthNtuple::TruthNtupleLooper*, unsigned int el_index);

      void setCharge(double);
      void setParentPdgid(int);

      double getCharge() const;
      int getParentPdgid() const;

    protected:
      double m_charge;
      int m_parent_pdgid;
  };

  // =============================================================================
  class Muon : public Particle
  {
    public:
      Muon();
      Muon(const TruthNtuple::TruthNtupleLooper*, unsigned int mu_index);

      void setCharge(double);
      void setParentPdgid(int);

      double getCharge() const;
      int getParentPdgid() const;

    protected:
      double m_charge;
      int m_parent_pdgid;
  };

  // =============================================================================
  class Jet : public Particle
  {
    public:
      Jet();
      Jet(const TruthNtuple::TruthNtupleLooper*, unsigned int jet_index);

      void setTheta(double);
      void setIsBJet(bool);

      double getTheta() const;
      bool getIsBJet() const;

    protected:
      double m_theta;
      bool m_is_b_jet;

  };
}
      
#endif
