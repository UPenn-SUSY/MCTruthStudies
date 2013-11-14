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

      unsigned int getIndex();
      double getPt();
      double getEta();
      double getPhi();
      double getE();
      double getPx();
      double getPy();
      double getPz();

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
      Electron(const TruthNtuple::TruthNtupleLooper*, unsigned int el_index);

      void setCharge(double);
      void setParentPdgid(int);

      double getCharge();
      int getParentPdgid();

    protected:
      double m_charge;
      int m_parent_pdgid;
  };

  // =============================================================================
  class Muon : public Particle
  {
    public:
      Muon(const TruthNtuple::TruthNtupleLooper*, unsigned int mu_index);

      void setCharge(double);
      void setParentPdgid(int);

      double getCharge();
      int getParentPdgid();

    protected:
      double m_charge;
      int m_parent_pdgid;
  };

  // =============================================================================
  class Jet : public Particle
  {
    public:
      Jet(const TruthNtuple::TruthNtupleLooper*, unsigned int jet_index);

      void setTheta(double);

      double getTheta();

    protected:
      double m_theta;

  };
}
      
#endif
