#ifndef OVERLAPREMOVAL_H
#define OVERLAPREMOVAL_H

#include <vector>

// =============================================================================
namespace TruthNtuple
{
  class Electron;
  class Muon;
  class Jet;
}

// =============================================================================
namespace TruthNtuple
{
  // =============================================================================
  class OverlapRemoval
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      OverlapRemoval();

      void setDrEE(double);
      void setDrEJ(double);
      void setDrJE(double);
      void setDrJM(double);
      void setDrEM(double);
      void setDrMM(double);

      void doOverlapRemoval( const std::vector<TruthNtuple::Electron*>& ref_el
                           , const std::vector<TruthNtuple::Muon*>&     ref_mu
                           , const std::vector<TruthNtuple::Jet*>&      ref_jet
                           , std::vector<TruthNtuple::Electron*>&       final_el
                           , std::vector<TruthNtuple::Muon*>&           final_mu
                           , std::vector<TruthNtuple::Jet*>&            final_jet
                           );


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      double m_dr_ee;
      double m_dr_ej;
      double m_dr_je;
      double m_dr_jm;
      double m_dr_em;
      double m_dr_mm;
  };
}

#endif
