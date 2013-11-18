#ifndef OVERLAPREMOVAL_H
#define OVERLAPREMOVAL_H

#include <iostream>
#include <vector>

// =============================================================================
namespace TruthNtuple
{
  class Particle;
  class Lepton;
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
      template <class T> void removeObjects( std::vector<T*>&
                                           , const std::vector<bool>&
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

namespace TruthNtuple
{
  // -----------------------------------------------------------------------------
  template <class T>
    void OverlapRemoval::removeObjects( std::vector<T*>& objects
                                     , const std::vector<bool>& to_remove
                                     )
    {
      if (objects.size() != to_remove.size()) {
        std::cout << "ERROR! object vector size does not equal to remove size"
                  << "\n  objects size: " << objects.size()
                  << "\n  to_remove size: " << to_remove.size()
                  << "\n";

      }

      typename std::vector<T*>::iterator objects_begin = objects.begin();
      for (size_t it = 0; it != to_remove.size(); ++it) {
        size_t this_index = to_remove.size() - it - 1;
        if (to_remove.at(this_index)) {
          objects.erase(objects_begin + this_index);
        }
      }
    }
}

#endif
