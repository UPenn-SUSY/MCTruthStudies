#include "include/OverlapRemoval.h"
#include "include/Calculators.h"
#include "include/ObjectDefs.h"

#include <vector>

// -----------------------------------------------------------------------------
TruthNtuple::OverlapRemoval::OverlapRemoval()
{
}

// -----------------------------------------------------------------------------
void TruthNtuple::OverlapRemoval::setDrEE(double val)
{
  m_dr_ee = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::OverlapRemoval::setDrEJ(double val)
{
  m_dr_ej = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::OverlapRemoval::setDrJE(double val)
{
  m_dr_je = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::OverlapRemoval::setDrJM(double val)
{
  m_dr_jm = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::OverlapRemoval::setDrEM(double val)
{
  m_dr_em = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::OverlapRemoval::setDrMM(double val)
{
  m_dr_mm = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::OverlapRemoval::doOverlapRemoval( const std::vector<TruthNtuple::Electron*>& ref_el
                                                  , const std::vector<TruthNtuple::Muon*>&     ref_mu
                                                  , const std::vector<TruthNtuple::Jet*>&      ref_jet
                                                  , std::vector<TruthNtuple::Electron*>&       final_el
                                                  , std::vector<TruthNtuple::Muon*>&           final_mu
                                                  , std::vector<TruthNtuple::Jet*>&            final_jet
                                                  )
{
  final_el.clear();
  final_mu.clear();
  final_jet.clear();

  final_el  = ref_el;
  final_mu  = ref_mu;
  final_jet = ref_jet;

  std::vector<int> el_to_remove;
  std::vector<int> mu_to_remove;
  std::vector<int> jet_to_remove;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // ee overlap removal
  el_to_remove.clear();
  el_to_remove.resize(final_el.size());
  // loop over electrons
  for (size_t el_it_1 = 0; el_it_1 != final_el.size(); ++el_it_1) {
    for (size_t el_it_2 = el_it_1+1; el_it_2 != final_el.size(); ++el_it_2) {
      // compute dr between two leptons
      // if dr < threshold, remove lower pt electron
      double dr = TruthNtuple::deltaR( final_el.at(el_it_1)
                                     , final_el.at(el_it_2)
                                     );
      if (dr < m_dr_ee) {
        if (final_el.at(el_it_1)->getPt() >= final_el.at(el_it_2)->getPt()) {
          el_to_remove.push_back(el_it_1);
        }
        else {
          el_to_remove.push_back(el_it_2);
        }
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // ej overlap removal
  jet_to_remove.clear();
  jet_to_remove.resize(final_jet.size());
  // loop over electrons jets
  for (size_t el_it = 0; el_it != final_el.size(); ++el_it) {
    for (size_t jet_it = 0; jet_it != final_jet.size(); ++jet_it) {
      // compute dr between two leptons
      // if dr < threshold, remove lower pt electron
      double dr = TruthNtuple::deltaR( final_el.at(el_it)
                                     , final_jet.at(jet_it)
                                     );
      if (dr < m_dr_ej) {
        jet_to_remove.push_back(jet_it);
        continue;
      }
    }
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // je overlap removal
  el_to_remove.clear();
  el_to_remove.resize(final_el.size());

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // jm overlap removal
  mu_to_remove.clear();
  mu_to_remove.resize(final_mu.size());

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // em overlap removal
  el_to_remove.clear();
  mu_to_remove.clear();
  el_to_remove.resize(final_el.size());
  mu_to_remove.resize(final_mu.size());

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // mm overlap removal
  mu_to_remove.clear();
  mu_to_remove.resize(final_mu.size());

}

