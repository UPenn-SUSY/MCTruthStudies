#include "include/OverlapRemoval.h"
#include "include/Calculators.h"
#include "include/ObjectDefs.h"

#include <iostream>
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // ee overlap removal
  {
    // define this vector in {}'s so it goes out of scope, and we can redefine it
    std::vector<bool> el_to_remove(final_el.size(), false);

    // loop over electrons
    for (size_t el_it_1 = 0; el_it_1 != final_el.size(); ++el_it_1) {
      for (size_t el_it_2 = el_it_1+1; el_it_2 != final_el.size(); ++el_it_2) {
        // compute dr between two objects
        // if dr < threshold, remove lower pt electron
        double dr = TruthNtuple::deltaR( final_el.at(el_it_1)
                                       , final_el.at(el_it_2)
                                       );
        if (dr < m_dr_ee) {
          if (final_el.at(el_it_1)->getPt() >= final_el.at(el_it_2)->getPt()) {
            el_to_remove.at(el_it_1) = true;
          }
          else {
            el_to_remove.at(el_it_2) = true;
          }
        }
      }
    }

    // remove electrons
    removeObjects(final_el, el_to_remove);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // ej overlap removal
  {
    // define this vector in {}'s so it goes out of scope, and we can redefine it
    std::vector<bool> jet_to_remove(final_jet.size(), false);

    // loop over electrons and jets
    for (size_t el_it = 0; el_it != final_el.size(); ++el_it) {
      for (size_t jet_it = 0; jet_it != final_jet.size(); ++jet_it) {
        // compute dr between two objects
        // if dr < threshold, remove jet
        double dr = TruthNtuple::deltaR( final_el.at(el_it)
                                       , final_jet.at(jet_it)
                                       );
        if (dr < m_dr_ej) {
          jet_to_remove.at(jet_it) = true;
          continue;
        }
      }
    }

    // remove jets
    removeObjects(final_jet, jet_to_remove);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // je overlap removal
  {
    // define this vector in {}'s so it goes out of scope, and we can redefine it
    std::vector<bool> el_to_remove(final_el.size(), false);

    // loop over jets and electrons
    for (size_t jet_it = 0; jet_it != final_jet.size(); ++jet_it) {
      for (size_t el_it = 0; el_it != final_el.size(); ++el_it) {
        // compute dr between two objects
        // if dr < threshold, remove electron
        double dr = TruthNtuple::deltaR( final_jet.at(jet_it)
                                      , final_el.at(el_it)
                                      );
        if (dr < m_dr_je) {
          el_to_remove.at(el_it) = true;
          continue;
        }
      }
    }

    // remove electrons
    removeObjects(final_el, el_to_remove);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // jm overlap removal
  {
    // define this vector in {}'s so it goes out of scope, and we can redefine it
    std::vector<bool> mu_to_remove(final_mu.size(), false);

    // loop over jets and electrons
    for (size_t jet_it = 0; jet_it != final_jet.size(); ++jet_it) {
      for (size_t mu_it = 0; mu_it != final_mu.size(); ++mu_it) {
        // compute dr between two objects
        // if dr < threshold, remove muon
        double dr = TruthNtuple::deltaR( final_jet.at(jet_it)
                                       , final_mu.at(mu_it)
                                       );
        if (dr < m_dr_jm) {
          mu_to_remove.at(mu_it) = true;
          continue;
        }
      }
    }

    // remove muons
    removeObjects(final_mu, mu_to_remove);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // em overlap removal
  {
    // define these vectors in {}'s so they go out of scope, and we can redefine them
    std::vector<bool> el_to_remove(final_el.size(), false);
    std::vector<bool> mu_to_remove(final_mu.size(), false);

    // loop over electrons and muons
    for (size_t el_it = 0; el_it != final_el.size(); ++el_it) {
      for (size_t mu_it = 0; mu_it != final_mu.size(); ++mu_it) {
        // compute dr between two objects
        // if dr < threshold, remove both
        double dr = TruthNtuple::deltaR( final_el.at(el_it)
                                       , final_mu.at(mu_it)
                                       );
        if (dr < m_dr_em) {
          el_to_remove.at(el_it) = true;
          mu_to_remove.at(mu_it) = true;
        }
      }
    }

    // remove electrons and muons
    removeObjects(final_el, el_to_remove);
    removeObjects(final_mu, mu_to_remove);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // mm overlap removal
  {
    // define this vector in {}'s so it goes out of scope, and we can redefine it
    std::vector<bool> mu_to_remove(final_mu.size(), false);

    // loop over muons
    for (size_t mu_it_1 = 0; mu_it_1 != final_el.size(); ++mu_it_1) {
      for (size_t mu_it_2 = mu_it_1+1; mu_it_2 != final_el.size(); ++mu_it_2) {
        // compute dr between two objects
        // if dr < threshold, remove both muons
        double dr = TruthNtuple::deltaR( final_el.at(mu_it_1)
                                       , final_el.at(mu_it_2)
                                       );
        if (dr < m_dr_mm) {
            mu_to_remove.at(mu_it_1) = true;
            mu_to_remove.at(mu_it_2) = true;
        }
      }
    }

    // remove muons
    removeObjects(final_mu, mu_to_remove);
  }
}

